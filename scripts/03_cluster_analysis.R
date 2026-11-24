#####################################################################################
# Project: BOX-PCR Genomic Fingerprinting
# Script: 03_cluster_analysis.R
# Author: Gustavo Goncalves
#
# Description:
# Calculate Dice similarity and perform UPGMA hierarchical clustering
# of BOX-PCR fingerprints using binary band presence/absence matrices.
#
# Input:
#   binary_matrix_all_*.csv
#   binary_matrix_calibrated_*.csv
#
# Output:
#   Dice similarity and distance matrices
#   UPGMA dendrograms
#   Similarity heatmaps
#   Clustering sensitivity and comparison tables
#
# Notes:
#   - Dice similarity is calculated directly from binary matrices.
#   - UPGMA clustering uses average linkage.
#   - Isolates without retained bands are excluded per analysis.
#   - No similarity cutoff is imposed to define isolate groups.
#####################################################################################

# Set up packages
library(ggplot2)

# Project directories
project_path <- getwd()
tables_path <- file.path(project_path, "results", "tables")
figures_path <- file.path(project_path, "results", "figures")

cluster_tables_path <- file.path(tables_path, "clustering")
cluster_figures_path <- file.path(figures_path, "clustering")

dir.create(cluster_tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_figures_path, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
tolerances <- c(1, 2, 5)
matrix_types <- c("all", "calibrated")

primary_tolerance <- 5
primary_matrix <- "all"

# Import binary matrix
read_binary_matrix <- function(file) {

  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  data <- read.csv2(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(data) < 2) {
    stop("No band classes found in: ", file)
  }

  if (anyNA(data) || anyDuplicated(data[[1]])) {
    stop("Missing values or duplicated isolates: ", file)
  }

  isolates <- as.character(data[[1]])

  if (anyNA(isolates) || any(!nzchar(isolates))) {
    stop("Invalid isolate identifiers: ", file)
  }

  # Convert band columns to numeric matrix
  binary <- as.matrix(data[, -1, drop = FALSE])
  suppressWarnings(storage.mode(binary) <- "numeric")

  rownames(binary) <- isolates

  if (anyNA(binary) || !all(binary %in% c(0, 1))) {
    stop("Matrix must contain only 0 and 1: ", file)
  }

  if (anyDuplicated(colnames(binary))) {
    stop("Duplicated band class names: ", file)
  }

  # Sort isolates numerically
  isolate_number <- suppressWarnings(
    as.integer(sub("^I", "", isolates))
  )

  if (anyNA(isolate_number)) {
    stop("Unexpected isolate identifiers: ", file)
  }

  binary <- binary[
    order(isolate_number),
    ,
    drop = FALSE
  ]

  return(binary)
}

# Calculate Dice similarity
dice_similarity <- function(binary) {

  if (nrow(binary) < 2) {
    stop("At least two isolates are required.")
  }

  if (any(rowSums(binary) == 0)) {
    stop("Zero-band isolates must be removed before Dice.")
  }

  # Shared bands between each isolate pair
  shared <- tcrossprod(binary)

  # Number of bands in each isolate
  n_bands <- rowSums(binary)

  # Dice similarity: 2a / (2a + b + c)
  denominator <- outer(n_bands, n_bands, "+")

  similarity <- 2 * shared / denominator

  diag(similarity) <- 1

  dimnames(similarity) <- list(
    rownames(binary),
    rownames(binary)
  )

  return(similarity)
}

# Export square matrix
export_matrix <- function(mat, filename) {

  output <- data.frame(
    isolate = rownames(mat),
    mat,
    check.names = FALSE
  )

  write.csv2(
    output,
    file.path(cluster_tables_path, filename),
    row.names = FALSE
  )
}

# Create UPGMA dendrogram
plot_dendrogram <- function(hc, title, filename) {

  png(
    filename = file.path(cluster_figures_path, filename),
    width = 3600,
    height = 2000,
    res = 300
  )

  par(mar = c(7, 5, 4, 2))

  plot(
    hc,
    main = title,
    xlab = "",
    sub = "",
    ylab = "Dice dissimilarity (1 - similarity)",
    cex = 0.65,
    hang = -1
  )

  dev.off()
}

# Create Dice similarity heatmap
plot_similarity_heatmap <- function(similarity, hc,
                                    title, filename) {

  isolate_order <- hc$labels[hc$order]

  ordered <- similarity[
    isolate_order,
    isolate_order,
    drop = FALSE
  ]

  heatmap_data <- as.data.frame(
    as.table(ordered),
    stringsAsFactors = FALSE
  )

  names(heatmap_data) <- c(
    "Isolate1",
    "Isolate2",
    "Similarity"
  )

  heatmap_data$Isolate1 <- factor(
    heatmap_data$Isolate1,
    levels = isolate_order
  )

  heatmap_data$Isolate2 <- factor(
    heatmap_data$Isolate2,
    levels = rev(isolate_order)
  )

  p <- ggplot(
    heatmap_data,
    aes(
      x = Isolate1,
      y = Isolate2,
      fill = Similarity
    )
  ) +
    geom_tile() +
    scale_fill_gradient(
      low = "white",
      high = "steelblue4",
      limits = c(0, 1)
    ) +
    coord_equal() +
    labs(
      title = title,
      x = "Isolate",
      y = "Isolate",
      fill = "Dice"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6
      ),
      axis.text.y = element_text(size = 6)
    )

  print(p)

  ggsave(
    filename = file.path(
      cluster_figures_path, filename
    ),
    plot = p,
    width = 11,
    height = 10,
    dpi = 300
  )
}

# Run clustering analyses
results <- list()
summary_results <- list()

for (tol in tolerances) {

  for (type in matrix_types) {

    analysis_id <- paste0(
      type, "_", tol, "pct"
    )

    cat("\nRunning analysis:", analysis_id, "\n")

    # Import binary matrix
    input_file <- file.path(
      tables_path,
      paste0("binary_matrix_", analysis_id, ".csv")
    )

    binary <- read_binary_matrix(input_file)

    # Exclude isolates without retained bands
    excluded <- rownames(binary)[
      rowSums(binary) == 0
    ]

    binary <- binary[
      rowSums(binary) > 0,
      ,
      drop = FALSE
    ]

    if (nrow(binary) < 2) {
      stop(
        analysis_id,
        ": fewer than two isolates with bands."
      )
    }

    # Remove unused band classes
    binary <- binary[
      ,
      colSums(binary) > 0,
      drop = FALSE
    ]

    # Dice similarity and dissimilarity
    similarity <- dice_similarity(binary)
    dissimilarity <- 1 - similarity

    # Hierarchical clustering (UPGMA)
    hc <- hclust(
      as.dist(dissimilarity),
      method = "average"
    )

    # Export similarity and distance matrices
    export_matrix(
      similarity,
      paste0("dice_similarity_", analysis_id, ".csv")
    )

    export_matrix(
      dissimilarity,
      paste0("dice_distance_", analysis_id, ".csv")
    )

    # Export dendrogram object
    saveRDS(
      hc,
      file.path(
        cluster_tables_path,
        paste0("upgma_", analysis_id, ".rds")
      )
    )

    # Export order of isolates in dendrogram
    dendrogram_order <- data.frame(
      order = seq_along(hc$order),
      isolate = hc$labels[hc$order]
    )

    write.csv2(
      dendrogram_order,
      file.path(
        cluster_tables_path,
        paste0("dendrogram_order_", analysis_id, ".csv")
      ),
      row.names = FALSE
    )

    # Dendrogram
    plot_dendrogram(
      hc,
      title = paste0(
        "BOX-PCR - UPGMA (Dice) - ",
        tol, "% - ", type
      ),
      filename = paste0(
        "upgma_", analysis_id, ".png"
      )
    )

    # Similarity heatmap
    plot_similarity_heatmap(
      similarity,
      hc,
      title = paste0(
        "BOX-PCR - Dice Similarity - ",
        tol, "% - ", type
      ),
      filename = paste0(
        "dice_heatmap_", analysis_id, ".png"
      )
    )

    # Summary statistics from unique isolate pairs
    pairwise_similarity <- similarity[
      upper.tri(similarity)
    ]

    summary_results[[analysis_id]] <- data.frame(
      analysis = analysis_id,
      tolerance_percent = tol,
      matrix_type = type,
      isolates = nrow(binary),
      band_classes = ncol(binary),
      excluded_isolates = length(excluded),
      mean_bands_per_isolate = mean(rowSums(binary)),
      mean_dice = mean(pairwise_similarity),
      median_dice = median(pairwise_similarity),
      min_dice = min(pairwise_similarity),
      max_dice = max(pairwise_similarity),
      identical_profiles = sum(
        pairwise_similarity == 1
      )
    )

    # Store results for comparisons
    results[[analysis_id]] <- list(
      binary = binary,
      similarity = similarity,
      dissimilarity = dissimilarity,
      hc = hc,
      excluded = excluded
    )
  }
}

# Clustering summary
clustering_summary <- do.call(
  rbind,
  summary_results
)

rownames(clustering_summary) <- NULL

write.csv2(
  clustering_summary,
  file.path(
    cluster_tables_path,
    "clustering_summary.csv"
  ),
  row.names = FALSE
)

# Compare dendrograms between tolerances
# Cophenetic distance correlation (Spearman).
# Each comparison uses only isolates shared by both analyses.
comparison_results <- list()

for (type in matrix_types) {

  ids <- paste0(type, "_", tolerances, "pct")

  combinations <- combn(ids, 2, simplify = FALSE)

  for (pair in combinations) {

    result1 <- results[[pair[1]]]
    result2 <- results[[pair[2]]]

    common <- intersect(
      rownames(result1$binary),
      rownames(result2$binary)
    )

    if (length(common) < 3) {
      next
    }

    # Cophenetic distances for shared isolates
    coph1 <- as.matrix(
      cophenetic(result1$hc)
    )[common, common, drop = FALSE]

    coph2 <- as.matrix(
      cophenetic(result2$hc)
    )[common, common, drop = FALSE]

    x <- coph1[upper.tri(coph1)]
    y <- coph2[upper.tri(coph2)]

    correlation <- if (
      length(unique(x)) > 1 &&
      length(unique(y)) > 1
    ) {
      cor(x, y, method = "spearman")
    } else {
      NA_real_
    }

    comparison_results[[length(comparison_results) + 1L]] <-
      data.frame(
        matrix_type = type,
        analysis_1 = pair[1],
        analysis_2 = pair[2],
        common_isolates = length(common),
        cophenetic_spearman = correlation
      )
  }
}

# Export clustering comparisons
clustering_comparison <- do.call(
  rbind,
  comparison_results
)

write.csv2(
  clustering_comparison,
  file.path(
    cluster_tables_path,
    "clustering_comparison.csv"
  ),
  row.names = FALSE
)

# Final summary
cat("\nBOX-PCR CLUSTER ANALYSIS SUMMARY\n")
cat("--------------------------------\n")
cat("Primary analysis:",
    paste0(primary_matrix, "_",
           primary_tolerance, "pct"), "\n\n")

print(clustering_summary)

cat("\nDENDROGRAM COMPARISON\n")
print(clustering_comparison)
