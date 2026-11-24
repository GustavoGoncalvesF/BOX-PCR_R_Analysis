#####################################################################################
# Project: BOX-PCR Genomic Fingerprinting
# Script: 02_band_matching.R
# Author: Gustavo Goncalves
#
# Description:
# Match BOX-PCR bands across isolates and gels using molecular sizes,
# optional reference-based alignment and relative size tolerances.
#
# Input:
#   bands_isolates.csv
#   bands_reference.csv
#   isolates_metadata.csv
#
# Output:
#   Reference alignment QC
#   Band matching assignments
#   Binary presence/absence matrices
#   Matching sensitivity summary
#
# Notes:
#   - Molecular sizes were estimated using Cubic Spline in Image Lab.
#   - Bands below 500 bp are extrapolated beyond the DNA ladder range.
#   - Reference-based alignment is exploratory.
#   - Band classes are defined by complete-linkage clustering.
#####################################################################################

# Set up packages
library(dplyr)
library(ggplot2)

# Project directories
project_path <- getwd()
tables_path <- file.path(project_path, "results", "tables")
figures_path <- file.path(project_path, "results", "figures")

dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
tolerances <- c(0.01, 0.02, 0.05)
primary_tolerance <- 0.02

apply_reference_correction <- TRUE
ladder_min_bp <- 500
ladder_max_bp <- 10000

# Import data (Brazilian CSV format)
bands <- read.csv2(
  file.path(tables_path, "bands_isolates.csv"),
  stringsAsFactors = FALSE
)

reference <- read.csv2(
  file.path(tables_path, "bands_reference.csv"),
  stringsAsFactors = FALSE
)

metadata <- read.csv2(
  file.path(tables_path, "isolates_metadata.csv"),
  stringsAsFactors = FALSE
)

# Initial validation
stopifnot(
  all(c("gel", "lane", "isolate", "band_number",
        "size_bp", "relative_front") %in% names(bands)),
  all(c("gel", "band_number", "size_bp") %in% names(reference)),
  all(c("isolate", "gel") %in% names(metadata)),
  all(is.finite(bands$size_bp)),
  all(bands$size_bp > 0),
  all(is.finite(reference$size_bp)),
  all(reference$size_bp > 0),
  !anyDuplicated(metadata$isolate),
  all(bands$isolate %in% metadata$isolate),
  setequal(unique(bands$gel), unique(reference$gel)),
  primary_tolerance %in% tolerances
)

# Reference alignment
# Each gel must contain five reference bands
reference <- reference %>%
  arrange(gel, band_number)

ref_counts <- reference %>%
  count(gel)

if (any(ref_counts$n != 5L)) {
  stop("Expected exactly five reference bands per gel.")
}

if (any(duplicated(reference[c("gel", "band_number")]))) {
  stop("Duplicated reference band numbers detected.")
}

if (!setequal(
  unique(reference$band_number),
  1:5
)) {
  stop("Unexpected I7-REF band numbering.")
}

# Geometric mean size of each reference band
reference_target <- reference %>%
  group_by(band_number) %>%
  summarise(
    target_bp = exp(mean(log(size_bp))),
    .groups = "drop"
  )

reference <- reference %>%
  left_join(reference_target, by = "band_number") %>%
  mutate(
    deviation_percent = 100 * (
      size_bp / target_bp - 1
    ),
    log_correction = log(target_bp) - log(size_bp)
  )

# Apply size correction using reference bands
# Correction is linearly interpolated in log-size space.
# Outside the reference range, the nearest correction is used.
correct_sizes <- function(sizes, ref_data) {

  ref_data <- ref_data %>%
    arrange(size_bp)

  correction <- approx(
    x = log(ref_data$size_bp),
    y = ref_data$log_correction,
    xout = log(sizes),
    method = "linear",
    rule = 2
  )$y

  exp(log(sizes) + correction)
}

# Preserve original estimates
bands <- bands %>%
  mutate(
    size_original_bp = size_bp,
    size_corrected_bp = size_bp,
    calibration_region = case_when(
      size_original_bp < ladder_min_bp ~ "Below ladder",
      size_original_bp > ladder_max_bp ~ "Above ladder",
      TRUE ~ "Within ladder"
    )
  )

# Optional inter-gel alignment
if (apply_reference_correction) {

  for (g in unique(bands$gel)) {

    idx <- bands$gel == g
    ref_gel <- reference %>%
      filter(gel == g)

    bands$size_corrected_bp[idx] <- correct_sizes(
      bands$size_original_bp[idx],
      ref_gel
    )
  }
}

# Reference QC export
write.csv2(
  reference,
  file.path(tables_path, "reference_alignment.csv"),
  row.names = FALSE
)

# Reference QC plot
p_ref <- ggplot(
  reference,
  aes(
    x = factor(band_number),
    y = deviation_percent,
    color = gel,
    group = gel
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_line() +
  geom_point(size = 2.5) +
  labs(
    title = "BOX-PCR - Reference band deviations",
    x = "Reference band",
    y = "Deviation from geometric mean (%)",
    color = "Gel"
  ) +
  theme_bw()

print(p_ref)

ggsave(
  file.path(figures_path, "reference_alignment.png"),
  plot = p_ref,
  width = 8,
  height = 5,
  dpi = 300
)

# Band matching using complete-linkage clustering
match_bands <- function(data, tolerance) {

  sizes <- data$size_corrected_bp

  # Absolute differences in log-size space
  distance <- dist(
    matrix(log(sizes), ncol = 1),
    method = "euclidean"
  )

  clustering <- hclust(
    distance,
    method = "complete"
  )

  groups <- cutree(
    clustering,
    h = log1p(tolerance)
  )

  result <- data %>%
    mutate(cluster = groups)

  # Order classes by decreasing median size
  class_info <- result %>%
    group_by(cluster) %>%
    summarise(
      median_bp = median(size_corrected_bp),
      min_bp = min(size_corrected_bp),
      max_bp = max(size_corrected_bp),
      n_bands = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(median_bp)) %>%
    mutate(
      band_class = sprintf(
        "B%03d",
        row_number()
      )
    )

  result <- result %>%
    left_join(
      class_info %>%
        select(cluster, band_class),
      by = "cluster"
    ) %>%
    select(-cluster)

  list(
    assignments = result,
    classes = class_info
  )
}

# Build binary matrix
create_binary_matrix <- function(assignments, metadata) {

  isolates <- metadata$isolate
  classes <- unique(assignments$band_class)

  matrix_binary <- matrix(
    0L,
    nrow = length(isolates),
    ncol = length(classes),
    dimnames = list(isolates, classes)
  )

  for (i in seq_len(nrow(assignments))) {

    matrix_binary[
      assignments$isolate[i],
      assignments$band_class[i]
    ] <- 1L
  }

  matrix_binary
}

# Run sensitivity analysis
summary_results <- list()

for (tol in tolerances) {

  tol_label <- paste0(
    round(tol * 100),
    "pct"
  )

  matching <- match_bands(bands, tol)

  assignments <- matching$assignments

  # Check whether multiple detected bands from the same
  # isolate were assigned to one class
  conflicts <- assignments %>%
    count(
      isolate, band_class,
      name = "n_bands"
    ) %>%
    filter(n_bands > 1)

  if (nrow(conflicts) > 0) {
    warning(
      tol_label, ": ",
      nrow(conflicts),
      " isolate/class combinations contain multiple bands."
    )
  }

  # Export matching assignments
  write.csv2(
    assignments,
    file.path(
      tables_path,
      paste0("band_assignments_", tol_label, ".csv")
    ),
    row.names = FALSE
  )

  # Export class summary
  write.csv2(
    matching$classes,
    file.path(
      tables_path,
      paste0("band_classes_", tol_label, ".csv")
    ),
    row.names = FALSE
  )

  # Export conflicts for manual evaluation
  write.csv2(
    conflicts,
    file.path(
      tables_path,
      paste0("matching_conflicts_", tol_label, ".csv")
    ),
    row.names = FALSE
  )

  # Matrix with all detected bands
  binary_all <- create_binary_matrix(
    assignments,
    metadata
  )

  write.csv2(
    data.frame(
      isolate = rownames(binary_all),
      binary_all,
      check.names = FALSE
    ),
    file.path(
      tables_path,
      paste0("binary_matrix_all_", tol_label, ".csv")
    ),
    row.names = FALSE
  )

  # Matrix excluding bands below 500 bp
  assignments_filtered <- assignments %>%
    filter(calibration_region == "Within ladder")

  binary_filtered <- create_binary_matrix(
    assignments_filtered,
    metadata
  )

  # Remove classes without retained bands
  if (ncol(binary_filtered) > 0) {
    binary_filtered <- binary_filtered[
      ,
      colSums(binary_filtered) > 0,
      drop = FALSE
    ]
  }

  write.csv2(
    data.frame(
      isolate = rownames(binary_filtered),
      binary_filtered,
      check.names = FALSE
    ),
    file.path(
      tables_path,
      paste0("binary_matrix_calibrated_", tol_label, ".csv")
    ),
    row.names = FALSE
  )

  # Store sensitivity summary
  summary_results[[tol_label]] <- data.frame(
    tolerance_percent = tol * 100,
    isolates = nrow(binary_all),
    classes_all = ncol(binary_all),
    classes_calibrated = ncol(binary_filtered),
    matching_conflicts = nrow(conflicts),
    isolates_without_calibrated_bands = sum(
      rowSums(binary_filtered) == 0
    )
  )
}

# Consolidate sensitivity results
matching_summary <- bind_rows(summary_results)

write.csv2(
  matching_summary,
  file.path(tables_path, "matching_summary.csv"),
  row.names = FALSE
)

# Print final summary
cat("\nBOX-PCR BAND MATCHING SUMMARY\n")
cat("--------------------------------\n")
cat("Isolates:", nrow(metadata), "\n")
cat("Experimental bands:", nrow(bands), "\n")
cat("Reference alignment:",
    apply_reference_correction, "\n")
cat("Primary tolerance:",
    primary_tolerance * 100, "%\n\n")

print(matching_summary)
