################################################################################
# BOX-PCR - Group characterization and PGP isolate selection
# Script: 06_isolate_selection.R
# Input: UPGMA, Dice, binary matrix, and isolate metadata
# Output: labeled group barplot, group summary, provisional selection of 30
################################################################################

library(ggplot2)
library(readxl)

# Configuration
project_path <- getwd()
tables_path <- file.path(project_path, "results", "tables")
output_path <- file.path(project_path, "results", "isolate_selection")
metadata_file <- file.path(project_path, "data", "metadata", "metadata-box.xlsx")
dissimilarity_cut <- 0.50     # 50% dissimilarity = 50% Dice similarity
n_target <- 30L
min_informative_bands <- 2L  # Flag only; no isolate is automatically discarded

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Import
hc <- readRDS(file.path(tables_path, "upgma.rds"))
dice <- as.matrix(read.csv2(file.path(tables_path, "dice_similarity.csv"),
                            row.names = 1, check.names = FALSE))
binary <- as.matrix(read.csv2(file.path(tables_path, "binary_matrix.csv"),
                              row.names = 1, check.names = FALSE))
metadata <- as.data.frame(read_excel(metadata_file))

required <- c("Isolate", "Morphotype", "BAG", "Species", "Accession")
if (!all(required %in% names(metadata))) {
  stop("Missing metadata columns: ",
       paste(setdiff(required, names(metadata)), collapse = ", "))
}
if (anyDuplicated(metadata$Isolate)) stop("Duplicated Isolate entries in metadata.")
if (anyDuplicated(rownames(dice)) ||
    !setequal(rownames(dice), colnames(dice)) ||
    !setequal(rownames(dice), rownames(binary)) ||
    !setequal(rownames(dice), hc$labels)) {
  stop("UPGMA, Dice and binary matrix have incompatible isolate identifiers.")
}

ids <- hc$labels
if (!all(ids %in% metadata$Isolate)) {
  stop("Isolates missing in metadata: ",
       paste(setdiff(ids, metadata$Isolate), collapse = ", "))
}
dice <- dice[ids, ids, drop = FALSE]
binary <- binary[ids, , drop = FALSE]
metadata <- metadata[match(ids, metadata$Isolate), required, drop = FALSE]
if (anyNA(metadata)) stop("Metadata contain missing values for included isolates.")
if (any(!is.finite(dice)) || any(dice < -1e-8 | dice > 1 + 1e-8) ||
    max(abs(dice - t(dice))) > 1e-8) stop("Invalid Dice similarity matrix.")

metadata$Morphotype <- as.character(metadata$Morphotype)
metadata$BAG <- as.character(metadata$BAG)
metadata$Species <- as.character(metadata$Species)
metadata$Accession <- as.character(metadata$Accession)

# Cut UPGMA and count bands
cluster <- cutree(hc, h = dissimilarity_cut)
cluster <- as.integer(cluster[ids])
number <- as.integer(sub("^I", "", ids))
members <- data.frame(
  Isolate = ids,
  Group = cluster,
  n_bands = rowSums(binary),
  metadata[, setdiff(required, "Isolate"), drop = FALSE],
  stringsAsFactors = FALSE
)
members <- members[order(members$Group, number), , drop = FALSE]
groups <- sort(unique(members$Group))

# Choose one representative per group: medoid by average within-group Dice.
# When possible, avoid representatives with fewer than 2 detected bands.
choose_representative <- function(group_ids) {
  candidates <- group_ids
  sufficient <- group_ids[rowSums(binary[group_ids, , drop = FALSE]) >=
                            min_informative_bands]
  if (length(sufficient)) candidates <- sufficient
  mean_sim <- vapply(candidates, function(x) {
    mean(dice[x, group_ids])
  }, numeric(1))
  band_n <- rowSums(binary[candidates, , drop = FALSE])
  numeric_id <- as.integer(sub("^I", "", candidates))
  candidates[order(-mean_sim, -band_n, numeric_id)][1]
}

summary_rows <- lapply(groups, function(g) {
  sub <- members[members$Group == g, , drop = FALSE]
  group_ids <- sub$Isolate
  rep_id <- choose_representative(group_ids)
  rep_row <- members[match(rep_id, members$Isolate), , drop = FALSE]
  pair_sim <- if (length(group_ids) > 1) {
    z <- dice[group_ids, group_ids, drop = FALSE]
    z[upper.tri(z)]
  } else NA_real_
  data.frame(
    Group = g,
    n_isolates = length(group_ids),
    Members = paste(group_ids, collapse = ", "),
    Representative = rep_id,
    Rep_n_bands = rep_row$n_bands,
    Mean_within_group_Dice = mean(pair_sim),
    Rep_mean_group_Dice = mean(dice[rep_id, group_ids]),
    Low_resolution_group = all(sub$n_bands < min_informative_bands),
    BAG = rep_row$BAG,
    Species = rep_row$Species,
    Accession = rep_row$Accession,
    Morphotype = rep_row$Morphotype,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL

# Provisional selection: start with one representative from every group and
# remove the 11 least costly for collection-wide Dice coverage. A removal is
# disallowed if it would eliminate the last representative of any metadata level.
# This is a sampling aid, NOT a biological validation of the 30 isolates.
selected <- summary$Representative
coverage <- function(representatives) {
  similarities <- dice[ids, representatives, drop = FALSE]
  mean(apply(similarities, 1, max))
}
keeps_categories <- function(representatives) {
  chosen <- metadata[match(representatives, metadata$Isolate), , drop = FALSE]
  for (v in c("BAG", "Species", "Accession", "Morphotype")) {
    if (!all(unique(metadata[[v]]) %in% unique(chosen[[v]]))) return(FALSE)
  }
  TRUE
}

if (n_target < 1L) stop("n_target must be >= 1.")
if (n_target > length(selected)) {
  warning("n_target exceeds the number of groups; retaining all groups.")
}
while (length(selected) > n_target) {
  options <- selected[vapply(selected, function(r) {
    keeps_categories(setdiff(selected, r))
  }, logical(1))]
  if (!length(options)) {
    warning("Cannot reach n_target without losing a metadata category.")
    break
  }
  scores <- vapply(options, function(r) coverage(setdiff(selected, r)), numeric(1))
  candidate_rows <- summary[match(options, summary$Representative), , drop = FALSE]
  # Tie-break: remove less informative representatives first, then smaller groups.
  ord <- order(-scores, candidate_rows$Rep_n_bands,
               candidate_rows$n_isolates, candidate_rows$Group)
  selected <- setdiff(selected, options[ord[1]])
}

summary$Selection <- ifelse(summary$Representative %in% selected,
                            "Selected (provisional)", "Reserve")
summary$Nearest_selected <- vapply(summary$Representative, function(rep_id) {
  candidates <- selected
  candidates[which.max(dice[rep_id, candidates])]
}, character(1))
summary$Dice_to_nearest_selected <- vapply(seq_len(nrow(summary)), function(i) {
  dice[summary$Representative[i], summary$Nearest_selected[i]]
}, numeric(1))
summary$Review <- ifelse(summary$Low_resolution_group,
                         "Low resolution: all fingerprints have 1 band",
                         ifelse(summary$Rep_n_bands <= 2,
                                "Review limited number of bands", ""))

members$Selection <- summary$Selection[match(members$Group, summary$Group)]
members$Representative <- summary$Representative[match(members$Group, summary$Group)]
members$Role <- ifelse(members$Isolate == members$Representative,
                       "Representative", "Other member")

# Export tables
write.csv2(members, file.path(output_path, "group_membership.csv"), row.names = FALSE)
write.csv2(summary, file.path(output_path, "group_summary.csv"), row.names = FALSE)
write.csv2(summary[summary$Selection == "Selected (provisional)", , drop = FALSE],
           file.path(output_path, "selected_30_provisional.csv"), row.names = FALSE)
write.csv2(summary[summary$Selection == "Reserve", , drop = FALSE],
           file.path(output_path, "reserve_groups.csv"), row.names = FALSE)

# Labeled stacked bar chart: one cell = one isolate
members$Group_label <- factor(paste0("G", members$Group),
                              levels = paste0("G", groups))
members$Position <- ave(seq_len(nrow(members)), members$Group, FUN = seq_along)
max_members <- max(summary$n_isolates)
p <- ggplot(members, aes(x = Group_label, y = Position - 0.5)) +
  geom_tile(aes(fill = Selection), width = 0.88, height = 0.91,
            color = "white", size = 0.25) +
  geom_text(aes(label = Isolate), size = 2.4) +
  scale_fill_manual(values = c("Selected (provisional)" = "#86B5A7",
                               "Reserve" = "#D9DDE3"), drop = FALSE) +
  scale_y_continuous(breaks = seq(0, max_members, 1),
                     limits = c(0, max_members + 0.6), expand = c(0, 0)) +
  labs(title = "BOX-PCR - UPGMA group composition",
       subtitle = paste0("Dice dissimilarity cut: ",
                         format(dissimilarity_cut, decimal.mark = ","),
                         " | ", length(groups), " groups | ",
                         length(selected), " preliminary representatives"),
       x = "UPGMA group", y = "Number of isolates", fill = "PGP selection") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "bottom")
ggsave(file.path(output_path, "group_composition_labeled.png"), p,
       width = max(15, 0.45 * length(groups)), height = 8.5, dpi = 300,
       limitsize = FALSE)
ggsave(file.path(output_path, "group_composition_labeled.pdf"), p,
       width = max(15, 0.45 * length(groups)), height = 8.5,
       limitsize = FALSE)

cat("\nBOX-PCR | isolate selection\n")
cat("Isolates:", length(ids), "| Groups:", length(groups), "\n")
cat("Suggested representatives:", length(selected), "\n")
cat("Low-resolution groups:", sum(summary$Low_resolution_group), "\n")
cat("Mean Dice similarity to closest selected representative:",
    round(coverage(selected), 3), "\n")
cat("Saved to:", output_path, "\n")
