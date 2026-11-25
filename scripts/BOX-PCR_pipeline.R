#####################################################################################
# Project: BOX-PCR Genomic Fingerprinting
# Script: BOX-PCR_pipeline.R
# Author: Gustavo Goncalves
#
# Description:
# Complete pipeline: Image Lab Excel import, I7-REF alignment, band matching
# (5% size tolerance), Dice similarity, and UPGMA clustering.
#
# Input:  data/raw/Gel*-Analysis-Table.xlsx
# Output: results/tables/ and results/figures/
#
# Lane layout (all gels):
#   1 = ladder; 2 = I7-REF; 3-18 = experimental isolates; 19 = ladder
#   GEL1 = I1-I16; GEL2 = I17-I32; ...; GEL6 = I81-I92
#   Experimental lanes without detected bands are omitted.
#
# Run from the RStudio project root: source("scripts/BOX-PCR_pipeline.R")
#####################################################################################

# Packages
library(readxl)
library(dplyr)
library(ggplot2)

# Paths: open BOX-PCR_R_Analysis.Rproj before running
project_path <- getwd()
raw_path <- file.path(project_path, "data", "raw")
tables_path <- file.path(project_path, "results", "tables")
figures_path <- file.path(project_path, "results", "figures")
dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

# Fixed parameters
n_isolates_total <- 92L
n_lanes <- 19L
columns_per_lane <- 11L
isolates_per_gel <- 16L
tolerance <- 0.05

# Find and order all Excel exports
files <- list.files(raw_path, pattern = "\\.xlsx$", full.names = TRUE,
                    ignore.case = TRUE)
files <- files[!grepl("^~\\$", basename(files))]
if (length(files) == 0L) stop("No .xlsx files found in data/raw.")

get_gel_number <- function(file) {
  m <- regmatches(basename(file),
                  regexec("gel[^0-9]*([0-9]+)", basename(file),
                          ignore.case = TRUE))[[1]]
  if (length(m) < 2L) stop("Cannot identify gel number: ", basename(file))
  as.integer(m[2])
}
gel_number <- vapply(files, get_gel_number, integer(1), USE.NAMES = FALSE)
ord <- order(gel_number)
files <- files[ord]
gel_number <- gel_number[ord]
if (anyDuplicated(gel_number)) stop("Duplicate gel numbers in data/raw.")
if (!identical(gel_number, seq_along(gel_number))) {
  stop("Gels must be consecutively numbered starting at GEL1.")
}
if (length(files) > ceiling(n_isolates_total / isolates_per_gel)) {
  stop("More gels than expected for 92 isolates.")
}

num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

# Import one gel; returns detected bands and metadata for all expected lanes
import_gel <- function(file, gel_no) {
  gel <- paste0("GEL", gel_no)
  raw <- read_excel(file, skip = 2, col_names = FALSE,
                    .name_repair = "minimal")
  # Expected number of lanes
  expected_lanes <- if (gel_no == 6L) 16L else 19L

  # Validate Image Lab export
  if (ncol(raw) != expected_lanes * columns_per_lane) {
    stop(
      gel, ": expected ", expected_lanes * columns_per_lane,
      " columns but found ", ncol(raw), "."
    )
  }

  first_id <- (gel_no - 1L) * isolates_per_gel + 1L
  n_expected <- min(isolates_per_gel,
                    n_isolates_total - first_id + 1L)
  if (n_expected < 1L) stop("No expected isolates for ", gel)

  metadata <- data.frame(
    gel = gel,
    lane = seq_len(n_expected) + 2L,
    isolate = paste0("I", seq.int(first_id, length.out = n_expected)),
    stringsAsFactors = FALSE
  )

  lane_tables <- lapply(seq.int(2L, n_expected + 2L), function(lane) {
    start <- (lane - 1L) * columns_per_lane + 1L
    x <- raw[, start:(start + columns_per_lane - 1L)]
    names(x) <- c("lane_export", "band_number", "band_label", "size_kb",
                  "relative_front", "adjusted_volume", "volume",
                  "absolute_quantity", "relative_quantity", "band_percent",
                  "lane_percent")
    x$band_number <- num(x$band_number)
    x <- x[!is.na(x$band_number), , drop = FALSE]
    if (nrow(x) == 0L) return(NULL)

    x$gel <- gel
    x$lane <- lane
    x$sample_type <- if (lane == 2L) "Reference" else "Isolate"
    x$isolate <- if (lane == 2L) "I7-REF" else paste0("I", first_id + lane - 3L)
    x$band_number <- as.integer(x$band_number)
    x$size_kb <- num(x$size_kb)
    x$size_bp <- x$size_kb * 1000
    x$relative_front <- num(x$relative_front)
    x$adjusted_volume <- num(x$adjusted_volume)
    x$volume <- num(x$volume)
    x$band_percent <- num(x$band_percent)
    x$lane_percent <- num(x$lane_percent)
    x[, c("gel", "lane", "sample_type", "isolate", "band_number",
          "size_bp", "size_kb", "relative_front", "adjusted_volume",
          "volume", "band_percent", "lane_percent")]
  })
  list(bands = dplyr::bind_rows(lane_tables), metadata = metadata)
}

# 1. Import and consolidate
imports <- lapply(seq_along(files), function(i) import_gel(files[i], gel_number[i]))
bands_all <- dplyr::bind_rows(lapply(imports, `[[`, "bands"))
metadata <- dplyr::bind_rows(lapply(imports, `[[`, "metadata"))

if (nrow(bands_all) == 0L) stop("No detected bands in input files.")
if (any(!is.finite(bands_all$size_bp)) || any(bands_all$size_bp <= 0)) {
  stop("Invalid molecular sizes detected. Review the Image Lab exports.")
}
if (any(!is.finite(bands_all$relative_front)) ||
    any(bands_all$relative_front < 0 | bands_all$relative_front > 1)) {
  stop("Invalid relative front values detected.")
}
if (anyDuplicated(bands_all[c("gel", "lane", "band_number")])) {
  stop("Duplicated band numbers within a lane.")
}

bands_isolates <- bands_all %>%
  filter(sample_type == "Isolate") %>%
  arrange(as.integer(sub("^I", "", isolate)), band_number)
bands_reference <- bands_all %>%
  filter(sample_type == "Reference") %>%
  arrange(gel, band_number)

metadata <- metadata %>%
  left_join(bands_isolates %>% count(gel, lane, isolate, name = "n_bands"),
            by = c("gel", "lane", "isolate")) %>%
  mutate(n_bands = ifelse(is.na(n_bands), 0L, n_bands)) %>%
  filter(n_bands > 0L) %>%
  arrange(as.integer(sub("^I", "", isolate)))

if (nrow(metadata) < 2L) stop("Fewer than two isolates with detected bands.")
if (anyDuplicated(metadata$isolate)) stop("Repeated isolate identifiers.")
stopifnot(all(bands_isolates$isolate %in% metadata$isolate))

# 2. I7-REF inter-gel alignment
counts <- bands_reference %>% count(gel, name = "n")
if (nrow(counts) != length(files) || any(counts$n != 5L) ||
    anyDuplicated(bands_reference[c("gel", "band_number")]) ||
    !all(vapply(split(bands_reference$band_number, bands_reference$gel),
                function(x) identical(sort(x), 1:5), logical(1)))) {
  stop("Each gel must have exactly five I7-REF bands numbered 1 to 5.")
}

target <- bands_reference %>%
  group_by(band_number) %>%
  summarise(target_bp = exp(mean(log(size_bp))), .groups = "drop")
reference_alignment <- bands_reference %>%
  left_join(target, by = "band_number") %>%
  mutate(deviation_percent = 100 * (size_bp / target_bp - 1),
         log_correction = log(target_bp) - log(size_bp))

correct_sizes <- function(size, ref_gel) {
  ref_gel <- ref_gel[order(ref_gel$size_bp), , drop = FALSE]
  delta <- approx(x = log(ref_gel$size_bp), y = ref_gel$log_correction,
                  xout = log(size), method = "linear", rule = 2)$y
  exp(log(size) + delta)
}

bands_isolates$size_original_bp <- bands_isolates$size_bp
bands_isolates$size_corrected_bp <- bands_isolates$size_bp
for (g in unique(bands_isolates$gel)) {
  idx <- bands_isolates$gel == g
  ref_gel <- reference_alignment[reference_alignment$gel == g, , drop = FALSE]
  bands_isolates$size_corrected_bp[idx] <-
    correct_sizes(bands_isolates$size_original_bp[idx], ref_gel)
}

# 3. Band matching: fixed 5% relative size tolerance
# Complete linkage ensures max/min band size within a class <= 1.05.
if (nrow(bands_isolates) < 2L) stop("Not enough bands to group.")
distance_bands <- dist(matrix(log(bands_isolates$size_corrected_bp), ncol = 1))
band_tree <- hclust(distance_bands, method = "complete")
bands_isolates$cluster <- cutree(band_tree, h = log1p(tolerance))

classes <- bands_isolates %>%
  group_by(cluster) %>%
  summarise(median_bp = median(size_corrected_bp),
            min_bp = min(size_corrected_bp),
            max_bp = max(size_corrected_bp),
            n_bands = n(), .groups = "drop") %>%
  arrange(desc(median_bp)) %>%
  mutate(band_class = sprintf("B%03d", row_number()))

assignments <- bands_isolates %>%
  left_join(classes %>% select(cluster, band_class), by = "cluster") %>%
  select(-cluster)

conflicts <- assignments %>%
  count(isolate, band_class, name = "n_bands") %>%
  filter(n_bands > 1L)
if (nrow(conflicts) > 0L) {
  write.csv2(conflicts, file.path(tables_path, "matching_conflicts.csv"),
             row.names = FALSE)
  stop("Multiple bands from the same isolate were assigned to one class. ",
       "Review results/tables/matching_conflicts.csv before proceeding.")
}

# 4. Binary presence/absence matrix
isolate_ids <- metadata$isolate
class_ids <- classes$band_class
binary <- matrix(0L, nrow = length(isolate_ids), ncol = length(class_ids),
                 dimnames = list(isolate_ids, class_ids))
for (i in seq_len(nrow(assignments))) {
  binary[assignments$isolate[i], assignments$band_class[i]] <- 1L
}
stopifnot(all(rowSums(binary) > 0),
          identical(as.integer(rowSums(binary)), as.integer(metadata$n_bands)))

# 5. Dice similarity and UPGMA
shared <- tcrossprod(binary)
n_bands <- rowSums(binary)
similarity <- 2 * shared / outer(n_bands, n_bands, "+")
diag(similarity) <- 1
dimnames(similarity) <- list(rownames(binary), rownames(binary))
dissimilarity <- 1 - similarity
hc <- hclust(as.dist(dissimilarity), method = "average")

# Export tables in Brazilian CSV format
write.csv2(bands_all, file.path(tables_path, "bands_all.csv"), row.names = FALSE)
write.csv2(bands_isolates, file.path(tables_path, "bands_isolates.csv"), row.names = FALSE)
write.csv2(bands_reference, file.path(tables_path, "bands_reference.csv"), row.names = FALSE)
write.csv2(metadata, file.path(tables_path, "isolates_metadata.csv"), row.names = FALSE)
write.csv2(reference_alignment, file.path(tables_path, "reference_alignment.csv"), row.names = FALSE)
write.csv2(assignments, file.path(tables_path, "band_assignments.csv"), row.names = FALSE)
write.csv2(classes, file.path(tables_path, "band_classes.csv"), row.names = FALSE)

write_square_csv <- function(x, name) {
  write.csv2(data.frame(isolate = rownames(x), x, check.names = FALSE),
             file.path(tables_path, name), row.names = FALSE)
}
write_square_csv(binary, "binary_matrix.csv")
write_square_csv(similarity, "dice_similarity.csv")
write_square_csv(dissimilarity, "dice_distance.csv")
saveRDS(hc, file.path(tables_path, "upgma.rds"))
write.csv2(data.frame(order = seq_along(hc$order),
                      isolate = hc$labels[hc$order]),
           file.path(tables_path, "dendrogram_order.csv"), row.names = FALSE)

# Figures
p_bands <- ggplot(metadata, aes(x = factor(isolate, levels = isolate_ids), y = n_bands)) +
  geom_col(fill = "steelblue4") +
  facet_wrap(~gel, scales = "free_x", nrow = 1) +
  labs(title = "BOX-PCR - Detected bands per isolate",
       x = "Isolate", y = "Number of bands") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
ggsave(file.path(figures_path, "bands_per_isolate.png"),
       plot = p_bands, width = max(12, 2.5 * length(files)), height = 5, dpi = 300)

p_ref <- ggplot(reference_alignment,
                aes(x = factor(band_number), y = deviation_percent,
                    color = gel, group = gel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line() + geom_point(size = 2) + theme_bw() +
  labs(title = "BOX-PCR - I7-REF consistency", x = "Reference band",
       y = "Deviation from shared reference (%)", color = "Gel")
ggsave(file.path(figures_path, "reference_alignment.png"),
       plot = p_ref, width = 8, height = 5, dpi = 300)

png(file.path(figures_path, "upgma_dice_5pct.png"),
    width = 4400, height = 2400, res = 300)
par(mar = c(7, 5, 4, 2))
plot(hc, main = "BOX-PCR - Dice / UPGMA (5%)", xlab = "", sub = "",
     ylab = "Dice dissimilarity", cex = if (nrow(binary) > 60) 0.50 else 0.70,
     hang = -1)
dev.off()

order_ids <- hc$labels[hc$order]
heatmap_data <- as.data.frame(as.table(similarity[order_ids, order_ids, drop = FALSE]))
names(heatmap_data) <- c("Isolate1", "Isolate2", "Dice")
heatmap_data$Isolate1 <- factor(heatmap_data$Isolate1, levels = order_ids)
heatmap_data$Isolate2 <- factor(heatmap_data$Isolate2, levels = rev(order_ids))
p_heatmap <- ggplot(heatmap_data, aes(Isolate1, Isolate2, fill = Dice)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 1)) +
  coord_equal() + theme_bw() +
  labs(title = "BOX-PCR - Dice similarity (5%)", x = "Isolate", y = "Isolate") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5))
ggsave(file.path(figures_path, "dice_heatmap_5pct.png"),
       plot = p_heatmap, width = 14, height = 12, dpi = 300)

# Final summary
pairwise <- similarity[upper.tri(similarity)]
cat("\nBOX-PCR COMPLETE PIPELINE\n")
cat("-----------------------------\n")
cat("Gels:", length(files), "\n")
cat("Isolates retained:", nrow(binary), "\n")
cat("Detected experimental bands:", nrow(bands_isolates), "\n")
cat("Band classes (5%):", ncol(binary), "\n")
cat("Mean Dice similarity:", round(mean(pairwise), 4), "\n")
cat("Identical binary-profile pairs:", sum(pairwise == 1), "\n")
cat("Output tables:", tables_path, "\n")
cat("Output figures:", figures_path, "\n")

# Similarity thresholds
similarities <- seq(0.50, 1.00, by = 0.05)

# Count groups at each threshold
data.frame(
  Similarity = similarities,
  Groups = sapply(similarities, function(s) {
    length(unique(cutree(hc, h = 1 - s)))
  })
)
