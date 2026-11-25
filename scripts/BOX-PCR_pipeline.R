#####################################################################################
# Project: BOX-PCR Genomic Fingerprinting
# Script: BOX-PCR_pipeline.R
# Author: Gustavo Goncalves
#
# Description:
# Integrated BOX-PCR analysis: Image Lab import, reference alignment,
# 5% band matching, Dice / UPGMA, PCoA, group fingerprints, individual strips
# and PGP selection. Selected-isolate figure assembly is handled separately.
#
# Input:  data/raw/Gel*-Analysis-Table.xlsx
# Output: results/tables/ and results/figures/
#
# Lane layout:
#   GEL1-GEL5: 1=M, 2=REF, 3:18=isolates, 19=M
#   GEL6:      1=M, 2=REF, 3:14=I81-I92, 15=NTC, 16=M
#   Experimental lanes without detected bands are omitted.
#
# Run from the RStudio project root: source("scripts/BOX-PCR_pipeline.R")
#####################################################################################

# Packages
library(readxl)
library(dplyr)
library(ggplot2)

# Optional stages run when their input files are available.
run_pcoa <- TRUE
run_fingerprint_panels <- TRUE
run_isolate_selection <- TRUE

# Fixed clustering / selection parameters
dissimilarity_cut <- 0.60  # 40% Dice similarity
n_target <- 25L

min_informative_bands <- 2L

# Fingerprint panel geometry (original TIFF pixels)
lane_width <- 28L
crop_top <- 155L
crop_bottom <- 560L

# Paths: open BOX-PCR_R_Analysis.Rproj before running
project_path <- getwd()
raw_path <- file.path(project_path, "data", "raw")
tables_path <- file.path(project_path, "results", "tables")
figures_path <- file.path(project_path, "results", "figures")
dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

# Fixed parameters
n_isolates_total <- 92L
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

conflicts_path <- file.path(tables_path, "matching_conflicts.csv")
if (file.exists(conflicts_path)) file.remove(conflicts_path)
conflicts <- assignments %>%
  count(isolate, band_class, name = "n_bands") %>%
  filter(n_bands > 1L)
if (nrow(conflicts) > 0L) {
  write.csv2(conflicts, conflicts_path,
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


# 6. PCoA of Dice dissimilarities -----------------------------------------------
make_pcoa <- function(distance_matrix, hc, sample_metadata, figures_path) {
  if (!all(hc$labels %in% sample_metadata$Isolate))
    stop("PCoA: missing isolate metadata.")
  ids <- hc$labels
  sample_metadata <- sample_metadata[match(ids, sample_metadata$Isolate), , drop = FALSE]
  for (variable in c("Morphotype", "BAG", "Species", "Accession")) {
    if (anyNA(sample_metadata[[variable]]))
      stop("PCoA: missing metadata for ", variable)
    sample_metadata[[variable]] <- factor(sample_metadata[[variable]])
  }
  pcoa_dir <- file.path(figures_path, "pcoa")
  dir.create(pcoa_dir, recursive = TRUE, showWarnings = FALSE)

  ord <- cmdscale(as.dist(distance_matrix[ids, ids, drop = FALSE]),
                  k = 2, eig = TRUE, add = TRUE)
  if (is.null(dim(ord$points)) || ncol(ord$points) < 2L)
    stop("PCoA requires two coordinates.")
  positive <- ord$eig[ord$eig > 0]
  axis_percent <- 100 * ord$eig[1:2] / sum(positive)
  coords <- data.frame(Isolate = ids,
                       PCoA1 = ord$points[, 1], PCoA2 = ord$points[, 2],
                       sample_metadata[, c("Morphotype", "BAG", "Species", "Accession")],
                       row.names = NULL)
  write.csv2(coords, file.path(pcoa_dir, "pcoa_coordinates.csv"), row.names = FALSE)
  xlab <- sprintf("PCoA1 (%.1f%%)", axis_percent[1])
  ylab <- sprintf("PCoA2 (%.1f%%)", axis_percent[2])

  for (variable in c("Species", "BAG", "Accession", "Morphotype")) {
    p <- ggplot(coords, aes(x = PCoA1, y = PCoA2,
                            color = .data[[variable]])) +
      geom_point(size = 2.7, alpha = 0.85) +
      labs(title = paste("BOX-PCR PCoA -", variable), x = xlab, y = ylab,
           color = variable) + theme_bw()
    ggsave(file.path(pcoa_dir, paste0("pcoa_", tolower(variable), ".png")),
           p, width = 8, height = 6, dpi = 300)
  }
  combined <- ggplot(coords, aes(PCoA1, PCoA2, color = Accession, shape = BAG)) +
    geom_point(size = 2.7, alpha = 0.85) +
    labs(title = "BOX-PCR PCoA - Accession and BAG", x = xlab, y = ylab) +
    theme_bw()
  ggsave(file.path(pcoa_dir, "pcoa_accession_bag.png"),
         combined, width = 8, height = 6, dpi = 300)
  message("PCoA: 5 figures exported.")
  invisible(coords)
}

# 7. Fingerprint panels using measured lane centers ----------------------------
make_fingerprint_panels <- function(hc, isolate_metadata, project_path,
                                    figures_path, dissimilarity_cut,
                                    lane_width, crop_top, crop_bottom) {
  if (!requireNamespace("magick", quietly = TRUE))
    stop("Fingerprint panels require the 'magick' package.")
  image_dir <- file.path(project_path, "data", "gel_images")
  output_dir <- file.path(figures_path, "fingerprint_groups")
  panels_dir <- file.path(output_dir, "aligned")
  dir.create(panels_dir, recursive = TRUE, showWarnings = FALSE)

  # Positions of the five REF bands in the original TIFF images.
  settings <- data.frame(
    gel = paste0("GEL", 1:6),
    n_lanes = c(19L, 19L, 19L, 19L, 19L, 16L),
    ref_1 = c(220, 203, 199, 207, 220, 205),
    ref_2 = c(371, 349, 343, 360, 370, 350),
    ref_3 = c(387, 365, 359, 375, 386, 366),
    ref_4 = c(456, 435, 429, 445, 457, 437),
    ref_5 = c(488, 469, 463, 478, 491, 472)
  )
  files <- setNames(file.path(image_dir, paste0(settings$gel, ".tif")), settings$gel)
  if (any(!file.exists(files)))
    stop("Missing TIFFs: ", paste(basename(files[!file.exists(files)]), collapse = ", "))
  images <- lapply(files, magick::image_read)
  centers <- read.csv2(file.path(image_dir, "lane_centers_measured.csv"),
                       stringsAsFactors = FALSE)
  stopifnot(all(c("gel", "lane", "center_x") %in% names(centers)),
            !anyDuplicated(paste(centers$gel, centers$lane)))

  m <- isolate_metadata[, c("gel", "lane", "isolate")]
  stopifnot(setequal(m$isolate, hc$labels), !anyDuplicated(m$isolate))
  m <- m[match(hc$labels, m$isolate), , drop = FALSE]
  m$lane <- as.integer(m$lane)
  m$group <- as.integer(cutree(hc, h = dissimilarity_cut)[m$isolate])
  leaf_order <- hc$labels[hc$order]
  m <- m[order(m$group, match(m$isolate, leaf_order)), , drop = FALSE]

  geometry <- setNames(lapply(seq_len(nrow(settings)), function(i) {
    s <- settings[i, ]; gel <- s$gel
    info <- magick::image_info(images[[gel]])
    if (info$width != 1726L || info$height != 832L)
      stop(gel, ": unexpected TIFF dimensions (expected 1726 x 832).")
    calibrated <- centers[centers$gel == gel, ]
    calibrated <- calibrated[match(seq_len(s$n_lanes), calibrated$lane), ]
    if (anyNA(calibrated$center_x) || anyNA(calibrated$lane))
      stop(gel, ": incomplete measured centers.")
    left <- as.integer(round(as.numeric(calibrated$center_x) - lane_width / 2))
    if (any(left < 0 | left + lane_width > info$width))
      stop(gel, ": crop outside image.")
    ref <- as.numeric(unlist(s[paste0("ref_", 1:5)], use.names = FALSE))
    if (any(diff(ref) <= 0) || min(ref) <= crop_top || max(ref) >= crop_bottom)
      stop(gel, ": invalid REF coordinates.")
    list(left = left, ref = ref)
  }), settings$gel)

  for (i in seq_len(nrow(m))) {
    lanes <- if (m$gel[i] == "GEL6") 3:14 else 3:18
    if (!(m$gel[i] %in% names(geometry)) || !(m$lane[i] %in% lanes))
      stop("Invalid gel/lane for ", m$isolate[i])
  }
  ref_targets <- apply(do.call(rbind, lapply(geometry, `[[`, "ref")), 2, median)
  header <- 32L; footer <- 25L; padding <- 7L; gutter <- 9L
  strip_width <- lane_width + 2L * padding
  strip_height <- header + (crop_bottom - crop_top) + footer

  align_lane <- function(img, ref) {
    src <- as.integer(c(crop_top, ref, crop_bottom))
    dst <- as.integer(round(c(crop_top, ref_targets, crop_bottom)))
    parts <- lapply(seq_len(length(src) - 1L), function(j) {
      section <- magick::image_crop(img, magick::geometry_area(
        lane_width, src[j + 1L] - src[j], 0L, src[j] - crop_top))
      magick::image_resize(section,
                           sprintf("%dx%d!", lane_width, dst[j + 1L] - dst[j]),
                           filter = "Triangle")
    })
    magick::image_append(magick::image_join(parts), stack = TRUE)
  }
  make_strip <- function(gel, lane, label) {
    piece <- magick::image_crop(images[[gel]], magick::geometry_area(
      lane_width, crop_bottom - crop_top, geometry[[gel]]$left[lane], crop_top))
    piece <- align_lane(piece, geometry[[gel]]$ref)
    canvas <- magick::image_blank(strip_width, strip_height, color = "white")
    canvas <- magick::image_composite(canvas, piece,
                                    offset = sprintf("+%d+%d", padding, header))
    canvas <- magick::image_annotate(canvas, label, size = 12, color = "black",
                                     gravity = "north", location = "+0+7")
    magick::image_annotate(canvas, gel, size = 10, color = "gray40",
                           gravity = "south", location = "+0+4")
  }
  assemble <- function(strips) {
    gap <- magick::image_blank(gutter, strip_height, color = "white")
    parts <- list()
    for (i in seq_along(strips)) {
      if (i > 1L) parts[[length(parts) + 1L]] <- gap
      parts[[length(parts) + 1L]] <- strips[[i]]
    }
    magick::image_append(magick::image_join(parts), stack = FALSE)
  }
  for (group in sort(unique(m$group))) {
    sub <- m[m$group == group, , drop = FALSE]
    strips <- lapply(seq_len(nrow(sub)), function(i)
      make_strip(sub$gel[i], sub$lane[i], sub$isolate[i]))
    magick::image_write(assemble(strips),
                        file.path(panels_dir, sprintf("group_%02d.png", group)),
                        format = "png")
  }
  write.csv2(m[, c("group", "isolate", "gel", "lane")],
             file.path(output_dir, "group_membership.csv"), row.names = FALSE)

  # Individual aligned strips for reproducible figure assembly.
  strips_dir <- file.path(output_dir, "strips_aligned")
  dir.create(strips_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(m))) {
    magick::image_write(
      make_strip(m$gel[i], m$lane[i], m$isolate[i]),
      file.path(strips_dir, paste0(m$isolate[i], ".png")), format = "png")
  }

  message("Individual aligned strips: ", nrow(m), " exported to ", strips_dir)

  ref_strips <- lapply(names(geometry), function(g) make_strip(g, 2L, g))
  magick::image_write(assemble(ref_strips),
                      file.path(output_dir, "REF_QC_aligned.png"), format = "png")
  message("Fingerprint panels: ", length(unique(m$group)), " groups exported.")
  invisible(m)
}

# 8. Group summary and preliminary selection for PGP ----------------------------
make_isolate_selection <- function(hc, similarity, binary, sample_metadata,
                                   output_path, dissimilarity_cut, n_target,
                                   min_informative_bands) {
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  ids <- hc$labels
  dice <- similarity[ids, ids, drop = FALSE]
  binary <- binary[ids, , drop = FALSE]
  sample_metadata <- sample_metadata[match(ids, sample_metadata$Isolate), , drop = FALSE]
  if (anyNA(sample_metadata)) stop("Missing isolate metadata for selection.")
  variables <- c("BAG", "Species", "Accession", "Morphotype")
  for (v in variables) sample_metadata[[v]] <- as.character(sample_metadata[[v]])

  cluster <- as.integer(cutree(hc, h = dissimilarity_cut)[ids])
  members <- data.frame(Isolate = ids, Group = cluster,
                        n_bands = rowSums(binary),
                        sample_metadata[, variables], row.names = NULL)
  members <- members[order(members$Group,
                           as.integer(sub("^I", "", members$Isolate))), ]
  groups <- sort(unique(members$Group))

  choose_representative <- function(group_ids) {
    candidates <- group_ids
    valid <- group_ids[rowSums(binary[group_ids, , drop = FALSE]) >=
                         min_informative_bands]
    if (length(valid)) candidates <- valid
    mean_sim <- vapply(candidates, function(x) mean(dice[x, group_ids]), numeric(1))
    band_n <- rowSums(binary[candidates, , drop = FALSE])
    num <- as.integer(sub("^I", "", candidates))
    candidates[order(-mean_sim, -band_n, num)][1]
  }
  summary_rows <- lapply(groups, function(g) {
    sub <- members[members$Group == g, ]
    group_ids <- sub$Isolate
    rep_id <- choose_representative(group_ids)
    rep_row <- members[match(rep_id, members$Isolate), ]
    within <- if (length(group_ids) > 1L) {
      z <- dice[group_ids, group_ids, drop = FALSE]
      z[upper.tri(z)]
    } else NA_real_
    data.frame(Group = g, n_isolates = length(group_ids),
               Members = paste(group_ids, collapse = ", "),
               Representative = rep_id, Rep_n_bands = rep_row$n_bands,
               Mean_within_group_Dice = mean(within),
               Rep_mean_group_Dice = mean(dice[rep_id, group_ids]),
               Low_resolution_group = all(sub$n_bands < min_informative_bands),
               BAG = rep_row$BAG, Species = rep_row$Species,
               Accession = rep_row$Accession, Morphotype = rep_row$Morphotype)
  })
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL

  selected <- summary$Representative
  coverage <- function(reps) {
    mean(apply(dice[ids, reps, drop = FALSE], 1, max))
  }
  keeps_categories <- function(reps) {
    chosen <- sample_metadata[match(reps, sample_metadata$Isolate), , drop = FALSE]
    all(vapply(variables, function(v)
      all(unique(sample_metadata[[v]]) %in% unique(chosen[[v]])), logical(1)))
  }
  if (n_target < 1L) stop("n_target must be >= 1.")
  if (n_target > length(selected))
    warning("n_target exceeds groups: retaining one representative per group.")
  if (!keeps_categories(selected))
    warning("Initial representatives do not cover every metadata category. Review group_summary.csv.")
  while (length(selected) > n_target) {
    options <- selected[vapply(selected, function(r)
      keeps_categories(setdiff(selected, r)), logical(1))]
    if (!length(options)) {
      warning("Cannot reach target without losing a metadata category.")
      break
    }
    scores <- vapply(options, function(r) coverage(setdiff(selected, r)), numeric(1))
    rows <- summary[match(options, summary$Representative), ]
    ord <- order(-scores, rows$Rep_n_bands, rows$n_isolates, rows$Group)
    selected <- setdiff(selected, options[ord[1]])
  }
  summary$Selection <- ifelse(summary$Representative %in% selected,
                              "Selected (provisional)", "Reserve")
  summary$Nearest_selected <- vapply(summary$Representative, function(rep_id) {
    selected[which.max(dice[rep_id, selected])]
  }, character(1))
  summary$Dice_to_nearest_selected <- vapply(seq_len(nrow(summary)), function(i)
    dice[summary$Representative[i], summary$Nearest_selected[i]], numeric(1))
  summary$Review <- ifelse(summary$Low_resolution_group,
                           "Low resolution: all fingerprints have 1 band",
                           ifelse(summary$Rep_n_bands <= 2,
                                  "Review limited number of bands", ""))
  members$Selection <- summary$Selection[match(members$Group, summary$Group)]
  members$Representative <- summary$Representative[match(members$Group, summary$Group)]
  members$Role <- ifelse(members$Isolate == members$Representative,
                         "Representative", "Other member")
  write.csv2(members, file.path(output_path, "group_membership.csv"), row.names = FALSE)
  write.csv2(summary, file.path(output_path, "group_summary.csv"), row.names = FALSE)
  write.csv2(summary[summary$Selection == "Selected (provisional)", ],
             file.path(output_path, sprintf("selected_%02d_provisional.csv", n_target)), row.names = FALSE)
  write.csv2(summary[summary$Selection == "Reserve", ],
             file.path(output_path, "reserve_groups.csv"), row.names = FALSE)

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
         subtitle = sprintf("Dice dissimilarity: %.2f | %d groups | %d proposed representatives",
                            dissimilarity_cut, length(groups), length(selected)),
         x = "UPGMA group", y = "Number of isolates", fill = "PGP selection") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          legend.position = "bottom")
  fig_width <- max(15, 0.45 * length(groups))
  ggsave(file.path(output_path, "group_composition_labeled.png"), p,
         width = fig_width, height = 8.5, dpi = 300, limitsize = FALSE)
  ggsave(file.path(output_path, "group_composition_labeled.pdf"), p,
         width = fig_width, height = 8.5, limitsize = FALSE)
  message("Selection: ", length(groups), " groups; ", length(selected),
          " provisional representatives.")
  invisible(summary)
}

# 9. Optional analyses -----------------------------------------------------------
metadata_file <- file.path(project_path, "data", "metadata", "metadata-box.xlsx")
image_dir <- file.path(project_path, "data", "gel_images")
center_file <- file.path(image_dir, "lane_centers_measured.csv")

sample_metadata <- NULL
if (run_pcoa || run_isolate_selection) {
  if (file.exists(metadata_file)) {
    sample_metadata <- as.data.frame(readxl::read_excel(metadata_file))
    needed <- c("Isolate", "Morphotype", "BAG", "Species", "Accession")
    if (!all(needed %in% names(sample_metadata)))
      stop("Metadata missing columns: ", paste(setdiff(needed, names(sample_metadata)), collapse = ", "))
    if (anyDuplicated(sample_metadata$Isolate)) stop("Duplicated Isolate in metadata.")
    if (!all(hc$labels %in% sample_metadata$Isolate))
      stop("Missing metadata for: ", paste(setdiff(hc$labels, sample_metadata$Isolate), collapse = ", "))
    sample_metadata <- sample_metadata[match(hc$labels, sample_metadata$Isolate), ]
  } else {
    warning("No data/metadata/isolate_metadata.xlsx: skipping PCoA and isolate selection.")
  }
}

if (run_pcoa && !is.null(sample_metadata))
  make_pcoa(dissimilarity, hc, sample_metadata, figures_path)

if (run_fingerprint_panels) {
  tif_files <- file.path(image_dir, paste0("GEL", 1:6, ".tif"))
  if (length(files) == 6L && all(file.exists(tif_files)) && file.exists(center_file)) {
    make_fingerprint_panels(hc, metadata, project_path, figures_path,
                            dissimilarity_cut, lane_width, crop_top, crop_bottom)
  } else {
    warning("Missing TIFFs / lane_centers_measured.csv (or fewer than 6 gels): skipping panels.")
  }
}

if (run_isolate_selection && !is.null(sample_metadata))
  make_isolate_selection(hc, similarity, binary, sample_metadata,
                         file.path(project_path, "results", "isolate_selection"),
                         dissimilarity_cut, n_target, min_informative_bands)

# 10. Final summary -------------------------------------------------------------
pairwise <- similarity[upper.tri(similarity)]
cat("\nBOX-PCR INTEGRATED PIPELINE\n")
cat("--------------------------------------\n")
cat("Gels:", length(files), "\n")
cat("Isolates retained:", nrow(binary), "\n")
cat("Detected experimental bands:", nrow(bands_isolates), "\n")
cat("Band classes (5%):", ncol(binary), "\n")
cat("Mean Dice similarity:", round(mean(pairwise), 4), "\n")
cat("Groups at Dice dissimilarity", dissimilarity_cut, ":",
    length(unique(cutree(hc, h = dissimilarity_cut))), "\n")
cat("Results:", file.path(project_path, "results"), "\n")

similarity_thresholds <- seq(0.50, 1.00, by = 0.05)
threshold_summary <- data.frame(
  Similarity = similarity_thresholds,
  Groups = vapply(similarity_thresholds,
                  function(s) length(unique(cutree(hc, h = 1 - s))), integer(1))
)
write.csv2(threshold_summary, file.path(tables_path, "group_counts_by_similarity.csv"),
           row.names = FALSE)
print(threshold_summary, row.names = FALSE)
