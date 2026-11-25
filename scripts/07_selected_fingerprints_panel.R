################################################################################
# Project: BOX-PCR Genomic Fingerprinting
# Script: 07_selected_fingerprints_panel.R
# Author: Gustavo Goncalves
# Description: Horizontal panel of selected BOX-PCR fingerprints.
# Input:  results/figures/fingerprint_groups/strips_aligned/I*.png
# Output: results/figures/selected_fingerprints/selected_25_panel.png
################################################################################

library(magick)

project_path <- getwd()
strips_dir <- file.path(project_path, "results", "figures",
                        "fingerprint_groups", "strips_aligned")
output_dir <- file.path(project_path, "results", "figures",
                        "selected_fingerprints")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Selected isolates (left to right)
selected_ids <- c(
  "I1", "I2", "I4", "I8", "I12",
  "I13", "I23", "I30", "I33", "I36",
  "I37", "I41", "I42", "I44", "I46",
  "I47", "I48", "I53", "I59", "I63",
  "I73", "I79", "I80", "I85", "I88"
)

# Figure layout
spacing <- 10L
header_height <- 65L
figure_title <- "BOX-PCR - Selected isolates (40% Dice similarity)"

# Read the original aligned strips (no resizing or contrast adjustment)
strip_files <- file.path(strips_dir, paste0(selected_ids, ".png"))
if (any(!file.exists(strip_files))) {
  stop("Missing strips: ", paste(selected_ids[!file.exists(strip_files)], collapse = ", "))
}
strips <- lapply(strip_files, image_read)
strip_info <- lapply(strips, image_info)
strip_width <- max(vapply(strip_info, function(x) x$width[1], numeric(1)))
strip_height <- max(vapply(strip_info, function(x) x$height[1], numeric(1)))

# Place strips on identical white tiles without scaling image pixels
strips <- lapply(strips, function(im) {
  d <- image_info(im)
  tile <- image_blank(strip_width, strip_height, color = "white")
  image_composite(tile, im, offset = sprintf("+%d+%d",
    (strip_width - d$width[1]) %/% 2L,
    (strip_height - d$height[1]) %/% 2L))
})

# All 25 isolates in a SINGLE horizontal row
parts <- list()
for (i in seq_along(strips)) {
  if (i > 1L) {
    parts[[length(parts) + 1L]] <-
      image_blank(spacing, strip_height, color = "white")
  }
  parts[[length(parts) + 1L]] <- strips[[i]]
}
panel <- image_append(image_join(parts), stack = FALSE)

# Full-width header avoids clipping of the title
panel_width <- image_info(panel)$width[1]
title_size <- min(28L, max(12L, floor(panel_width / (nchar(figure_title) * 0.65))))
header <- image_blank(panel_width, header_height, color = "white")
header <- image_annotate(header, figure_title, size = title_size,
                         color = "black", gravity = "center")
final_panel <- image_append(image_join(list(header, panel)), stack = TRUE)

output_file <- file.path(output_dir, "selected_25_panel.png")
image_write(final_panel, output_file, format = "png")
cat("Saved:", output_file, "\n")
cat("Dimensions:", image_info(final_panel)$width, "x",
    image_info(final_panel)$height, "pixels\n")
