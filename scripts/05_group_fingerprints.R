################################################################################
# BOX-PCR | 05_group_fingerprints.R
# Panels of BOX-PCR fingerprints by UPGMA group (Dice dissimilarity = 0.50).
# Horizontal crop: measured lane centers; vertical alignment: I7-REF.
################################################################################

library(magick)

# Paths and parameters ------------------------------------------------------------
project_path <- getwd()
image_dir <- file.path(project_path, "data", "gel_images")
tables_dir <- file.path(project_path, "results", "tables")
output_dir <- file.path(project_path, "results", "figures", "fingerprint_groups")
panels_dir <- file.path(output_dir, "aligned")
dir.create(panels_dir, recursive = TRUE, showWarnings = FALSE)

threshold_dissimilarity <- 0.50
lane_width <- 33L
crop_top <- 155L
crop_bottom <- 560L
header <- 32L
footer <- 25L
padding <- 7L
gutter <- 9L

# Five REF-band positions (top to bottom) measured in the original TIFFs.
settings <- data.frame(
  gel = paste0("GEL", 1:6),
  n_lanes = c(19L, 19L, 19L, 19L, 19L, 16L),
  ref_1 = c(220, 203, 199, 207, 220, 205),
  ref_2 = c(371, 349, 343, 360, 370, 350),
  ref_3 = c(387, 365, 359, 375, 386, 366),
  ref_4 = c(456, 435, 429, 445, 457, 437),
  ref_5 = c(488, 469, 463, 478, 491, 472)
)

# Read data -----------------------------------------------------------------------
files <- setNames(file.path(image_dir, paste0(settings$gel, ".tif")), settings$gel)
if (any(!file.exists(files))) stop("Missing TIFF files: ", paste(basename(files[!file.exists(files)]), collapse = ", "))
images <- lapply(files, image_read)

hc <- readRDS(file.path(tables_dir, "upgma.rds"))
meta <- read.csv2(file.path(tables_dir, "isolates_metadata.csv"), stringsAsFactors = FALSE)
centers <- read.csv2(file.path(image_dir, "lane_centers_measured.csv"), stringsAsFactors = FALSE)

stopifnot(inherits(hc, "hclust"),
          all(c("isolate", "gel", "lane") %in% names(meta)),
          all(c("gel", "lane", "center_x") %in% names(centers)),
          !anyDuplicated(meta$isolate),
          !anyDuplicated(paste(centers$gel, centers$lane)),
          setequal(meta$isolate, hc$labels))

meta <- meta[match(hc$labels, meta$isolate), , drop = FALSE]
meta$lane <- as.integer(meta$lane)
meta$group <- as.integer(cutree(hc, h = threshold_dissimilarity)[meta$isolate])
leaf_order <- hc$labels[hc$order]
meta <- meta[order(meta$group, match(meta$isolate, leaf_order)), , drop = FALSE]

# Geometry and REF alignment ------------------------------------------------------
geometry <- setNames(lapply(seq_len(nrow(settings)), function(i) {
  s <- settings[i, ]
  gel <- s$gel
  info <- image_info(images[[gel]])
  if (info$width != 1726 || info$height != 832)
    stop(gel, ": unexpected TIFF dimensions")

  calibrated <- centers[centers$gel == gel, , drop = FALSE]
  calibrated <- calibrated[match(seq_len(s$n_lanes), calibrated$lane), , drop = FALSE]
  if (anyNA(calibrated$center_x) || anyNA(calibrated$lane))
    stop(gel, ": incomplete measured lane centers")

  x <- as.numeric(calibrated$center_x)
  left <- as.integer(round(x - lane_width / 2))
  if (any(left < 0 | left + lane_width > info$width))
    stop(gel, ": lane crop outside image")

  ref <- as.numeric(unlist(s[paste0("ref_", 1:5)], use.names = FALSE))
  if (any(diff(ref) <= 0) || min(ref) <= crop_top || max(ref) >= crop_bottom)
    stop(gel, ": invalid REF coordinates")

  list(left = left, ref = ref)
}), settings$gel)

for (i in seq_len(nrow(meta))) {
  g <- meta$gel[i]
  valid_lanes <- if (g == "GEL6") 3:14 else 3:18
  if (!(g %in% names(geometry)) || !(meta$lane[i] %in% valid_lanes))
    stop("Invalid gel/lane for ", meta$isolate[i])
}

ref_targets <- apply(do.call(rbind, lapply(geometry, `[[`, "ref")), 2, median)

align_lane <- function(img, ref) {
  src <- as.integer(c(crop_top, ref, crop_bottom))
  dst <- as.integer(round(c(crop_top, ref_targets, crop_bottom)))
  parts <- lapply(seq_len(length(src) - 1L), function(j) {
    section <- image_crop(img, geometry_area(
      lane_width, src[j + 1L] - src[j], 0L, src[j] - crop_top
    ))
    image_resize(section, sprintf("%dx%d!", lane_width, dst[j + 1L] - dst[j]),
                 filter = "Triangle")
  })
  image_append(image_join(parts), stack = TRUE)
}

# Extract and assemble ------------------------------------------------------------
strip_width <- lane_width + 2L * padding
strip_height <- header + (crop_bottom - crop_top) + footer

make_strip <- function(gel, lane, label) {
  piece <- image_crop(images[[gel]], geometry_area(
    lane_width, crop_bottom - crop_top, geometry[[gel]]$left[lane], crop_top
  ))
  piece <- align_lane(piece, geometry[[gel]]$ref)
  canvas <- image_blank(strip_width, strip_height, color = "white")
  canvas <- image_composite(canvas, piece, offset = sprintf("+%d+%d", padding, header))
  canvas <- image_annotate(canvas, label, size = 12, color = "black",
                           gravity = "north", location = "+0+7")
  image_annotate(canvas, gel, size = 10, color = "gray40",
                 gravity = "south", location = "+0+4")
}

assemble <- function(strips) {
  gap <- image_blank(gutter, strip_height, color = "white")
  parts <- list()
  for (i in seq_along(strips)) {
    if (i > 1L) parts[[length(parts) + 1L]] <- gap
    parts[[length(parts) + 1L]] <- strips[[i]]
  }
  image_append(image_join(parts), stack = FALSE)
}

for (group in sort(unique(meta$group))) {
  members <- meta[meta$group == group, , drop = FALSE]
  strips <- lapply(seq_len(nrow(members)), function(i) {
    make_strip(members$gel[i], members$lane[i], members$isolate[i])
  })
  image_write(assemble(strips), file.path(panels_dir, sprintf("group_%02d.png", group)))
}

# Outputs -------------------------------------------------------------------------
write.csv2(meta[, c("group", "isolate", "gel", "lane")],
           file.path(output_dir, "group_membership.csv"), row.names = FALSE)

ref_strips <- lapply(names(geometry), function(g) make_strip(g, 2L, g))
image_write(assemble(ref_strips), file.path(output_dir, "REF_QC_aligned.png"))

message("Exported ", length(unique(meta$group)), " aligned fingerprint groups to: ", panels_dir)
