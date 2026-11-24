#####################################################################################
# Project: BOX-PCR Genomic Fingerprinting
# Script: 01_import_data.R
# Author: Gustavo Goncalves
#
# Description:
# Automatically import and consolidate band detection tables exported
# from Bio-Rad Image Lab.
#
# Lane layout:
#   Lane 1: DNA Ladder
#   Lane 2: Reference (I7-REF)
#   Lanes 3-18: Experimental isolates
#   Lane 19: DNA Ladder
#
# Notes:
#   - Input files are automatically identified in data/raw.
#   - Isolates without detected bands are classified as No amplification.
#   - Ladder measurements are excluded from the output.
#   - Reference measurements are retained for inter-gel comparisons.
#
# Input:  data/raw/*.xlsx
# Output: bands_all.csv, bands_isolates.csv,
#         bands_reference.csv, isolates_metadata.csv
#####################################################################################

# Set up packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Project directories
project_path <- getwd()
raw_path <- file.path(project_path, "data", "raw")
tables_path <- file.path(project_path, "results", "tables")
figures_path <- file.path(project_path, "results", "figures")

dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

# Experimental configuration
total_isolates <- 92L
n_lanes <- 19L
columns_per_lane <- 11L
n_isolates_per_gel <- n_lanes - 3L

# Identify input files
gel_files <- list.files(
  path = raw_path,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Exclude temporary Excel files
gel_files <- gel_files[
  !grepl("^~\\$", basename(gel_files))
]

if (length(gel_files) == 0) {
  stop("No Excel files found in: ", raw_path)
}

# Extract gel numbers from file names
get_gel_number <- function(file) {
  name <- basename(file)
  match <- regmatches(
    name,
    regexec("gel[^0-9]*([0-9]+)", name,
            ignore.case = TRUE)
  )[[1]]

  if (length(match) < 2) {
    stop("Could not identify gel number: ", name)
  }

  as.integer(match[2])
}

gel_numbers <- vapply(
  gel_files, get_gel_number, integer(1)
)

# Sort gels numerically
order_idx <- order(gel_numbers)
gel_files <- gel_files[order_idx]
gel_numbers <- gel_numbers[order_idx]

if (anyDuplicated(gel_numbers)) {
  stop("Duplicate gel numbers detected.")
}

# Gel numbers must be sequential from GEL1
if (!identical(
  unname(gel_numbers),
  seq_len(length(gel_numbers))
)) {
  stop("Gel numbers must be sequential, starting at GEL1.")
}

if (length(gel_files) * n_isolates_per_gel >
    total_isolates + n_isolates_per_gel - 1L) {
  stop("More gels than expected for the configured isolates.")
}

# Convert numeric measurements
to_number <- function(x) {
  suppressWarnings(
    as.numeric(gsub(",", ".", as.character(x)))
  )
}

# Import a single gel
import_gel <- function(file, gel_number) {

  gel <- paste0("GEL", gel_number)

  raw <- read_excel(
    path = file,
    skip = 2,
    col_names = FALSE,
    .name_repair = "minimal"
  )

  expected_columns <- n_lanes * columns_per_lane

  if (ncol(raw) != expected_columns) {
    stop(
      gel, ": expected ", expected_columns,
      " columns, found ", ncol(raw)
    )
  }

  # Number of experimental isolates in this gel
  first_isolate <- (
    (gel_number - 1L) * n_isolates_per_gel
  ) + 1L

  n_experimental <- min(
    n_isolates_per_gel,
    max(0L, total_isolates - first_isolate + 1L)
  )

  # Create metadata for all expected isolates
  metadata <- data.frame(
    gel = gel,
    lane = seq_len(n_experimental) + 2L,
    isolate = paste0(
      "I",
      seq.int(first_isolate, length.out = n_experimental)
    )
  )

  # Extract each lane
  lane_tables <- lapply(seq_len(n_lanes), function(lane) {

    start_col <- (lane - 1L) * columns_per_lane + 1L
    end_col <- start_col + columns_per_lane - 1L

    lane_data <- raw[, start_col:end_col]

    names(lane_data) <- c(
      "lane_export",
      "band_number",
      "band_label",
      "size_kb",
      "relative_front",
      "adjusted_volume",
      "volume",
      "absolute_quantity",
      "relative_quantity",
      "band_percent",
      "lane_percent"
    )

    # Skip DNA ladders and unused final lanes
    if (lane == 1L || lane == n_lanes ||
        (lane >= 3L && lane > n_experimental + 2L)) {
      return(NULL)
    }

    sample_type <- if (lane == 2L) {
      "Reference"
    } else {
      "Isolate"
    }

    isolate <- if (lane == 2L) {
      "I7-REF"
    } else {
      paste0("I", first_isolate + lane - 3L)
    }

    lane_data %>%
      mutate(
        gel = gel,
        lane = lane,
        sample_type = sample_type,
        isolate = isolate,
        band_number = to_number(band_number),
        size_kb = to_number(size_kb),
        relative_front = to_number(relative_front),
        adjusted_volume = to_number(adjusted_volume),
        volume = to_number(volume),
        band_percent = to_number(band_percent),
        lane_percent = to_number(lane_percent)
      ) %>%
      filter(!is.na(band_number)) %>%
      mutate(
        band_number = as.integer(band_number),
        size_bp = size_kb * 1000
      ) %>%
      select(
        gel, lane, sample_type, isolate,
        band_number, size_bp, size_kb,
        relative_front, adjusted_volume, volume,
        band_percent, lane_percent
      )
  })

  list(
    bands = bind_rows(lane_tables),
    metadata = metadata
  )
}

# Import all gels
gel_data <- lapply(seq_along(gel_files), function(i) {
  import_gel(
    file = gel_files[i],
    gel_number = gel_numbers[i]
  )
})

# Consolidate all detected bands
bands_all <- bind_rows(
  lapply(gel_data, function(x) x$bands)
)

# Consolidate isolate metadata
isolates_metadata <- bind_rows(
  lapply(gel_data, function(x) x$metadata)
) %>%
  left_join(
    bands_all %>%
      filter(sample_type == "Isolate") %>%
      count(gel, lane, isolate, name = "n_bands"),
    by = c("gel", "lane", "isolate")
  ) %>%
  mutate(
    n_bands = replace_na(n_bands, 0L)
  ) %>%
  filter(n_bands > 0) %>%
  arrange(
    as.integer(sub("^I", "", isolate))
  )

# Experimental bands
bands_isolates <- bands_all %>%
  filter(sample_type == "Isolate") %>%
  arrange(gel, lane, band_number)

# Reference bands for inter-gel comparisons
bands_reference <- bands_all %>%
  filter(sample_type == "Reference") %>%
  arrange(gel, band_number)

# Quality control
stopifnot(
  !anyDuplicated(isolates_metadata$isolate),
  all(!is.na(bands_all$size_bp) &
        bands_all$size_bp > 0),
  all(!is.na(bands_all$relative_front) &
        bands_all$relative_front >= 0 &
        bands_all$relative_front <= 1)
)

# Export tables
write.csv2(
  bands_all,
  file.path(tables_path, "bands_all.csv"),
  row.names = FALSE
)

write.csv2(
  bands_isolates,
  file.path(tables_path, "bands_isolates.csv"),
  row.names = FALSE
)

write.csv2(
  bands_reference,
  file.path(tables_path, "bands_reference.csv"),
  row.names = FALSE
)

write.csv2(
  isolates_metadata,
  file.path(tables_path, "isolates_metadata.csv"),
  row.names = FALSE
)

# Summary
cat("\nBOX-PCR IMPORT SUMMARY\n")
cat("--------------------------------\n")
cat("Gels imported:", length(gel_files), "\n")
cat("Isolates retained:", nrow(isolates_metadata), "\n")
cat("Experimental bands:", nrow(bands_isolates), "\n")
cat("Reference bands:", nrow(bands_reference), "\n")

print(
  isolates_metadata %>%
    group_by(gel) %>%
    summarise(
      isolates = n(),
      total_bands = sum(n_bands),
      .groups = "drop"
    )
)

# QC plot: detected bands per isolate
p <- ggplot(
  isolates_metadata,
  aes(
    x = factor(
      isolate,
      levels = isolates_metadata$isolate
    ),
    y = n_bands
  )
) +
  geom_col(fill = "steelblue4") +
  facet_wrap(~ gel, scales = "free_x", nrow = 1) +
  labs(
    title = "BOX-PCR - Detected bands per isolate",
    x = "Isolate",
    y = "Number of detected bands"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

print(p)

ggsave(
  filename = file.path(figures_path, "bands_per_isolate.png"),
  plot = p,
  width = 14,
  height = 5,
  dpi = 300
)
