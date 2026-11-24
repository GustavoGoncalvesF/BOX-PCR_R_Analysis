#####################################################################################
# Project: BOX-PCR Fingerprint Analysis
# Script: 04_fingerprint_pcoa.R
#
# Description:
# Principal Coordinates Analysis (PCoA) based on Dice dissimilarity,
# with visualization of isolate metadata.
#####################################################################################

# Load packages
library(ggplot2)
library(readxl)

# Project directory
project_path <- getwd()

# Output directory
figures_path <- file.path(project_path, "results", "figures", "pcoa")
dir.create(figures_path, recursive = TRUE, showWarnings = FALSE)

# Load Dice dissimilarity matrix
dice <- read.csv2(
  file.path(project_path, "results", "tables", "dice_distance.csv"),
  row.names = 1,
  check.names = FALSE
)

# Load metadata
metadata <- read_excel(
  file.path(project_path, "data", "metadata", "metadata-box.xlsx")
)

# Convert to data frame
metadata <- as.data.frame(metadata)

# Match metadata with isolates in the distance matrix
isolates <- rownames(dice)

stopifnot(
  identical(rownames(dice), colnames(dice)),
  !anyDuplicated(metadata$Isolate),
  all(isolates %in% metadata$Isolate)
)

metadata <- metadata[
  match(isolates, metadata$Isolate),
]

# Set variables as factors
metadata$BAG <- factor(metadata$BAG)
metadata$Species <- factor(metadata$Species)
metadata$Accession <- factor(metadata$Accession)
metadata$Morphotype <- factor(metadata$Morphotype)

# PCoA based on Dice dissimilarity
pcoa <- cmdscale(
  as.dist(as.matrix(dice)),
  k = 2,
  eig = TRUE,
  add = TRUE
)

# Percentage of variation represented by each axis
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]

variance <- 100 * positive_eigenvalues / sum(positive_eigenvalues)

# Build plotting data frame
pcoa_data <- data.frame(
  Isolate = isolates,
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  BAG = metadata$BAG,
  Species = metadata$Species,
  Accession = metadata$Accession,
  Morphotype = metadata$Morphotype
)

# Plot function
plot_pcoa <- function(variable) {

  p <- ggplot(
    pcoa_data,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = .data[[variable]]
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    labs(
      title = paste("BOX-PCR PCoA -", variable),
      x = paste0("PCoA1 (", round(variance[1], 1), "%)"),
      y = paste0("PCoA2 (", round(variance[2], 1), "%)"),
      color = variable
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )

  ggsave(
    filename = file.path(
      figures_path,
      paste0("pcoa_", tolower(variable), ".png")
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )

  print(p)

  invisible(p)
}

# Generate PCoA plots
pcoa_species <- plot_pcoa("Species")
pcoa_bag <- plot_pcoa("BAG")
pcoa_accession <- plot_pcoa("Accession")
pcoa_morphotype <- plot_pcoa("Morphotype")

# Save PCoA coordinates
write.csv2(
  pcoa_data,
  file.path(figures_path, "pcoa_coordinates.csv"),
  row.names = FALSE
)
