# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Run this script once before running the analysis pipeline.
## Installs all required R packages.

cran_packages <- c(
  "tidyverse",
  "cowplot",
  "ggtext",
  "gridtext",
  "patchwork",
  "MASS",
  "broom",
  "vroom",
  "stringr",
  "ggsignif",
  "foreach"
)

install.packages(cran_packages)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Biostrings",
  "Rsamtools",
  "rtracklayer"
))
