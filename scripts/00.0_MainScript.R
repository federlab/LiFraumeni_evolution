# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Li-Fraumeni evolution main script

## Create output directories ####
output_dirs <- c(
  "results",
  "results/Manuscript_figures",
  "results/Manuscript_figures/Fig_1",
  "results/Manuscript_figures/Fig_2",
  "results/Manuscript_figures/Fig_3",
  "results/Manuscript_figures/Fig_4",
  "results/Manuscript_figures/Fig_S1",
  "results/Manuscript_figures/Fig_S2",
  "results/Manuscript_figures/Fig_S3",
  "results/Manuscript_figures/Fig_S4",
  "results/Manuscript_figures/Fig_S5",
  "results/Manuscript_figures/Fig_S6",
  "results/Manuscript_figures/Fig_S7",
  "results/Manuscript_figures/Fig_S8",
  "results/Manuscript_figures/Fig_S9",
  "results/Manuscript_figures/revisions",
  "results/Manuscript_tables"
)
invisible(lapply(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

## Load Modules ####
library("tidyverse")
library("cowplot")
library("ggtext")
library("tidyr")
library("rtracklayer")
library("dplyr")
library("ggplot2")
library("gridtext")
library("tibble")
library("foreach")
library("MASS")
library("broom")
library("patchwork")
library("stringr")
library("vroom")
library("Rsamtools")
library("Biostrings")



# Load inputs list
inputs <- read_delim("processing_config.txt", delim="=", col_names = FALSE, comment = "#", skip_empty_rows = TRUE) %>%
  pivot_wider(names_from = X1, values_from = X2) %>% 
  type_convert()

## sample information
sample_data <- read_csv(inputs$extra_data_file)

## Load Functions####
source("scripts/fn_load_maf.R")
source("scripts/fn_load_bed.R")
source("scripts/fn_check_bed_overlap.R")
source("scripts/color_scale.R")

## Shared constants (tissue lists, subject groupings, gene panels, palettes) ####
source("scripts/00.1_SharedConstants.R")

# variant classification table renames subvariants to "intron", "exon", etc.
variant_clasification_table = read_delim(inputs$varClassTrans) 

## Load AlphaMissense file ####
alphamissense <- read_delim(inputs$AlphaMissense, delim = "\t", skip = 3) %>% 
  dplyr::rename('CHROM' = "#CHROM")

## Load IntOGen file
intOGen <- read_delim(inputs$intOGen, delim = "\t") %>% 
  dplyr::select(chr,pos,alt, boostDM_score, boostDM_class) %>% 
  mutate(chr = paste0("chr",chr))

# Read in BED for MUT
targ_bed12_MUT <- loadBed(inputs$bed12_MUT)
targ_bed6_MUT <- loadBed(inputs$bed6_MUT)

# Convert genome.mut to depth file
source("scripts/data_mut_to_depth.R")
source("scripts/combine_maf.R")
# Read in Depth
source("scripts/data_load_depth.R")
## filter out samples that do not contain depth requirements
source("scripts/data_filter_depth.R")

# Read in MAF
source("scripts/data_load_maf.R")
source("scripts/data_filter_maf.R")

# assign coding regions and repeat masking filter
source("scripts/annotate_coding_depth.R")

# filter variants that existed at high VAF in other samples
source("scripts/cross_contamination_filter.R")

######## Depth QC (Fig 1, S1, Table S4)
source("scripts/qc_depth_by_panel.R")
source("scripts/qc_depth_visualization.R")

######## CHIP blood analysis (Fig 2, S2, S3, S4, S5)
## MF regression: 2BC
source("scripts/chip_mf_coding_noncoding.R")
## multivariate regression coefficients: 2DE
source("scripts/chip_mf_multivar_coef.R")
## mutagenesis MF and regression excluding CTx subjects: S4
source("scripts/supp_chip_mf_no_ctx.R")
## per-gene regression plots: S2
coding = FALSE
source("scripts/supp_regression_genes.R")
coding = TRUE
source("scripts/supp_regression_genes.R")
## Table S3
source("scripts/supp_table_overlap.R")
## clone size by VAF: 3A (sourced here; provides objects used by chip_mutation_count)
source("scripts/tp53_skyscraper_blood.R")
## mutation counts per patient: 2A
source("scripts/chip_mutation_count.R")
## mutation signatures: 2F
source("scripts/chip_mutsigs.R")
## SBS strand spectra: 2G
source("scripts/chip_sbsg_spectra.R")
## dN/dS all CHIP genes: 2H
source("scripts/chip_dnds_target_sites.R")
source("scripts/chip_dnds_target_annotate.R")
source("scripts/chip_dnds_all_genes.R")
## mutation burden regression: S3
source("scripts/supp_chip_mutation_burden.R")
## CC>TT DNVs with SBS signatures in blood: S5
source("scripts/supp_blood_dnv_mutsigs.R")

######## TP53 blood analysis (Fig 3, S6)
## MF regression: 3BC
source("scripts/tp53_mf_coding_noncoding.R")
## multivariate regression coefficients: 3DE
source("scripts/tp53_mf_multivar_coef.R")
## lollipop: 3B
source("scripts/tp53_lollipop_blood.R")
## coding/noncoding MF ratio: 3F
source("scripts/tp53_mf_ratio.R")
## dN/dS by LFS/chemo group: 3J
source("scripts/tp53_dnds_blood.R")
## DBD AlphaMissense pathogenicity: 3C/D & 4F/G
source("scripts/tp53_binding_domain.R")
## AlphaMissense annotation + grouped plot: 3K, S6
source("scripts/tp53_annotate_alphamissense.R")
source("scripts/supp_tp53_alphamissense.R")

######## Tissue analysis (Fig 4, S7, S8, S9)
## clone size by VAF across tissues: 4A
source("scripts/tissue_skyscraper.R")
## SNV/indel proportions: 4H
source("scripts/tissue_variant_types.R")
## mutation sharing heatmap: 4A/4G
source("scripts/tissue_mutation_overlap.R")
## cross-tissue contamination filter: S8 (must precede tissue_dnv, which overwrites skyscraper_prep)
source("scripts/supp_contamination.R")
## mutation overlap all tissues: S9 (must precede tissue_dnv, which overwrites skyscraper_prep)
source("scripts/supp_tissue_overlap.R")
## dinucleotide variants: 4B
source("scripts/tissue_dnv.R")
## coding vs noncoding MF: 4C
source("scripts/tissue_mf_coding_noncoding.R")
## mutation signatures: 4D
source("scripts/tissue_mutsigs.R")
## lollipop: 4E
source("scripts/tissue_lollipop.R")
## dN/dS across tissue groups: 4I
source("scripts/tissue_dnds.R")
## 181 LFS mutation frequency: S7
source("scripts/supp_lfs_181_freq.R")
## phasing: 4K
## Note: an external python script must be run between these two:
##   tissue_phasing_prepare.R → phasing_tp53_181_indels.py → tissue_phasing.R
source("scripts/tissue_phasing_prepare.R")
source("scripts/tissue_phasing.R")
