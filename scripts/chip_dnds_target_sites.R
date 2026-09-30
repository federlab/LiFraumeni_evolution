# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Create and process a .bed file for use in VEP (all possible mutations)
## use coding_targets.bed from bed tools script
## for use in dN/dS analysis


## import cds regions from cds_targets_bed.R
bed_targets_path <- file.path("inputs/BEDs/coding_targets.bed")

cols <- c("chr", "start", "end")
bed_targets <- read_delim(bed_targets_path, col_names = cols)

## expand individual sites from .bed target ranges
expand_all_sites <- bed_targets %>% 
  rowwise() %>%
  mutate(pos = list(seq(start+1, end))) %>%
  unnest(pos) %>%
  dplyr::select(chr, pos)
expand_all_sites

## reference genome
fa_path <- file.path("inputs/ref_genomes/hg38.fa")
fa <- Rsamtools::FaFile(fa_path)

## pull the reference base for each position
ref_bases <- getSeq(
  fa,
  GRanges(
    seqnames = expand_all_sites$chr,
    ranges = IRanges(start = expand_all_sites$pos, end = expand_all_sites$pos)
  )
)
## attach reference base to the expanded dataframe
expand_all_sites$ref <- as.character(ref_bases)

## generate all possible nucleotide variants
bases <- c("A", "C", "G", "T")
vep_df_unmasked <- expand_all_sites %>%
  mutate(alt = map(ref, ~setdiff(bases, .x))) %>%
  unnest(alt) %>%
  mutate(mutation = paste0(ref, "/", alt),
         strand = "+")

## create masking file
mask_path <- file.path("inputs/BEDs/LiFraumeni.mask.bed")

cols <- c("chr", "start", "end")
mask_file <- read_delim(mask_path, col_names = cols)
mask_file_expanded <- mask_file %>% 
  dplyr::select(chr, start, end) %>%
  rowwise() %>%
  mutate(pos = list(seq(start+1, end))) %>%
  unnest(pos) %>%
  dplyr::select(chr, pos) %>% 
  distinct(chr, pos, .keep_all = TRUE) %>% 
  mutate(inRepeatMask = TRUE)

## apply masking to vep file
vep_masked <- vep_df_unmasked %>%
  left_join(mask_file_expanded) %>% 
  mutate(inRepeatMask = if_else(is.na(inRepeatMask), FALSE, inRepeatMask)) %>%
  filter(!inRepeatMask) %>% 
  mutate(pos2 = pos) %>%
  dplyr::select(chr, pos, pos2, mutation, strand)

## write .tsv for use in vep
vep_path <- file.path("results/allmutations.tsv")

write_delim(vep_masked, vep_path, delim='\t', col_names = FALSE)
