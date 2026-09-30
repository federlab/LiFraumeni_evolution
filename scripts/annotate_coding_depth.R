# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## 1. Generate a coding/non-coding file for target sequences
## For each gene sequenced, find the transcript used to call the variants
## Pull the transcript annotations from the .gff3 file
## 2. Use those annotations to call coding/non-coding in the depth file
## Then use target bed file with masking to apply final coding/non-coding
## 3. Use annotations to call coding/non-coding in mutations file

## generate list of transcripts used for variant calling

transcripts <- filt_maf %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>% 
  group_by(Hugo_Symbol, Transcript_ID) %>% 
  summarise(count = n()) %>% 
  dplyr::select(-count)
transcripts_list <- transcripts$Transcript_ID
transcripts_list <- c(transcripts_list, "ENST00000291552") # U2AF1 

## filter for annotations using transcripts in gff file
gff <- import("inputs/gencode.v38.annotation.gff3.gz")
gff <- as_tibble(gff)
annotations <- gff %>% 
  mutate(transcript_id_base = sub("\\..*", "", transcript_id)) %>%
  filter(transcript_id_base %in% transcripts_list) %>% 
  filter(type %in% c("transcript", "CDS") )%>% 
  dplyr::select(seqnames, strand, start, end, width, type, gene_name, exon_number)
write.csv(annotations, "results/CHIP_cds_annotations.csv")


############################################################################
## annotate depth file within CDS region (coding/non-coding)
############################################################################

## expand the positions of the transcript region for genic/non-genic
annotations_expand <- annotations %>% filter(type == "transcript") %>%
  dplyr::select(seqnames, start, end, gene_name) %>%
  rowwise() %>%
  mutate(Pos = list(seq(start, end))) %>%
  unnest(Pos) %>%
  dplyr::select(seqnames, Pos, gene_name) %>% 
  mutate(Pos = Pos) %>%
  distinct(seqnames, Pos, gene_name, .keep_all = TRUE)


## join genic regions
joined <- DP_filt %>%
  dplyr::rename(End = end) %>%
  filter(Samp != "DevDNA1_S1.1") %>%
  left_join(annotations_expand, by = c("Chr" = "seqnames", "End" = 'Pos'))

## add 2bp to cds regions for splice sites
annotations_splice <- annotations %>%
  filter(type == "CDS") %>%
  mutate(exon_number = as.integer(exon_number)) %>%
  group_by(gene_name) %>%
  mutate(
    first_exon = min(exon_number, na.rm = TRUE),
    last_exon  = max(exon_number, na.rm = TRUE),
    start_adj = case_when(
      strand == "+" & exon_number == first_exon ~ start,
      strand == "+" & exon_number == last_exon  ~ start - 2,
      strand == "-" & exon_number == first_exon ~ start - 2,
      strand == "-" & exon_number == last_exon  ~ start,
      TRUE ~ start - 2
    ),
    end_adj = case_when(
      strand == "+" & exon_number == first_exon ~ end + 2,
      strand == "+" & exon_number == last_exon  ~ end,
      strand == "-" & exon_number == first_exon ~ end,
      strand == "-" & exon_number == last_exon  ~ end + 2,
      TRUE ~ end + 2
    ),
    start_adj = pmax(start_adj, 0)
  ) %>%
  ungroup() 

## expand cds regions for join

cds_expand <- annotations_splice %>%
  dplyr::select(seqnames, start, end, gene_name, exon_number, start_adj, end_adj) %>%
  rowwise() %>%
  mutate(Pos = list(seq(start_adj, end_adj))) %>%
  unnest(Pos) %>%
  dplyr::select(seqnames, Pos, gene_name, exon_number) %>% 
  mutate(Pos = Pos) %>%
  distinct(seqnames, Pos, gene_name, exon_number, .keep_all = TRUE) 

# Now for CDS: only check for genic
joined_cds_depth <- joined %>%
  left_join(
    cds_expand,
    by = c("Chr" = "seqnames", "End" = "Pos", "gene_name" = "gene_name"))

###############################################################################
######### Add target bed file annotations and masking
##############################################################################

## get targets from running cds_targets_bed.R
bed_targets_path <- file.path("inputs/BEDs",
                              "coding_targets.bed")
cols <- c("chr", "start", "end")
bed_targets <- read_delim(bed_targets_path, col_names = cols)

## expand individual sites from .bed target ranges
expand_all_sites <- bed_targets %>% 
  rowwise() %>%
  mutate(Pos = list(seq(start+1, end))) %>%
  unnest(Pos) %>%
  dplyr::select(chr, Pos)


### pull repeat masking

mask_path <- file.path("inputs/BEDs", "LiFraumeni.mask.bed")
cols <- c("chr", "start", "end")
mask_file <- read_delim(mask_path, col_names = cols)
mask_file_expanded <- mask_file %>% 
  dplyr::select(chr, start, end) %>%
  rowwise() %>%
  #mutate(Pos = list(seq(start+1, end))) %>%
  mutate(Pos = list(seq(start, end))) %>%
  unnest(Pos) %>%
  dplyr::select(chr, Pos) %>% 
  distinct(chr, Pos, .keep_all = TRUE) %>% 
  mutate(inRepeatMask = TRUE)

## apply masking to bed targets
coding_sites_masked <- expand_all_sites %>%
  #left_join(mask_file_expanded) %>% 
  full_join(mask_file_expanded) %>%
  mutate(inRepeatMask = if_else(is.na(inRepeatMask), FALSE, inRepeatMask))

## apply targets/mask to depth file
final_masked_depth <- joined_cds_depth %>%
  left_join(coding_sites_masked, by = c("Chr" = "chr", "End" = "Pos"))

############################################################################
## annotate mutations within genic/CDS region/and masking
############################################################################

maf_masked_coding <- filt_maf %>%
  left_join(annotations_expand, by = c("Chromosome" = "seqnames", "Start_Position" = 'Pos')) %>%
  left_join(
    cds_expand,
    by = c("Chromosome" = "seqnames", "Start_Position" = "Pos", "gene_name" = "gene_name")) %>%
  left_join(coding_sites_masked, 
            by = c("Chromosome" = "chr", "Start_Position" = "Pos")) %>%
  filter(Tumor_Sample_Barcode != "DevDNA1_S1.1") 


testBed_MUT_add <- testBed_MUT %>%
  dplyr::rename(MUT_region = Gene)
maf_masked_coding <- maf_masked_coding %>%
  left_join(testBed_MUT_add, by = c("Chromosome" = "Chr", "Start_Position"= "Pos")) %>%
  left_join(testBed_MUT_add, by = c("Chromosome" = "Chr", "End_Position"= "Pos"),
            suffix=c("_StartPosition","_EndPosition")) 


########## apply metadata to SNV filtering table
snp_filtering <- snp_table %>% 
  left_join(annotations_expand, by = c("Chromosome" = "seqnames", "Start_Position" = 'Pos')) %>%
  left_join(
    cds_expand,
    by = c("Chromosome" = "seqnames", "Start_Position" = "Pos", "gene_name" = "gene_name")) %>%
  left_join(coding_sites_masked, 
            by = c("Chromosome" = "chr", "Start_Position" = "Pos")) %>%
  filter(Tumor_Sample_Barcode != "DevDNA1_S1.1") %>%
  filter(!inRepeatMask | is.na(inRepeatMask))

