# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

### Output a Supplemental table with gene/exon and region sequenced and proportion
### of the region sequenced

columns <- c("Chromosome", 
             "Start", 
             "End", 
             "Gene", 
             "Exon", 
             "Strand",
             "X",
             "Sequenced_base_pairs", 
             "Total_base_pairs", 
             "Proportion_sequenced")

coverage_df <- read_delim("inputs/BEDs/cds_baits_overlap.bed", col_names = columns)

cds_annotations <- read_delim("results/CHIP_cds_annotations.csv")

cds_annotations <- cds_annotations %>%
  filter(type == "CDS") %>%
  dplyr::select(seqnames, start, end, gene_name, exon_number)
  

joined <- cds_annotations %>%
  mutate(join_start = start -1) %>%
  left_join(coverage_df, by = c("seqnames" = "Chromosome",
                                       "join_start" = "Start",
                                       "end" = "End")) %>%


  mutate(
    Chromosome = seqnames,
    Start = join_start,
    End = end,
    Transcript = Gene,
    Gene = gene_name,
    Exon = exon_number
  ) %>%
  dplyr::select(Chromosome, Start, End, Gene, Exon, Transcript, Sequenced_base_pairs, Total_base_pairs, Proportion_sequenced) %>%
  filter(Sequenced_base_pairs >0) %>%
  mutate(Transcript = sub("_.*", "", Transcript))

out_file <- "results/exon_proportions.tsv"
write_tsv(
  joined,
  out_file,
  col_names = TRUE
)

