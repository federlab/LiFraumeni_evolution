# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## combine all alphamissense scores with all genes

am_all <- vroom("inputs/alphamissense/AlphaMissense_hg38.tsv.gz", delim="\t", skip = 3, progress = TRUE) %>%
  dplyr::rename(Chromosome = "#CHROM",
                Start_Position = POS,
                Reference_Allele = REF,
                Tumor_Seq_Allele2 = ALT)

joined <- maf_masked_coding %>%
  left_join(am_all %>% 
              dplyr::select(Chromosome, Start_Position, Reference_Allele,
                     Tumor_Seq_Allele2, am_pathogenicity),
            by = c("Chromosome", "Start_Position",
                   "Reference_Allele", "Tumor_Seq_Allele2"))

write_delim(joined, "results/all_genes_am.tsv")
