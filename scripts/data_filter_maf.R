# Copyright (c) Brendan Kohrn and University of Washington

filt_maf <- MAF_table %>%
  # Apply filters
  filter(t_depth >= inputs$min_depth) %>%
  ## Per position Ns filter
  filter(NF <= inputs$maxNs) %>%
  ## Per-position VAF filter
  filter(!(VAF>= 0.3)) %>%
  mutate(path_am_order = case_when(
    am_class == "likely_pathogenic" ~ 1,
    am_class == "ambiguous" ~ 2,
    am_class == "likely_benign" ~ 3, 
    T ~ NA
  )) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  arrange(desc(t_alt_count), path_am_order) %>% 
  mutate(SampCodingOrder = rank(dplyr::desc(t_alt_count), ties.method = "first")) %>% 
  group_by(NULL)

## table of mutations greater than VAF threshold
snp_table <- MAF_table %>% 
  filter(NF <= inputs$maxNs & t_depth >= inputs$min_depth) %>%
  filter(VAF >= 0.3)
