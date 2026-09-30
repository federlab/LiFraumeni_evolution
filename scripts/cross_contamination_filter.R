# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

### SNV cross-contamination analysis

ann_wide <- tibble::tibble(
  Subject    = c("UW volunteer 1","UW volunteer 2","UW volunteer 3","UW volunteer 4",
                 "UW volunteer 5","UW volunteer 6","UW volunteer 7",
                 "Patient","Family member A","Family member C","Family member B"),
  SampleCode = c("CON01","CON02","CON03","CON04","CON05","CON06","CON07",
                 "LFS01","LFS02","LFS03","REL01"),
  Age        = c(25,30,27,25,37,60,76,34,39,69,61),
  LFS        = c("No","No","No","No","No","No","No","Yes","Yes","Yes","No"),
  CTx        = c("No","No","No","No","No","No","Yes","Yes","No","No","No"),
  Depth      = c("High","High","High","High","High","High","High","High","High","High","High")
) %>%
  mutate(
    LFS   = factor(LFS, levels = c("No","Yes")),
    CTx   = factor(CTx, levels = c("No","Yes")),
    Depth = factor(Depth, levels = c("Low","High"))
  )

maf_masked_coding
unique_somatic <- maf_masked_coding %>%
  left_join(ann_wide %>% dplyr::select(Subject, SampleCode)) %>%
  mutate(SampleLabel = paste(SampleCode, Tissue, sep = ":")) %>%
  distinct(Tumor_Sample_Barcode, SampleLabel, mutPosition)

unique_snp <- snp_filtering %>%
  left_join(ann_wide %>% dplyr::select(Subject, SampleCode)) %>%
  mutate(SampleLabel = paste(SampleCode, Tissue, sep = ":")) %>%
  distinct(Tumor_Sample_Barcode, SampleLabel, mutPosition)
length(intersect(unique(unique_snp$mutPosition), unique(unique_somatic$mutPosition)))



overlap_df <- expand_grid(
  germ_sample = unique(unique_snp$Tumor_Sample_Barcode),
  som_sample  = unique(unique_somatic$Tumor_Sample_Barcode)
) %>%
  rowwise() %>%
  mutate(
    germ_muts = list(unique_snp$mutPosition[unique_snp$Tumor_Sample_Barcode == germ_sample]),
    som_muts  = list(unique_somatic$mutPosition[unique_somatic$Tumor_Sample_Barcode == som_sample]),
    germ_label = unique(unique_snp$SampleLabel[unique_snp$Tumor_Sample_Barcode == germ_sample]),
    som_label  = unique(unique_somatic$SampleLabel[unique_somatic$Tumor_Sample_Barcode == som_sample]),
    n_germ     = length(germ_muts),
    overlap    = sum(germ_muts %in% som_muts),
    proportion = ifelse(n_germ > 0, overlap / n_germ, NA_real_)
  ) %>%
  ungroup()

overlap_df
snp_overlap <- ggplot(overlap_df, aes(x = germ_label, y = som_label, fill = proportion)) +
  geom_tile(color = "grey80") +
  scale_fill_viridis_c() +
  labs(
    x = "Germline variants",
    y = "Somatic variants",
    fill = "Proportion\noverlap"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 8), 
    axis.title.y = element_text(size = 8),
    panel.grid = element_blank()
  )
snp_overlap


pair_overlap_ids <- expand_grid(
  germ_sample = unique(unique_snp$Tumor_Sample_Barcode),
  som_sample  = unique(unique_somatic$Tumor_Sample_Barcode)
) %>%
  rowwise() %>%
  mutate(
    germ_muts   = list(unique_snp$mutPosition[unique_snp$Tumor_Sample_Barcode == germ_sample]),
    som_muts    = list(unique_somatic$mutPosition[unique_somatic$Tumor_Sample_Barcode == som_sample]),
    overlap_ids = list(intersect(germ_muts, som_muts)),
    n_overlap   = length(overlap_ids)
  ) %>%
  ungroup()

overlaps_long <- pair_overlap_ids %>%
  filter(n_overlap > 0) %>%
  dplyr::select(germ_sample, som_sample, overlap_ids) %>%
  unnest(overlap_ids) %>%
  dplyr::rename(mutPosition = overlap_ids) %>%
  # attach labels for both sides
  left_join(unique_snp %>% dplyr::select(Tumor_Sample_Barcode, SampleLabel, mutPosition),
            by = c("germ_sample" = "Tumor_Sample_Barcode", "mutPosition" = "mutPosition")) %>%
  dplyr::rename(germ_label = SampleLabel) %>%
  left_join(unique_somatic %>% dplyr::select(Tumor_Sample_Barcode, SampleLabel, mutPosition),
            by = c("som_sample"  = "Tumor_Sample_Barcode", "mutPosition" = "mutPosition")) %>%
  dplyr::rename(som_label = SampleLabel) %>%
  arrange(germ_label, som_label, mutPosition)

overlap_exclude <- overlaps_long %>% 
  distinct(mutPosition) %>% 
  mutate(SNP_exclude = TRUE)

## exclude snp mutations from all samples 
maf_masked_coding <- maf_masked_coding %>% 
  left_join(overlap_exclude) %>% 
  mutate(SNP_exclude = if_else(is.na(SNP_exclude), FALSE, SNP_exclude)) %>%
  filter(!SNP_exclude)


nm8_flagged <- maf_masked_coding %>%
  filter(str_detect(FILTER, "NM8\\.0")) %>% 
  group_by(mutPosition) %>%  
  filter(n_distinct(Samp) >= 3) %>% 
  ungroup()

removed_df <- maf_masked_coding %>%
  semi_join(
    nm8_flagged %>% dplyr::select(mutPosition, Samp),
    by = c("mutPosition", "Samp")
  )

maf_masked_coding <- maf_masked_coding %>%
  anti_join(
    nm8_flagged %>% dplyr::select(mutPosition, Samp),
    by = c("mutPosition", "Samp")
  )

