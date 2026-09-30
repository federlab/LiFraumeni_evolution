# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

#### coding non coding in Tissues


sample_id_mapping_path <- "inputs/sampleID_mapping.txt"

sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))


###################################################
########### Pull sequencing depths in coding regions
###################################################

GEOM_POINT_SIZE = 1.5

## by gene
chip_depths <- final_masked_depth %>% 
  filter(!is.na(exon_number)) %>%
  filter(gene_name == "TP53") %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(subject %in% c("Patient"))

chip_depths_non_coding <- final_masked_depth %>%
  filter(!is.na(gene_name)) %>%
  filter(gene_name == "TP53") %>%
  filter(is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(subject %in% c("Patient"))

## mutFreq by gene
mutFreq_counts <- 
  maf_masked_coding %>%
  filter(gene_name == "TP53") %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(Subject %in% c("Patient")) %>%
  #filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol, age) %>% 
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )


mutFreq_combined <- chip_depths %>%
  left_join(mutFreq_counts, 
            by = c("Samp" = "Tumor_Sample_Barcode",
                   "gene_name" = "Hugo_Symbol")) %>%
  mutate(
    n_muts   = replace_na(n_muts, 0),
    mutReads = replace_na(mutReads, 0),
    mutFreq   = n_muts   / denominator,
    mutBurden = mutReads / denominator,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
  dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

#### write non-coding
mutFreq_counts_non_coding <- 
  maf_masked_coding %>%
  filter(!is.na(gene_name)) %>%
  filter(gene_name == "TP53") %>%
  filter(is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  filter(Subject %in% c("Patient")) %>%
  #filter(Tissue %in% family_patient_blood_samples) %>%
  filter(!Variant_Classification == "Splice_Site") %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol, age) %>%
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )



mutFreq_combined_non_coding <- chip_depths_non_coding %>%
  left_join(mutFreq_counts_non_coding, 
            by = c("Samp" = "Tumor_Sample_Barcode",
                   "gene_name" = "Hugo_Symbol")) %>%
  mutate(
    n_muts = replace_na(n_muts, 0),
    mutReads = replace_na(mutReads, 0),
    mutFreq  = n_muts   / denominator,
    mutBurden = mutReads / denominator,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
  dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

## write mutation frequencies by gene/subject

mutFreq_combined_coding <- mutFreq_combined %>% mutate(coding = "coding")
mutFreq_combined_non_coding <- mutFreq_combined_non_coding %>% mutate(coding = "non-coding-CHIP")
MF_CHIP_genic <- rbind(mutFreq_combined_coding, mutFreq_combined_non_coding)


########## Add mutagenesis region depths and mutation counts to non_coding genes
## first aggregate genes 
mutFreq_prep <- MF_CHIP_genic %>% 
  group_by(Subject, tissue, age, LFS, CTx, coding) %>%
  summarise(denominator = sum(denominator), 
            n_muts = sum(n_muts),
            mutReads = sum(mutReads)) %>%
  mutate(mutFreq = n_muts / denominator,
         mutBurden = mutReads / denominator)

## add mutagenesis regions depth
mutagenesis_depths <- final_masked_depth %>%
  filter(str_starts(Gene, "region")) %>% 
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  #filter(tissue %in% family_patient_blood_samples) %>% 
  mutate(coding = "non-coding-MUT")

## find mutagenesis mutation counts
mutFreq_counts_non_coding_mutagenesis <- 
  maf_masked_coding %>%
  filter(!is.na(MUT_region_StartPosition) | !is.na(MUT_region_EndPosition)) %>%
  filter(Subject %in% c("Patient")) %>%
  #filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  group_by(Subject, Tissue) %>%
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )

MF_MUT <- mutFreq_counts_non_coding_mutagenesis %>%
  left_join(mutagenesis_depths, by = c("Subject" = "subject", "Tissue" = "tissue"))

MF_MUT

apply_data <- mutFreq_prep %>% filter(coding == 'coding') %>% dplyr::select(Subject, tissue, age, LFS, CTx)
MF_MUT <- MF_MUT %>% left_join(apply_data, by = c("Tissue" = "tissue", 'Subject')) %>%
  mutate(mutFreq = n_muts/denominator,
         mutBurden = mutReads/denominator) %>%
  dplyr::select(-Samp) %>%
  dplyr::rename(tissue = Tissue)

mutFreq_prep <- rbind(mutFreq_prep, MF_MUT)
mutFreq_prep
## calculate total non-coding (CHIP + MUT)
total_non_coding <- mutFreq_prep %>% filter(coding != "coding") %>%
  group_by(Subject, tissue, age, LFS, CTx) %>%
  mutate(denominator_total = sum(denominator), 
         n_muts_total = sum(n_muts), 
         mutReads_total = sum(mutReads),
         mutFreq_total = n_muts_total/denominator_total,
         mutBurden_total = mutReads_total/denominator_total) %>%
  dplyr::select(-denominator, -n_muts, -mutReads, -mutFreq, -mutBurden) %>%
  dplyr::rename(denominator = denominator_total, n_muts = n_muts_total, 
                mutReads = mutReads_total, mutFreq = mutFreq_total, 
                mutBurden = mutBurden_total) %>%
  mutate(coding = "non-coding-total") %>%
  distinct()

total_non_coding
mutFreq_prep <- rbind(mutFreq_prep, total_non_coding)


###############################################################################
### Mutation frequency
###############################################################################

MF_tissue_groups <- mutFreq_prep %>% 
  left_join(tissue_categories, by=c("tissue" = "Tissue"))

MF_tissue_groups_tp53 <- MF_tissue_groups %>% filter(coding == "coding" | coding == "non-coding-CHIP")

## save supp table
MF_tissue_groups_tp53_table <- MF_tissue_groups_tp53 %>%
  dplyr::rename(Tissue = tissue,
                Age = age,
                MF = mutFreq) %>% 
  dplyr::select(-mutReads, -mutBurden,-Category) %>% 
  mutate(coding = if_else(coding == "non-coding-CHIP", "non-coding", coding))
write_delim(MF_tissue_groups_tp53_table, "results/MF_tissues_table.csv", delim=',')


MF_ratio_tissue_groups <- MF_tissue_groups %>% filter(coding == "coding" | coding == "non-coding-CHIP") %>%
  group_by(Category, coding) %>%
  summarise(denominator = sum(denominator), 
         n_muts = sum(n_muts) 
         ) %>% 
  mutate(MF = n_muts/denominator) %>% 
  dplyr::select(Category, coding, MF) %>%
  tidyr::pivot_wider(
    names_from = coding,
    values_from = MF
  ) %>%
  mutate(MF_ratio = coding / `non-coding-CHIP`) 
MF_ratio_tissue_groups

write_delim(MF_ratio_tissue_groups, "results/MF_ratio_tissue_groups.csv", delim=',')

MF_ratio_tissues <- MF_tissue_groups %>% filter(coding == "coding" | coding == "non-coding-CHIP") %>%
  group_by(tissue, coding) %>%
  summarise(denominator = sum(denominator), 
            n_muts = sum(n_muts) 
  ) %>% 
  mutate(MF = n_muts/denominator) %>% 
  dplyr::select(tissue, coding, MF) %>%
  tidyr::pivot_wider(
    names_from = coding,
    values_from = MF
  ) %>%
  mutate(MF_ratio = coding / `non-coding-CHIP`) 
MF_ratio_tissues

MF_tissue_groups


## from tissue_skyscraper.R
x_order_saved <- levels(skyscraper_prep$Tissue_ordered)




MF_tissue_groups <- MF_tissue_groups %>%
  mutate(tissue = factor(tissue, levels = x_order_saved))




tissue_coding_non_MF <- ggplot(MF_tissue_groups %>% filter(coding == "coding" | coding == "non-coding-CHIP"), aes(x = tissue, y = mutFreq, fill = coding)) +
  geom_col(position = position_dodge(width = 1)) +
  #scale_y_log10() +
  scale_fill_manual(
    values = c("coding" = "black",
               "non-coding-CHIP" = "#999999"),
    labels = c("coding" = "coding", "non-coding-CHIP" = "non-coding")
  ) +
  labs(
    x = "Tissue",
    y = "MF"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size =8),
    axis.title.y = element_text(size=8),
    axis.text.y = element_text(size=8),
    axis.title.x = element_blank()
  )


MF_base <- MF_tissue_groups %>%
  filter(coding %in% c("coding", "non-coding-CHIP"))

custom_label <- function(x) {
  sapply(x, function(t) {
    abbr <- tissue_abbreviations$Tissue_abbr[match(t, tissue_abbreviations$Tissue)]
    if (is.na(abbr)) abbr <- t
    if (t %in% cancer_samples) {
      paste0("<span style='color:red;'>", abbr, "</span>")
    } else {
      abbr
    }
  })
}


tissues_to_exclude <- c("Skin, non-sun-exposed", "Thyroid", "Mainstem bronchus", "Colon")
MF_tissue_groups_TP53_exclude_low_depth <- MF_tissue_groups_tp53_table %>%
  filter(!(Tissue %in% tissues_to_exclude)) %>%
  left_join(tissue_categories) %>%
  mutate(Category = if_else(Tissue == "Skin", "Solid", Category)) %>%
  mutate(Category = factor(Category, levels = c("Blood", "Solid", "Cancer"))) %>%
  group_by(Category, coding) %>%
  arrange(MF) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)))

MF_tissue_groups_TP53_exclude_low_depth

TP53_MF_tissues <- ggplot(data = MF_tissue_groups_TP53_exclude_low_depth, aes(x = Tissue, y = MF, fill = coding, color = coding)) +
  #geom_point(size = 2, shape = 1) +
  geom_point(size = 2, shape = 21, alpha = 0.75) +
  theme_minimal() +
  scale_fill_manual(
    values = c("coding" = "#999999",
               "non-coding" = "#000000")
  ) +
  scale_color_manual(
    values = c("coding" = "#999999",
               "non-coding" = "#000000")
  ) +
  scale_y_continuous(
    labels = function(x) ifelse(x == 0, "0", scales::scientific(x))
  ) +
  scale_x_discrete(labels = custom_label) +
  labs(y = "*TP53* MF") +
  theme(axis.title.x = element_blank(),
        axis.text.x.bottom = element_markdown(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_markdown(size = 8),
        #legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 8, margin = margin(0,0,0,1)), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(2,0,5,45),
        legend.key.spacing.x = unit(3, "pt"),
        legend.position = c(0.3,0.8))

TP53_MF_tissues

ggsave("results/Manuscript_figures/Fig_4/tissue_mf_coding_noncoding.png", TP53_MF_tissues, width = 3.75, height = 2, units = "in", dpi = 300)


