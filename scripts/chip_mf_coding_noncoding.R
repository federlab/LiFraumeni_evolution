# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

# ## Li-Fraumeni mutation frequencies all chip genes


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
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)

chip_depths_non_coding <- final_masked_depth %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>% 
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)

## mutFreq by gene
mutFreq_counts <- 
  maf_masked_coding %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
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
    age = age_map[subject]) %>%
    dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

#### write non-coding
mutFreq_counts_non_coding <- 
  maf_masked_coding %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number)) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(!Variant_Classification == "Splice_Site") %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>% 
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
    age = age_map[subject]) %>%
  dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")


## write mutation frequencies by gene/subject
mutFreq_combined_coding <- mutFreq_combined %>% mutate(coding = "coding")
mutFreq_combined_non_coding <- mutFreq_combined_non_coding %>% mutate(coding = "non-coding-CHIP")
MF_CHIP_genic <- rbind(mutFreq_combined_coding, mutFreq_combined_non_coding)
write.csv(MF_CHIP_genic, "results/CHIP_muts_by_gene.csv")


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
  filter(tissue %in% family_patient_blood_samples) %>% 
  mutate(coding = "non-coding-MUT")

## find mutagenesis mutation counts
mutFreq_counts_non_coding_mutagenesis <- 
  maf_masked_coding %>%
  filter(!is.na(MUT_region_StartPosition) | !is.na(MUT_region_EndPosition)) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  group_by(Subject, Tissue) %>%
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )

MF_MUT <- mutFreq_counts_non_coding_mutagenesis %>%
  left_join(mutagenesis_depths, by = c("Subject" = "subject", "Tissue" = "tissue"))

apply_data <- mutFreq_prep %>% filter(coding == 'coding') %>% dplyr::select(Subject, tissue, age, LFS, CTx)
MF_MUT <- MF_MUT %>% left_join(apply_data, by = c("Tissue" = "tissue", 'Subject')) %>%
  mutate(mutFreq = n_muts/denominator,
         mutBurden = mutReads/denominator) %>%
  dplyr::select(-Samp) %>%
  dplyr::rename(tissue = Tissue)

mutFreq_prep <- rbind(mutFreq_prep, MF_MUT)
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

mutFreq_prep <- rbind(mutFreq_prep, total_non_coding)


###############################################################################
### Mutation frequency
###############################################################################


## offsets are used for plotting manual jitter of individuals with same age
offsets <- mutFreq_prep %>%
  ungroup() %>%
  group_by(Subject) %>%
  mutate(age_offset = runif(1,-2.1,2.1), 
         age_j = age + age_offset)

mutFreq_prep <- offsets %>%
  mutate(
    shape_group = case_when(
      LFS == "LFS" & CTx == "CTx" ~ "LFS/CTx",   # both LFS and CTx
      LFS == "LFS"             ~ "LFS/no-CTx",  # just LFS
      LFS == "non-LFS" & CTx == "CTx" ~ "non-LFS/CTx",
      TRUE                 ~ "non-LFS/no-CTx"      # if you want a default
    )
  )

patient_history_order = c("LFS/no-CTx", "LFS/CTx", "non-LFS/no-CTx", "non-LFS/CTx")
mutFreq_prep$shape_group = factor(mutFreq_prep$shape_group, levels = patient_history_order)

###############################################################################
### Mutation frequency plot
###############################################################################

## function to turn "e-6" to "10^-6"
fancy_scientific <- function(l) { 
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text = l)
}

LFS_color <- "#882255"
nonLFS_color <- "#44aa99"

# Define colors by LFS status

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS/no-CTx" = 1,
  "non-LFS/CTx"     = 16,
  "LFS/no-CTx"     = 1,
  "LFS/CTx" = 16
)

## coding model

lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "coding"))
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutFreq_coding <- ggplot(mutFreq_prep %>% filter(coding == "coding"), aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'coding'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_y_log10(limits = c(1.3e-7, 4e-7)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Coding\nmutation frequency"
  ) +
  annotate(
    "text",
    x = 60, y = 1.4e-07,
    label = paste0("italic(p) == ", signif(pval, 2)),,
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal() 

mutFreq_coding

###############################################################################
### Non-coding plot
###############################################################################

# Define colors by LFS status

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS/no-CTx" = 1,
  "non-LFS/CTx"     = 16,
  "LFS/no-CTx"     = 1,
  "LFS/CTx" = 16
)

#non_coding_set = 0 # non-coding-total
non_coding_set = 1 # non-coding-MUT
#non_coding_set = 2 # non-coding-CHIP

if(non_coding_set == 0){
  lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "non-coding-total"))
  plot_mutFreq = mutFreq_prep %>% filter(coding == "non-coding-total")
  ylimits = c(1.3e-7, 4.2e-7)
  ylabel = "Non-coding\nCHIP + MUT MF"
  file_out = "chip_mf_total_noncoding.png"
}
if(non_coding_set == 1){
  lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "non-coding-MUT"))
  plot_mutFreq = mutFreq_prep %>% filter(coding == "non-coding-MUT")
  ylimits = c(1.3e-7, 5e-7)
  ylabel = "Mutagenesis\nmutation frequency"
  file_out = "chip_mf_mut_only.png"
}
if(non_coding_set == 2){
  lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "non-coding-CHIP"))
  plot_mutFreq = mutFreq_prep %>% filter(coding == "non-coding-CHIP")
  ylimits = c(1e-7, 3e-7)
  ylabel = "non-coding\nCHIP MF"
  file_out = "chip_mf_chip_only.png"
}

model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

plot_mutFreq

mutFreq_non_coding <- ggplot(plot_mutFreq, aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = plot_mutFreq,
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  #scale_y_log10(limits = c(1.5e-6, 3e-5), labels=fancy_scientific) +
  scale_y_log10(limits = ylimits) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = ylabel
  ) +
  annotate(
    "text",
    x = 60, y = 1.4e-07,
    label = paste0("italic(p) == ", signif(pval, 2)),,
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal()


################################################################################
########### 
################################################################################

mutFreq_coding <- mutFreq_coding + theme(legend.position = "none")
mutFreq_non_coding <- mutFreq_non_coding + theme(legend.position = "none")

legend_shared <- get_legend(
  mutFreq_coding + 
    guides(
      color = guide_legend(nrow = 2, title.position = "top"),  # force 1 row
      shape = guide_legend(nrow = 2, title.position = "top")   # force 1 row
    ) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text  = element_text(size = 8),
          legend.margin = margin(r=1, l=1, b = -20, t = -25),
          legend.key.size   = unit(0.5, "lines"),
          legend.spacing.x  = unit(0.2, "cm"),
          legend.spacing.y  = unit(0.2, "cm")
    )
)

mutFreq_non_coding <- mutFreq_non_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  ylab(expression(atop(NA, atop(textstyle("Mutagenesis"), textstyle("mutation frequency"))))) +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_coding <- mutFreq_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  ylab(expression(atop(NA, atop(textstyle("CHIP"), textstyle("mutation frequency"))))) +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

plots <- plot_grid(
  plot_grid(mutFreq_coding, mutFreq_non_coding, ncol = 2, align = "v"), 
  ncol = 1)

combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))

ggsave(paste0("results/Manuscript_figures/Fig_2/", file_out), combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)
