# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

# ## Li-Fraumeni mutation frequencies TP53 coding and non-coding




################################################################################
############################# Coding vs non-coding #############################
################################################################################

###################################################
########### mutation frequencies
###################################################



tp53_depth_coding <- final_masked_depth %>%
  filter(gene_name == "TP53") %>%
  filter(!is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  group_by(Samp) %>% 
  summarise(denominator_coding = sum(DP))

tp53_depth_noncoding <- final_masked_depth %>%
  filter(gene_name == "TP53") %>%
  filter(is.na(exon_number)) %>%
  group_by(Samp) %>% 
  summarise(denominator_noncoding = sum(DP))

tp53_depth_split <- tp53_depth_coding %>%
  left_join(tp53_depth_noncoding)
tp53_depth_split

#### mutation frequency # of mutations/denominator
mutFreq_prep <-
  maf_masked_coding %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53") %>%
  mutate(age = age_map[Subject]) %>%
  mutate(plot_coding = if_else(!is.na(exon_number), "coding", "non-coding"))

### filter for plot
mutFreq_prep <- mutFreq_prep %>%
  left_join(tp53_depth_split) %>%
  group_by(plot_coding, Subject, denominator_coding, denominator_noncoding, age) %>%
  summarise(n_muts = n()) %>%
  mutate(mutFreq = if_else(plot_coding == "coding", n_muts/denominator_coding, n_muts/denominator_noncoding)) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx"))

offsets <- mutFreq_prep %>%
  ungroup() %>%
  group_by(Subject) %>%
  mutate(age_offset = runif(1,-2.1,2.1), 
         age_j = age + age_offset) %>%
  dplyr::select(-n_muts, -mutFreq)

mutFreq_prep <- mutFreq_prep %>%
  left_join(offsets) %>%
  mutate(
    shape_group = case_when(
      LFS == "LFS"     & CTx == "CTx"    ~ "LFS/CTx",
      LFS == "LFS"     & CTx == "non-CTx" ~ "LFS/no-CTx",
      LFS == "non-LFS" & CTx == "CTx"    ~ "non-LFS/CTx",
      LFS == "non-LFS" & CTx == "non-CTx" ~ "non-LFS/no-CTx"
    )
  )
patient_history_order <- c("LFS/no-CTx", "LFS/CTx", "non-LFS/no-CTx", "non-LFS/CTx")
mutFreq_prep$shape_group <- factor(mutFreq_prep$shape_group, levels = patient_history_order)

###############################################################################
### TP53 MF coding/non-coding ratio
###############################################################################

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS/no-CTx" = 1,
  "non-LFS/CTx"     = 16,
  "LFS/no-CTx"     = 1,
  "LFS/CTx" = 16
)

ctx_shapes <- c(
  "non-CTx" = 1,   # open circle
  "CTx"    = 16   # filled circle
)
LFS_colors <- c(
  "non-LFS" = "#44AA99",  # non-LFS color
  "LFS"     = "#882255"  # LFS color
)

MF_comparison_prep <- mutFreq_prep %>% dplyr::select(-age_offset, -age_j)
MF_ratio_prep <- MF_comparison_prep %>%
  dplyr::select(Subject, age, LFS, CTx, shape_group, mutFreq, plot_coding) %>%
  pivot_wider(
    names_from  = plot_coding,
    values_from = mutFreq
  ) %>%
  mutate(
    MF_ratio = coding / `non-coding`
  ) 
MF_ratio_prep
MF_coding_ratio <- ggplot(MF_ratio_prep,
                          aes(x = LFS, y = MF_ratio)) +
  geom_jitter(aes(color = shape_group, shape = shape_group), width = 0.15, height = 0, size = 1.5, stroke = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_y_log10(limits = c(3e-1, 31)) +
  #scale_y_log10() +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = NULL,
    y = "*TP53* MF<br>coding / non-coding"
  ) +
  theme_minimal() +
  guides(color = guide_legend(ncol = 2)) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "none",
    #legend.position   = "top",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8, margin = margin(r=0, l=1)),
    legend.margin     = margin(r=0, l=-5, b=1, t=0),
    legend.box.margin = margin(r=0, l=-35.5, b=0, t=0),
    legend.key.size   = unit(8, "pt"),
    #legend.spacing.x  = unit(-50, "pt"),
    #legend.spacing.y  = unit(3, "cm"),
    axis.text.x       = element_text(size = 8),
    axis.text.y       = element_text(size = 8),
    axis.title.y = element_markdown(size = 8)
  )
#ggsave("results/MF_coding_ratio_tp53_ms.png", MF_coding_ratio, width = 1.5, height = 2.2, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_mf_ratio.png", MF_coding_ratio, width = 1.5, height = 2.2, units = "in", dpi = 300)
