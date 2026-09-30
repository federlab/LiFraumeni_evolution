# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Supplemental Figure S4: mutagenesis MF excluding chemotherapy-exposed subjects
## Left:  mutagenesis (non-coding MUT) mutation frequency by age
## Right: multivariate regression coefficients (MF ~ age + LFS)
##
## Requires mutFreq_prep with shape_group and age_j (chip_mf_coding_noncoding.R)
## and 00.1_SharedConstants.R (ctx_subjects, shape_group_colors, GEOM_POINT_SIZE).
## Uses a filtered copy so mutFreq_prep is unchanged for downstream scripts.

mutFreq_prep_noCTx <- mutFreq_prep %>% filter(!(Subject %in% ctx_subjects))

###############################################################################
### Mutation frequency plot
###############################################################################

shape_group_shapes <- c(
  "non-LFS/no-CTx" = 1,
  "non-LFS/CTx"     = 16,
  "LFS/no-CTx"     = 1,
  "LFS/CTx" = 16
)

plot_mutFreq_noCTx <- mutFreq_prep_noCTx %>% filter(coding == "non-coding-MUT")
lm_model_noCTx <- lm(mutFreq ~ age, data = plot_mutFreq_noCTx)
pval_noCTx <- coef(summary(lm_model_noCTx))[2, 4]

mutFreq_non_coding_noCTx <- ggplot(plot_mutFreq_noCTx, aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = plot_mutFreq_noCTx,
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  scale_y_log10(limits = c(1.3e-7, 5e-7)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(x = "Age") +
  ylab(expression(atop(NA, atop(textstyle("Mutagenesis"), textstyle("mutation frequency"))))) +
  annotate(
    "text",
    x = 60, y = 1.4e-07,
    label = paste0("italic(p) == ", signif(pval_noCTx, 2)),
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_non_coding_noCTx
ggsave("results/Manuscript_figures/Fig_S4/chip_mf_mut_only_noCTx_samples.png", mutFreq_non_coding_noCTx, width = 2, height = 1.5, units = "in", dpi = 300)

###############################################################################
### Multivariate regression coefficients (CTx term dropped)
###############################################################################

mutFreq_subject_non_coding_noCTx <- mutFreq_prep_noCTx %>%
  mutate(LFS_n = if_else(LFS == "LFS", 1, 0),
         age_decades = age / 10,
         depth_scaled = denominator / 10000000) %>%
  mutate(MF_scale = n_muts/depth_scaled) %>%
  filter(coding == "non-coding-MUT")

model_freq_non_coding_noCTx <- glm(
  MF_scale ~ 1 + age_decades + LFS_n,
  data = mutFreq_subject_non_coding_noCTx
)
summary(model_freq_non_coding_noCTx)

coef_freq_non_coding_noCTx <- tidy(model_freq_non_coding_noCTx, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))

term_labels <- c("age_decades" = "Age", "LFS_n" = "LFS")
coef_plot_non_coding_noCTx <- ggplot(coef_freq_non_coding_noCTx,
                                     aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-1.6, 1.63), breaks = seq(-1.5, 1.5, 0.5)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.2, size = 5, color = "black") +
  labs(x = expression("Effect size"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5),
        plot.margin = margin(1,1,1,1)) +
  coord_cartesian(clip = "off")

coef_plot_non_coding_noCTx
ggsave("results/Manuscript_figures/Fig_S4/chip_mf_multivar_coef_mut_only_noCTx_samples.png", coef_plot_non_coding_noCTx, width = 1.5, height = 1.5, units = "in", dpi = 300)
