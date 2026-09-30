# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

### multiple regression analysis TP53 only

## generate mutFreq_combined from chip_mf_coding_noncoding.R

# helper to encode and fit
mutFreq_prep_CHIP_encode <- mutFreq_prep %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  )

mutFreq_prep_CHIP_encode
##########################################
# CHIP Frequency
##########################################


mutFreq_subject <- mutFreq_prep_CHIP_encode %>%
  mutate(age_decades = age / 10,
         depth_scaled = denominator / 10000000) %>%
  mutate(MF_scale = n_muts/depth_scaled)


## coding
mutFreq_subject_coding <- mutFreq_subject %>%
  filter(coding == "coding")
mutFreq_subject_non_coding <- mutFreq_subject %>%
  filter(coding == "non-coding-CHIP")

mutFreq_subject_non_coding

model_freq_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_coding
)
summary(model_freq_coding)

model_freq_non_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ 1 + age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_non_coding
)
summary(model_freq_non_coding)
mutFreq_subject_non_coding

coef_freq_CHIP_coding <- tidy(model_freq_coding, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))
coef_freq_CHIP_coding
term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
CHIP_freq_plot_coding <- ggplot(coef_freq_CHIP_coding,
                                aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-3, 9), breaks = seq(-3, 9, 3)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  #labs(x = expression("Effect size\n(mutations / "  * 10^7 * " bases)"), y = NULL) +
  labs(x = expression("Effect size"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

#ggsave("results/MF_tp53_multi_coding_ms.png", CHIP_freq_plot_coding, width = 1.5, height = 1.5, units = "in", dpi = 300)

#########
##non-coding

coef_freq_CHIP_non_coding <- tidy(model_freq_non_coding, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))

term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")

CHIP_freq_plot_non_coding <- ggplot(coef_freq_CHIP_non_coding,
                                    aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-2.2, 2.5), breaks = seq(-2, 2, 1)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
  #labs(x = expression("Effect size\n(mutations / "  * 10^7 * " bases)"), y = NULL) +
  labs(x = expression("Effect size"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

#ggsave("results/MF_tp53_multi_coding_ms.png", CHIP_freq_plot, width = 1.5, height = 1.5, units = "in", dpi = 300)

####### freq + burden plot

CHIP_freq_plot_coding <- CHIP_freq_plot_coding + 
  theme(plot.margin = margin(1,1,1,1)) + 
  coord_cartesian(clip = "off")

CHIP_freq_plot_non_coding <- CHIP_freq_plot_non_coding + 
  theme(plot.margin = margin(1,1,1,1)) +
  coord_cartesian(clip = "off")

combined_plot <- CHIP_freq_plot_coding | CHIP_freq_plot_non_coding 

#ggsave("results/tp53_mf_multivar_coef.png", combined_plot, width = 3.5, height = 1, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_mf_multivar_coef.png", combined_plot, width = 3.5, height = 1, units = "in", dpi = 300)

# save supp table

model_to_row <- function(model, gene, coding) {
  
  coding_tidy <- broom::tidy(model)
  
  coding_wide <- coding_tidy %>%
    dplyr::select(term, estimate, std.error, p.value) %>%
    dplyr::mutate(
      est_se = sprintf("%.4g (%.4g)", estimate, std.error),
      p      = p.value
    ) %>% 
    dplyr::select(term, est_se, p) %>% 
    dplyr::filter(term %in% c("(Intercept)", "age_decades", "LFS_n", "CTx_n")) %>% 
    tidyr::pivot_wider(
      names_from  = term,
      values_from = c(est_se, p),
      names_glue  = "{term}.{.value}"
    ) %>% 
    dplyr::mutate(
      model  = "MF ~ 1 + decades + LFS + CTx",
      gene   = gene,
      coding = coding,
      .before = 1
    ) %>% 
    dplyr::rename(
      `(Intercept).Estimate` = `(Intercept).est_se`,
      `decades.Estimate`     = `age_decades.est_se`,
      `LFS.Estimate`         = `LFS_n.est_se`,
      `CTx.Estimate`         = `CTx_n.est_se`,
      `decades.p`            = `age_decades.p`,
      `LFS.p`                = `LFS_n.p`,
      `CTx.p`                = `CTx_n.p`
    ) %>%
    dplyr::select(
      model, gene, coding,
      `(Intercept).Estimate`, `(Intercept).p`,
      `decades.Estimate`,     `decades.p`,
      `LFS.Estimate`,         `LFS.p`,
      `CTx.Estimate`,         `CTx.p`
    )
  
  coding_wide
}
row_coding <- model_to_row(model_freq_coding,      gene = "TP53", coding = "coding")
row_noncoding <- model_to_row(model_freq_non_coding,   gene = "TP53", coding = "non-coding")
model_table <- bind_rows(row_coding, row_noncoding)

write.table(
  model_table,
  file = "results/MF_multivariate_TP53_summary.csv",
  sep = ",",
  row.names = FALSE
)

