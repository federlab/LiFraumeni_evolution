# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington


###############################################################################
### Mutation burden
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

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS/no-CTx" = 1,
  "non-LFS/CTx"     = 16,
  "LFS/no-CTx"     = 1,
  "LFS/CTx" = 16
)

## coding model

lm_model <- lm(mutBurden ~ age, data = mutFreq_prep %>% filter(coding == "coding"))
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutBurd_coding <- ggplot(mutFreq_prep %>% filter(coding == "coding"), aes(x = age_j, y = mutBurden)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'coding'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_y_log10(limits = c(1e-7, 1e-5)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Coding\nmutation burden"
  ) +
  annotate(
    "text",
    x = 60, y = 1.4e-07,
    label = paste0("italic(p) == ", signif(pval, 2)),,
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal() 

mutBurd_coding


################################################################################
########### 
################################################################################

mutBurd_coding <- mutBurd_coding + theme(legend.position = "none")

legend_shared <- get_legend(
  mutBurd_coding + 
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


mutBurd_coding <- mutBurd_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  ylab(expression(atop(NA, atop(textstyle("CHIP"), textstyle("mutation burden"))))) +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))


combined_plots <- plot_grid(mutBurd_coding, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))



################################################################################
########### 
################################################################################

mutFreq_prep_CHIP_encode <- mutFreq_prep %>%
  mutate(
    LFS_n = if_else(LFS == "LFS", 1, 0),
    CTx_n = if_else(CTx == "CTx", 1, 0)
  )
mutFreq_subject <- mutFreq_prep_CHIP_encode %>%
  mutate(age_decades = age / 10,
         depth_scaled = denominator / 10000000) %>%
  mutate(MB_scale = mutReads/depth_scaled)

## coding
mutFreq_subject_coding <- mutFreq_subject %>%
  filter(coding == "coding")

model_freq_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MB_scale ~ age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_coding
)

coef_freq_CHIP_coding <- tidy(model_freq_coding, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))


term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
CHIP_burden_plot_coding <- ggplot(coef_freq_CHIP_coding,
                                aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-25, 45), breaks = seq(-30, 50, 10)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.2, size = 5, color = "black") +
  labs(x = expression("Effect size (mutant bases / "  * 10^7 * " bases)"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

combined_plot <- combined_plots + CHIP_burden_plot_coding +
  plot_layout(widths = c(2,1))

combined_plot <- combined_plot + theme(plot.margin = margin(0,0,10,0))

legend_shared + theme(plot.margin = margin(0,0,0,0))
combined_plot <- plot_grid(mutBurd_coding, CHIP_burden_plot_coding, legend_shared, ncol = 2, 
                            rel_heights = c(1.1,0.15))

ggsave(paste0("results/Manuscript_figures/Fig_S3/supp_chip_mutation_burden.png"), combined_plot, width = 6, height = 3.5, units = "in", dpi = 300)


## save supp table
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
      model  = "MB ~ 1 + decades + LFS + CTx",
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
row_coding <- model_to_row(model_freq_coding,      gene = "CHIP", coding = "coding")

write.table(
  row_coding,
  file = "results/MB_multivariate_CHIP_summary.csv",
  sep = ",",
  row.names = FALSE
)



