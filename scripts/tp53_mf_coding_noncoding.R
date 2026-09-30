# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

# ## Li-Fraumeni mutation frequencies tp53 coding/non-coding



sample_id_mapping_path <- "inputs/sampleID_mapping.txt"
sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))


###################################################
########### Pull sequencing depths in coding regions
###################################################

GEOM_POINT_SIZE = 1.5


### generate from chip_mf_coding_noncoding.R
MF_CHIP_genic

###############################################################################
### Mutation frequency
###############################################################################


mutFreq_prep <- MF_CHIP_genic %>% filter(Hugo_Symbol == "TP53")

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
mutFreq_prep
lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "coding"))
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutFreq_coding <- ggplot(mutFreq_prep %>% filter(coding == "coding"), aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'coding'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_y_log10(limits = c(1e-7, 1.2e-6)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Coding\nmutation frequency"
  ) +
  annotate(
    "text",
    x = 60, y = 4.5e-07,
    label = paste0("italic(p) == ", signif(pval, 3)),,
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal() 

mutFreq_coding
###############################################################################
### Non-coding plot
###############################################################################

lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "non-coding-CHIP"))
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutFreq_prep %>% filter(coding == "non-coding-CHIP")
mutFreq_non_coding <- ggplot(mutFreq_prep %>% filter(coding == "non-coding-CHIP"), aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'non-coding-CHIP'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  #scale_y_log10(limits = c(1.5e-6, 3e-5), labels=fancy_scientific) +
  scale_y_log10(limits = c(3.29e-8, 4.2e-7)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Non-coding\nmutation frequency"
  ) +
  annotate(
    "text",
    x = 60, y = 5e-08,
    label = paste0("italic(p) == ", signif(pval, 3)),,
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal()

mutFreq_non_coding
################################################################################
########### Plot frequency and Burden combined
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
  #ylab(expression(atop(NA, atop(italic(TP53) * " non-coding","mutation frequency")))) +
  labs(y = "*TP53* non-coding<br>mutation frequency") +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        #axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
        axis.title.y = element_markdown(size = 8, margin= margin(l=0), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_coding <- mutFreq_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  #ylab(expression(atop(NA, atop(italic(TP53) * " coding","mutation frequency")))) +
  labs(y = "*TP53* coding<br>mutation frequency") +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        #axis.title.y = element_markdown(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.title.y = element_markdown(size = 8, margin= margin(l=0), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_coding
plots <- plot_grid(
  plot_grid(mutFreq_coding, mutFreq_non_coding, ncol = 2, align = "v"), 
  ncol = 1)

combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))

ggsave("results/Manuscript_figures/Fig_3/tp53_mf_coding_noncoding.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)

