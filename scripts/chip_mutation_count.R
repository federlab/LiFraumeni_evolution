# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Figure 2A
## Mutation counts by gene in coding regions of the CHIP panel
## generate x-axis labels in tp53_skyscraper_blood.R first

## devDNA and non-buffy coat blood samples for comparison
samples_all_exclude <- c("DNA03980_S17.1", "DNA03965_S30_S31.1", "DNA03966_S32.1", "DevDNA1_S1.1")

#######################################################
##### Number of mutations per subject
#######################################################

mutationCounts <- 
  maf_masked_coding %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>% # only chip panel
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>% # exclude non-coding 
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>% #exclude repeat masking
  group_by(Subject, Hugo_Symbol) %>%
  summarise(total_mut_count = n(), .groups = "drop") %>%
  mutate(Subject = factor(
  Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
  ))

## ann_wide from tp53_skyscraper_blood.R
mutationCounts <- mutationCounts %>%
  left_join(ann_wide) %>%
  ungroup() %>%
  mutate(SampleCode = fct_reorder(SampleCode, Age, .fun = min))


gene_colors = c("#800000", 
                "#9A6324", 
                "#808000", 
                "#469990", 
                "#000075", 
                "#e6194B", 
                "#f58231", 
                "#ffe119", 
                "#bfef45", 
                "#3cb44b", 
                "#42d4f4", 
                "#4363d8", 
                "#911eb4", 
                "#f032e6", 
                "#a9a9a9", 
                "#fabed4", 
                "#ffd8b1", 
                "#fffac8", 
                "#aaffc3", 
                "#dcbeff", 
                "#000000", 
                "#DDDDDD", 
                "#AA8977", 
                "#446633", 
                "#b0e0e6", 
                "#ffb6c1", 
                "#dda0dd", 
                "#add8e6"  
                            )

## reorder by total mutations across subjects
gene_order <- mutationCounts %>%
  group_by(Hugo_Symbol) %>%
  summarise(total_all = sum(total_mut_count, na.rm = TRUE)) %>%
  arrange(desc(total_all)) %>%
  pull(Hugo_Symbol)

mutationCount_plot <- mutationCounts %>% 
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = gene_order)) %>%
  ggplot(aes(x = SampleCode,y = total_mut_count, fill = Hugo_Symbol)) + 
  geom_col() +
  theme_minimal()  +
  scale_fill_manual(values = gene_colors) +
  scale_y_continuous(limits = c(0,200), expand =c(0,0)) +
  labs(y = "Number of CHIP mutations", fill = "Gene") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text  = element_text(size = 8),
        axis.title.y = element_markdown(size=8),
        legend.margin = margin(l=-8),
        legend.key.size = unit(8, "pt"),
        panel.border = element_rect(color = "black"),
        plot.margin = margin(t=0, r=0, l=0,b=4)) +
  guides(fill = guide_legend(ncol = 2))
show(mutationCount_plot)


## generate from tp53_skyscraper_blood.R

p_age <- p_age + theme(
  legend.margin = margin(l=-8, t=5, r=0, b=0),
  legend.key.spacing = unit(0.7, "pt"),
  legend.title = element_text(margin=margin(b=2)),
  plot.margin = margin(b=0))

p_depth <- p_depth + theme(
  legend.margin = margin(l=1, t=5, r=0, b=0),
  legend.key.spacing = unit(0.7, "pt"),
  legend.title = element_text(margin=margin(b=2)),
  plot.margin = margin(b=0))

p_lfs <- p_lfs + theme(
  legend.margin = margin(l=1, t=5, r=0, b=0),
  legend.key.spacing = unit(0.7, "pt"),
  legend.title = element_text(margin=margin(b=2)),
  plot.margin = margin(b=0))

# Collect legends for the annotation rows (on the right)
annot_block <- (p_age / p_depth / p_lfs) +
  plot_layout(heights = c(0.1, 0.1, 0.1), guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "horizontal",
    legend.justification = c(0.5, 0.5),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.spacing.x = unit(0.1, "lines"),
    axis.text.y = element_text(margin=margin(0,0,0,0)),
    plot.margin = margin(b=0)
  )

combined_plot_legend <- (mutationCount_plot / annot_block) +
  plot_layout(heights = c(5,0.75)) + 
  theme(plot.margin = margin(b=0))

combined_plot_legend



mutationCount_plot_crop <- mutationCount_plot + theme(axis.title.y = element_blank())
p_age_crop <- p_age + theme(axis.title.y = element_blank())
p_depth_crop <- p_depth + theme(axis.title.y = element_blank())
p_lfs_crop <- p_lfs + theme(axis.title.y = element_blank())
  
annot_block_crop <- (p_age_crop / p_depth_crop / p_lfs_crop) +
  plot_layout(heights = c(0.1, 0.1, 0.1), guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "horizontal",
    legend.justification = c(0.5, 0.5),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.spacing.x = unit(0.1, "lines"),
    axis.text.y = element_text(margin=margin(0,0,0,0)),
    plot.margin = margin(0,0,0,0)
  )
  
combined_plot_legend_crop <- (mutationCount_plot_crop / annot_block_crop) +
  plot_layout(heights = c(5,0.75)) + 
  theme(plot.margin = margin(l=0,r=0,b=0))
combined_plot_legend_crop

ggsave("results/Manuscript_figures/Fig_2/chip_mutation_count_cropped.png", combined_plot_legend_crop, width = 3.85, height = 3, units = "in", dpi = 300)

