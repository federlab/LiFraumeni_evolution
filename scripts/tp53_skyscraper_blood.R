# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Figure 3A
## generate x-axis legend for age, depth, lfs and ctx status

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


blood_types <- c("PBMC","Buffy coat")
subject_barcode <- maf_masked_coding %>%
  filter(Subject %in% ann_wide$Subject) %>%
  filter(Tissue %in% blood_types) %>%
  group_by(Subject) %>%
  summarise(Tumor_Sample_Barcode = dplyr::first(Tumor_Sample_Barcode), .groups = "drop")

ann_wide <- ann_wide %>%
  left_join(subject_barcode, by = "Subject")

sample_barcodes <- ann_wide$Tumor_Sample_Barcode

## coding depth summaries by subject (patient buffy coat only)
coding_depth_summary <- final_masked_depth %>% 
  filter(!is.na(exon_number)) %>% 
  filter(Samp %in% sample_barcodes) %>%
  group_by(Samp) %>%
  summarise(meanDP = mean(DP))

subject_depth <- coding_depth_summary %>%
  left_join(subject_barcode, by = c("Samp" = "Tumor_Sample_Barcode")) %>%
  filter(Subject %in% ann_wide$Subject) %>%
  group_by(Subject) %>%
  summarise(Depth = mean(meanDP, na.rm = TRUE), .groups = "drop")

ann_wide <- ann_wide %>%
  dplyr::select(-Depth) %>%                # drop the old High/Low labels
  left_join(subject_depth, by = "Subject")

ann_wide

skyscraper_prep <- maf_masked_coding %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter((!is.na(exon_number)) | (Variant_Classification %in% c("Splice_Site"))) %>%
  filter(!inRepeatMask | (Variant_Classification %in% c("Splice_Site"))) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  mutate(
    color_group = case_when(
      am_class == "likely_benign"     & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign"     & t_alt_count >  1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count >  1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous"         & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous"         & t_alt_count >  1 ~ "ambiguous_LC"
    )
  ) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  arrange(desc(t_alt_count), path_am_order) %>% 
  mutate(SampCodingOrder = rank(dplyr::desc(t_alt_count), ties.method = "first"))

axis_table <- ann_wide %>%
  mutate(LFS_CTx = case_when(
    LFS == "Yes" & CTx == "Yes" ~ "LFS/CTx",
    LFS == "Yes" & CTx == "No"  ~ "LFS/no-CTx",
    LFS == "No"  & CTx == "No"  ~ "non-LFS/no-CTx",
    LFS == "No"  & CTx == "Yes" ~ "non-LFS/CTx"
  ))

skyscraper_prep <- skyscraper_prep %>% 
  left_join(axis_table, by=c("Subject", "Tumor_Sample_Barcode")) %>%
  ungroup() %>%
  mutate(SampleCode = fct_reorder(SampleCode, Age, .fun = min))

skyscraper <- skyscraper_prep %>%
  ggplot(aes(x = SampleCode, y = SampCodingOrder, fill = color_group, label = t_alt_count)) +
  geom_tile() +
  geom_text(size = 8*25.4/72.27) +
  geom_tile(aes(x = SampleCode, y = SampCodingOrder), inherit.aes = FALSE,
            data = skyscraper_prep, color = "black", fill = NA, linewidth = .1) +
  scale_y_continuous(
    expand = c(0,0),
    breaks = c(0.5,5.5,10.5,15.5,20.5),
    labels = c(0,5,10,15,20),
    limits = c(0,23),
    name = "Number of <i>TP53</i> mutations"
  ) +
  Col.amClass.fill +
  theme_bw() +
  
  theme(
    axis.text.x  = element_blank(), 
    #axis.text.x  = element_text(angle = 90),
    axis.ticks.x = element_blank(),  
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "#DDDDDD", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=8, color = "black", margin=margin(r=0, unit = "pt")),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7)
  )
skyscraper


axis_table <- axis_table %>% mutate(SampleCode = fct_reorder(SampleCode, Age, .fun = min))

base_theme <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

p_age <- ggplot(axis_table, aes(SampleCode, " ", fill = Age)) +
  geom_tile(height = 1) +
  scale_fill_gradient(low = "#DDDDDD", high = "#555555") +
  labs(y = "Age", fill = "Age") +
  base_theme + theme(axis.text.x = element_blank(),
                     axis.text.y = element_text(color = "grey20", size = 8, margin=margin(0,0,0,0)),
                     axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 0)),
                     legend.key.size = unit(8, "pt"))


p_lfs <- ggplot(axis_table, aes(SampleCode, " ")) +
  geom_tile(aes(fill = LFS_CTx, colour = LFS_CTx),
            height = 1, width = 0.8, linewidth = 0.8) +
  scale_fill_manual(values = c(
    "LFS/CTx"          = "#882255",     # filled maroon
    "LFS/no-CTx"       = "#DDDDDD", # outlined maroon
    "non-LFS/CTx"      = "#44aa99",     # filled teal
    "non-LFS/no-CTx"   = "#DDDDDD"  # outlined teal
  ), na.value = "transparent", drop = FALSE) +
  scale_colour_manual(values = c(
    "LFS/CTx"          = "#882255",
    "LFS/no-CTx"       = "#882255",
    "non-LFS/CTx"      = "#44aa99",
    "non-LFS/no-CTx"   = "#44aa99"
  ), drop = FALSE) +
  labs(y = "LFS/CTx", fill = "LFS/CTx", colour = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8),
    axis.text.y  = element_text(color = "grey20", size=8, margin = margin(0,0,0,0)),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 0)),
    legend.key.size = unit(8, "pt")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        linewidth = 0.8,
        colour = c("#882255", "#882255", "#44aa99", "#44aa99")
      )
    ),
    colour = "none"
  )
p_lfs

p_depth <- ggplot(axis_table, aes(SampleCode, " ", fill = Depth)) +
  geom_tile(height = 1) +
  scale_fill_gradient(high = "#9944AA", low = "#DDDDDD") +
  labs(y = "Depth", fill = "Depth") +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 8, margin=margin(0,0,0,0)),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 0)),
        legend.key.size = unit(8, "pt"))

axis_table
skyscraper_inset <- skyscraper +
  theme(
    #legend.position = c(1, 0.95),
    legend.justification = c(0, 1),
    #legend.background = element_rect(fill = "white", colour = "#DDDDDD"),
    legend.background = element_blank(),
    legend.key.size = unit(8, "pt"),
    legend.title = element_text(size = 8, margin=margin(b=2)),
    legend.text  = element_text(size = 8),
    plot.margin = margin(t=0, r=0, b=0, l=0)
  )

skyscraper_inset <- skyscraper_inset +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.position = c(0,1.03),
        axis.title.y = element_blank()) 

p_age   <- p_age   + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                           legend.position = c(1.08,7))
  
p_depth <- p_depth + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                           legend.position = c(1.29,8.3))

p_lfs   <- p_lfs   + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                           legend.position = c(1.23,2.75))


annot_block <- (p_age / p_depth / p_lfs) +
  plot_layout(heights = c(1, 1, 1, 1)) &
  theme(
    #legend.position = "right",
    #legend.box = "horizontal",
    #legend.box.margin = margin(0,0,0,-10),
    #legend.justification = c(0,-0.5),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    #legend.spacing.x = unit(-0.2, "lines"),
    axis.text.y = element_text(margin=margin(0,-10,0,-10)),
    #axis.title.y = element_text(margin=margin(0,0,0,0))
    axis.title.y = element_blank()
  )

final_plot <- (skyscraper_inset / annot_block) +
  plot_layout(heights = c(2.5, 1), guides = 'keep') &
  theme(plot.margin = margin(2,38,-5,9))
final_plot

#ggsave("results/tp53_skyscraper_blood.png", final_plot, width = 3.75, height = 3.25, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_skyscraper_blood.png", final_plot, width = 3.75, height = 3.25, units = "in", dpi = 300)
