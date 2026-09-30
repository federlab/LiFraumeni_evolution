# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Li-Fraumeni skyscraper plot

abbreviations <- c("WB", "Buffy", "Plas", "BM", "Bucc", "Thyr", "Bron", "Lung",
                   "Eso1", "Eso2", "Gast1", "Gast2", "CardM", "Spln", "Liver", "Colon",
                   "Omen", "Perit", "Renal", "Testis", "SkelM", "Skin", "SkinNS",
                   "MedMet", "LungMet", "EsoCa1", "EsoCa2", "LivMet1", "LivMet2")

tissue_labels <- abbreviations
names(tissue_labels) <- tissue_order


skyscraper_prep <-
  maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>% 
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  mutate(
    color_group = case_when(
      am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC"
    )
  )   %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  #group_by(Tumor_Sample_Barcode, coding, Hugo_Symbol) %>% 
  arrange(desc(t_alt_count), path_am_order) %>% 
  mutate(SampCodingOrder = rank(dplyr::desc(t_alt_count), ties.method = "first"))

skyscraper_prep <- skyscraper_prep %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))


custom_label <- function(x) {
  sapply(x, function(t) {
    abbrev <- tissue_labels[[t]]  # look up abbreviation
    if (t %in% cancer_samples) {
      paste0("<span style='color:red;'>", abbrev, "</span>")
    } else {
      abbrev
    }
  })
}

skyscraper_prep <- skyscraper_prep %>% 
  mutate(group_type = case_when(Tissue %in% cancer_samples ~ 'cancer',
                           Tissue %in% blood_samples ~ 'blood',
                           TRUE ~ 'normal'), 
         group_type = factor(group_type, levels = c("blood", "normal", "cancer")))

tissue_max <- skyscraper_prep %>%
  group_by(group_type, Tissue) %>% 
  summarise(max_sampcoding = max(SampCodingOrder), .groups = "drop") %>%
  group_by(group_type) %>%
  arrange(max_sampcoding, .by_group = TRUE) %>%   # within cancer groups
  mutate(
    Tissue_ordered = factor(Tissue, levels = unique(Tissue))
  ) %>%
  ungroup()
tissue_max
skyscraper_prep <- skyscraper_prep %>%
  left_join(tissue_max %>% dplyr::select(group_type, Tissue, Tissue_ordered), by = c("group_type", "Tissue"))

#save tissue order for other plots:
tissue_levels_all_plots <- levels(skyscraper_prep$Tissue_ordered)

skyscraper <- skyscraper_prep %>%
  ggplot(
    aes(
      x = Tissue_ordered,
      y = SampCodingOrder,
      fill = color_group,
      label = t_alt_count
    )
  ) +
  geom_tile(width = 1, height = 1) +
  geom_tile(data = skyscraper_prep,
            aes(
              x = Tissue_ordered,
              y = SampCodingOrder, 
              height = 1, 
              width = 1
            ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
  scale_y_continuous(expand = c(0,0),
                     #breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
                     #labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
                     limits = c(0,53),
                     name = '\\# of <i>TP53</i> mutations') +
  Col.amClass.fill +
  theme_bw() + 
  scale_x_discrete(labels = custom_label, expand = c(0,0)) +
  theme(#axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 8),
        panel.spacing = unit(0,"pt"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x.top = element_markdown(family = "sans", size=8, color = "black"),
        axis.title.y = element_markdown(family = "sans", size=8, color = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.25,0.7),
  )


skyscraper_slide <- skyscraper +
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 8, margin=margin(t=0)),
    axis.text.y = element_text(size = 8, margin=margin(r=1)),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, color = "black"),
    axis.title.y = element_markdown(family = "sans", size=8, color = "black"),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    ## legend with duplex reads
    #legend.position = c(0.3,0.7),
    ## legend without duplex reads 
    legend.position = c(0.4,0.68),
    
    legend.text = (element_text(size = 8)),
    legend.title = (element_text(size=8, margin=margin(0,0,2,0))),
    legend.key.size = unit(8, "pt"),
    legend.margin = margin(0,0,0,0)
  ) 
#+ geom_text(size=8*25.4/72.27) ## with or without duplex read counts

show(skyscraper_slide)


ggsave("results/Manuscript_figures/Fig_4/tissue_skyscraper.png", skyscraper_slide, width = 4, height = 2.5, units = "in", dpi = 300)



