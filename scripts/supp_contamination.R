# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

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

blood_tissues <- c("Whole blood", "Plasma", "Buffy coat")
hema_tissues <- c("Bone marrow", "Whole blood", "Plasma", "Buffy coat")

skyscraper_blood_contamination <- skyscraper_prep %>%
  group_by(mutPosition, HGVSp_Short) %>%
  mutate(
    #in_blood = any(Tissue %in% blood_tissues & t_alt_count > 0)
    #in_nonblood = any(!(Tissue %in% blood_tissues) & t_alt_count > 0)
    in_blood = any(Tissue %in% hema_tissues & t_alt_count > 0)
  ) %>%
  ungroup() %>%
  filter(in_blood) %>%
  dplyr::select(-in_blood)

mut_info <- skyscraper_blood_contamination %>%
  group_by(mutPosition, HGVSp_Short) %>%
  summarise(
    n_tissues  = n_distinct(Tissue[t_alt_count > 0]),
    in_buffy   = any(Tissue == "Buffy coat"  & t_alt_count > 0),
    in_plasma  = any(Tissue == "Plasma"      & t_alt_count > 0),
    in_whole   = any(Tissue == "Whole blood" & t_alt_count > 0),
    in_marrow  = any(Tissue == "Bone marrow" & t_alt_count > 0),
    .groups    = "drop"
  ) %>%
  mutate(
    blood_group = case_when(
      in_buffy              ~ "Buffy coat",
      !in_buffy & in_plasma ~ "Plasma",
      !in_buffy & !in_plasma & in_whole ~ "Whole blood",
      !in_buffy & !in_plasma & !in_whole & in_marrow ~ "Bone marrow",
      TRUE                  ~ "Other"
    ),
    blood_group = factor(
      blood_group,
      levels = c("Bone marrow", "Whole blood", "Plasma", "Buffy coat" )
    )
  ) %>%
  arrange(blood_group, n_tissues) 

mutation_order <- mut_info$HGVSp_Short

skyscraper_blood_ordered <- skyscraper_blood_contamination %>%
  mutate(HGVSp_Short = factor(HGVSp_Short, levels = mutation_order))

blood_contamination <- skyscraper_blood_ordered %>%
  ggplot(
    aes(
      x = Tissue_ordered,
      y = HGVSp_Short,
      fill = color_group,
      label = t_alt_count
    )
  ) +
  
  geom_tile(width = 1, height = 1) +
  geom_text(size = 8*25.4/72.27) +
  Col.amClass.fill +
  guides(fill = guide_legend(title = "Pathogenicity class", nrow = 2)) +
  theme_bw() + 
  scale_x_discrete(labels = custom_label, expand = c(0,0)) +
  theme(#axis.text.x = element_markdown(angle = 90, vjust = 0.5,size = 8),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, color = "black"),
    #axis.title.y = element_markdown(family = "sans", size=8, color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    #legend.position = c(0.25,0.7),
    legend.position = "top"
  )

blood_contamination_slide <- blood_contamination +
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1, size = 8, margin=margin(t=0)),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    
    legend.text = element_text(size = 8, margin=margin(l=2)),
    legend.title = element_text(size=8, margin=margin(0,4,0,0)),
    legend.key.size = unit(8, "pt"),
    legend.key.spacing.x = unit(4, "pt"),
    legend.margin = margin(0,0,0,-75)
  )


show(blood_contamination_slide)
ggsave("results/Manuscript_figures/Fig_S8/blood_contamination_tissues.png", blood_contamination_slide, width = 5, height = 6, units = "in", dpi = 300)
