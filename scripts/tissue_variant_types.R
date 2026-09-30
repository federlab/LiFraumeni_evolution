# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Li-Fraumeni skyscraper plot


################################################################################
######### Mutation Types
################################################################################



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

mutation_type_prep <- skyscraper_prep %>%
  mutate(Mutation_type = factor(Mutation_type, levels = c(
    "Silent",
    "Splice",
    "Indel",
    "Nonsense_Mutation",
    "Missense_Mutation"
  )))


variant_colors <- c(
  "Indel"             = "#CC6677",
  "Missense_Mutation" = "#DDCC77", 
  "Nonsense_Mutation" = "#117733", 
  "Silent"            = "#88CCEE", 
  "Splice"            = "#999933"  
)

variant_type_proportion <- ggplot(mutation_type_prep, aes(x = Tissue_ordered, fill = Mutation_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = variant_colors,     
    labels = c(
    "Indel"             = "Indel",
    "Missense_Mutation" = "Missense",
    "Nonsense_Mutation" = "Nonsense",
    "Silent"            = "Silent",
    "Splice"            = "Splice"
  )
  ) +
  labs(
    x = "Tissue",
    #y = "Proportion",
    y = "Prop.",
    fill = "Variant Classification"
  ) +
  theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  scale_y_continuous(breaks = c(0, 1), labels = c("0","1")) +
  theme(
        axis.text.x.bottom = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size=8),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=8, margin=margin(r=2)),
        legend.key.size = unit(8,"pt"),
        legend.key.spacing.x = unit(3,"pt"),
        legend.margin = margin(-15,0,0,0),
        axis.text.y.left = element_text(size=8, margin = margin(r=2)),
        plot.margin=margin(2,4,2,4))

variant_type_proportion
ggsave("results/Manuscript_figures/Fig_4/tissue_variant_types.png", variant_type_proportion, width = 4, height = .5, units= "in", dpi = 300)

