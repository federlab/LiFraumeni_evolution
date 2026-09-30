# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

################################################################################
############# Mutation tissue overlap - Amino acid changes
################################################################################

df_agg <- maf_masked_coding %>%
  filter(Tissue %in% tissue_order) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(Subject == "Patient") %>% 
  filter(!inRepeatMask) %>%
  filter(!is.na(prot.pos)) %>%
  mutate(
    color_group = case_when(
      am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC")) %>%
  group_by(prot.pos, Tissue, VAF, color_group, protein_variant) %>%
  summarise(max_t_alt_count = max(t_alt_count, na.rm = TRUE)) %>%
  filter(!is.na(protein_variant)) %>%
  ungroup()

lookup_prot_pos <- df_agg %>%
  filter(!is.na(prot.pos)) %>%
  dplyr::select(protein_variant, prot.pos) %>%
  distinct()


df_complete <- df_agg %>%
  #complete(prot.pos, Tissue, fill = list(max_t_alt_count = 0, VAF=0))  # Fill missing combinations with 0
  complete(protein_variant,Tissue = factor(tissue_order, levels = tissue_order),fill = list(max_t_alt_count = 0, VAF = 0)) %>%
  #complete(protein_variant, Tissue, fill = list(max_t_alt_count = 0, VAF=0)) %>% # Fill missing combinations with 0
  left_join(lookup_prot_pos, by = "protein_variant") %>%
  mutate(prot.pos = prot.pos.y) %>%
  dplyr::select(-prot.pos.x, -prot.pos.y)
df_agg <- df_complete %>% arrange(desc(VAF))

mutation_count <- df_agg %>%
  group_by(protein_variant) %>%
  summarise(mutation_count = sum(max_t_alt_count > 0), .groups = "drop")

mutation_tissue_count <- df_agg %>%
  filter(max_t_alt_count > 0) %>%
  group_by(protein_variant) %>%
  summarise(mutation_tissue_count = n_distinct(Tissue),.groups = "drop")

codon_tissue_count <- df_agg %>%
  filter(max_t_alt_count > 0) %>%
  group_by(prot.pos) %>%
  summarise(codon_tissue_count = n_distinct(Tissue),.groups = "drop")

# Merge the tissue count with df_agg
df_agg <- df_agg %>%
  left_join(mutation_count, by = "protein_variant") %>%
  left_join(mutation_tissue_count, by = "protein_variant") %>%
  left_join(codon_tissue_count, by = "prot.pos")

# Order HGVSc based on descending number of tissue matches
df_agg <- df_agg %>%
  mutate(protein_variant = fct_reorder(factor(protein_variant), mutation_tissue_count, .desc = TRUE))
df_agg

## compile list of tissues by protein variant
df_unique <- df_agg %>%
  filter(max_t_alt_count > 0) %>%
  group_by(protein_variant) %>%
  summarise(
    tissue_count = n_distinct(Tissue),  # Count number of unique tissues
    tissues = list(unique(Tissue)),  # Store unique tissues as a list
    .groups = "drop"
  ) %>%
  arrange(desc(tissue_count))  # Sort by tissue count descending


# Convert HGVSc and Tissue to factors to ensure proper ordering in the heatmap
df_agg$protein_variant <- as.factor(df_agg$protein_variant)
df_agg$Tissue <- as.factor(df_agg$Tissue)

df_agg <- df_agg %>% 
  #arrange(desc(prot.pos)) %>%
  arrange(desc(protein_variant)) %>%
  #mutate(prot.pos = factor(prot.pos, levels = unique(prot.pos)))
  mutate(protein_variant = factor(protein_variant, levels = unique(protein_variant)))


df_agg <- df_agg %>%
  mutate(Tissue = factor(Tissue, levels = levels(skyscraper_prep$Tissue_ordered)))

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

## all blood
df_agg_all_blood <- df_agg %>%
  group_by(protein_variant) %>%
  filter(any(max_t_alt_count > 0 & Tissue %in% c("Plasma", "Whole blood", "Buffy coat"))) %>%
  filter(Tissue %in% c("Plasma", "Whole blood", "Buffy coat")) %>%
  filter(max_t_alt_count > 0)
## buffy coat only
df_agg_buffycoat <- df_agg %>%
  group_by(protein_variant) %>%
  filter(any(max_t_alt_count > 0 & Tissue == "Buffy coat"))

## filter by commonly shared codons
df_agg_top <- df_agg %>% filter(codon_tissue_count > 2)


df_agg_top <- df_agg_top %>%
  mutate(codon_number = as.numeric(stringr::str_extract(protein_variant, "\\d+")))
df_agg_top <- df_agg_top %>%
  arrange(
    codon_tissue_count,
    codon_number,
    mutation_tissue_count
  )

df_agg_top <- df_agg_top %>%
  mutate(protein_variant = factor(protein_variant, levels = unique(protein_variant)))

### make boxes around codons
df_agg_top <- df_agg_top %>%
  mutate(y_pos = as.numeric(protein_variant))

codon_groups <- df_agg_top %>%
  group_by(codon_number) %>%
  summarise(
    ymin = min(y_pos) - 0.5,  # extend half-unit below bottom tile
    ymax = max(y_pos) + 0.5,  # extend half-unit above top tile
    xmin = 0.5,               # across entire x-axis
    xmax = length(unique(df_agg_top$Tissue)) + 0.5
  )


heatmap_codon <- ggplot(df_agg_top, aes(x = Tissue, y = protein_variant, fill = color_group, label = max_t_alt_count)) +
  geom_tile(color="black", size = 0.1) +
  geom_rect(
    data = codon_groups,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE
  ) +
  Col.amClass.fill +
  labs(x = "Tissue",
       y = "Codon",
       fill = "t_alt_count") +
  theme_minimal() +
  geom_text( size=8*25.4/72.27) +
  scale_x_discrete(labels = custom_label) +
  theme(
    axis.text.x.bottom = element_markdown(angle = 90,size = 8, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(8, 'pt')
  ) + guides(fill = guide_legend(nrow=2, title = "Pathogenicity class" ))
show(heatmap_codon)

ggsave("results/Manuscript_figures/Fig_S9/mutation_heatmap_most_shared_codons.png", heatmap_codon, width = 5, height = 3.5, units = "in")

