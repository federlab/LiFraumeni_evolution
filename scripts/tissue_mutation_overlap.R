# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

base_theme <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

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
  #group_by(prot.pos, Tissue, VAF, color_group, isHotspot, protein_variant) %>%
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

## 248 only
df_agg_top <- df_agg_top %>% filter(codon_number == 248) 
df_agg_top %>% filter(color_group == "likely_benign_LC") %>% arrange(desc(max_t_alt_count))

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
    legend.position = "none",
    plot.margin = margin(t=2,b=2, r = 2, l = 2)
  ) + guides(fill = guide_legend(nrow=2))
show(heatmap_codon)

ggsave("results/Manuscript_figures/Fig_4/tissue_mutation_overlap_heatmap.png", heatmap_codon, width = 3.75, height = 1, units = "in", dpi = 300)

###########
# Stack with the tissue_skyscraper.R plot
###########

shared_x_scale <- scale_x_discrete(
  limits = tissue_order,
  labels = custom_label,
  expand = c(0,0)
)

heatmap_codon <- heatmap_codon + labs(y=" ") + theme(axis.text.x = element_blank(), 
                                       plot.margin = margin(l=2, r=2, t=0, b=-3))
skyscraper_slide <- skyscraper_slide + theme(axis.title.y = element_markdown(margin=margin(r=-10), hjust = 0.5),
                                                                         plot.margin = margin(l=2, r=2, t=2, b=2))


## with duplex reads
stacked_plot <- skyscraper_slide/heatmap_codon + 
  plot_layout(heights=c(3,1), axes="collect_x") &
  theme(plot.margin=margin(2,2,-1,2))

## without duplex reads
stacked_plot <- skyscraper_slide/heatmap_codon + 
  plot_layout(heights=c(1,1), axes="collect_x") &
  theme(plot.margin=margin(2,2,-1,2))

## without duplex reads and 248 only
stacked_plot <- skyscraper_slide/heatmap_codon + 
  #plot_layout(heights=c(2.4,0.6), axes="collect_x") &
  plot_layout(heights=c(2.5,0.5)) &
  theme(plot.margin=margin(2,2,-2,2),
        axis.title.y.left = element_blank())

stacked_plot <- ggdraw(stacked_plot) +
  draw_label(
    expression(paste("# of ", italic("TP53"), " mutations")),
    x = 0.05, y = .65, angle = 90,
    #vjust = 0, hjust = 0.5,
    size = 8
  )


#### add depths
sample_id_mapping_path <- "inputs/sampleID_mapping.txt"

sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))

tissue_depths <- final_masked_depth %>% 
  #filter(InBed == TRUE & Samp %in% samples_with_chip) %>%
  filter(gene_name == "TP53") %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue) %>%
  summarise(mean_DP = mean(DP), .groups = "drop") %>%
  filter(subject == "Patient") %>%
  dplyr::rename(Tissue = "tissue")

### add depth
axis_table <- skyscraper_prep %>%
  distinct(Tissue, Subject, Tissue_ordered) %>%
  left_join(tissue_depths, by = "Tissue")
skyscraper_slide <- skyscraper_slide + theme(axis.text.x = element_blank())

p_depth <- ggplot(axis_table, aes(Tissue, " ", fill = mean_DP)) +
  geom_tile(height = 1) +
  scale_x_discrete(limits = levels(skyscraper_prep$Tissue_ordered),
                   labels = custom_label) +
  scale_fill_gradient(high = "#9944AA", low = "#DDDDDD", breaks = c(5000,15000,25000)) +
  labs(y = "Depth", fill = "Depth") +
  base_theme +
  theme(axis.text.x = element_blank(),
    #axis.text.x = element_markdown(size =8, angle = 90),
    axis.text.y = element_text(color = "grey20", size = 8, margin=margin(0,0,0,0)),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 0)),
    legend.key.size = unit(8, "pt"),
    legend.text = element_text(size = 8, margin = margin(l=1)),
    legend.title = element_text(size = 8, margin = margin(b=2)),
    legend.margin = margin(0,0,0,-5))

tissue_depths_denominator <- final_masked_depth %>% 
  filter(gene_name == "TP53") %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>%
  filter(subject == "Patient") %>%
  dplyr::rename(Tissue = "tissue")
mf_mb_scale_prep <- skyscraper_prep %>%
  group_by(Tissue) %>%
  summarise(mutcount = n(), 
            burden = sum(t_alt_count)) %>%
  left_join(tissue_depths_denominator, by = "Tissue") %>%
  mutate(MF = mutcount/denominator,
         MB = burden/denominator)

mf_scale <- ggplot(mf_mb_scale_prep, aes(Tissue, " ", fill = MF)) +
  geom_tile(height = 1) +
  scale_x_discrete(limits = levels(skyscraper_prep$Tissue_ordered),
                   labels = custom_label) +
  scale_fill_gradient(high = "black", low = "#DDDDDD", breaks = c(1e-7, 5e-7, 1e-6, 5e-6), labels = c(1e-7, 5e-7, 1e-6, 5e-6)) +
  labs(y = "MF", fill = "MF") +
  base_theme +
  theme(axis.text.x = element_blank(),
        #axis.text.x = element_markdown(size =8, angle = 90),
        axis.text.y = element_text(color = "grey20", size = 8, margin=margin(0,0,0,0)),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 0)),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 8, margin = margin(l=1)),
        legend.title = element_text(size = 8, margin = margin(b=2)),
        legend.margin = margin(0,0,0,-5))

skyscraper_slide <- skyscraper_slide + 
  theme(axis.title.y = element_blank(),
        #legend.position = c(.4,.75), # 4 in width
        legend.position = c(.4,.75), # 3.75 in width
        axis.text.x = element_blank())
p_depth <- p_depth + 
  #coord_cartesian(clip = "off") +
  theme(axis.title.y = element_blank(),
        #legend.position = c(1.11, 3.5)) # 4 in width
        legend.position = c(1.1, 10)) # 3.75 in width
mf_scale <- mf_scale + 
  coord_cartesian(clip = "off") +
  theme(axis.title.y = element_blank(),
        #legend.position = c(1.3,4.9)) # 4 in width
        legend.position = c(1.095,5.35)) # 3.75 in width
#mb_scale <- mb_scale + theme(axis.title.y = element_blank())
variant_type_proportion <- variant_type_proportion + 
  theme(axis.title.y = element_blank(),
        #legend.position = c(1.135,0.6)) # 4 in width
        legend.position = c(1.122,0.6)) # 3.75 in width

stacked_plot <- (skyscraper_slide / p_depth / mf_scale / variant_type_proportion ) +
  plot_layout(heights = c(3, 0.2, 0.2, 0.5), # 4 in width
              guides = 'keep') &
  #theme(plot.margin = margin(2,38,2,6)) # 4 in width
  theme(plot.margin = margin(2,25,2,7)) # 3.75 in width


stacked_plot
ggsave("results/Manuscript_figures/Fig_4/tissue_mutation_overlap.png", stacked_plot, width = 3.75, height = 3.25, units = "in", dpi = 300)


