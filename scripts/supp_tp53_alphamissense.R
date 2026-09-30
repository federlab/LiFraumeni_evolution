# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington




## all possible mutations
all_possible_muts_path <- "inputs/all_possible_sites_annotated.tsv.gz"
all_possible_muts <- read_delim(all_possible_muts_path, delim="\t")

## assign alphamissense to all_possible missense mutations and exclude the rest
all_missense <- all_possible_muts %>%
  filter(GENE == "TP53") %>%
  left_join(alphamissense, by = c("POS" = "POS", "REF" = "REF", "ALT" = "ALT")) %>%
  filter(!is.na(am_pathogenicity))

#######################
######### All genes alphamissense supp figure
########################

all_genes <- read_delim("results/all_genes_am.tsv")


blood_filtering <- all_genes %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(!is.na(am_pathogenicity.x) | !is.na(am_pathogenicity.y)) %>%
  mutate(am_pathogenicity.y = if_else(is.na(am_pathogenicity.y), am_pathogenicity.x, am_pathogenicity.y))

blood_filtering$group <- "Observed"
all_missense$group <- "All possible"
combine_df <- bind_rows(blood_filtering, all_missense)

lfs_subjects
ctx_subjects
am_groups <- combine_df %>% 
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>% 
  mutate(LFS = if_else(is.na(Subject), "All\npossible", LFS)) %>%
  mutate(CTX = if_else(Subject %in% ctx_subjects, "CTX", "no-CTX")) %>% 
  mutate(CTX = if_else(is.na(Subject), "", CTX)) %>%
  mutate(lfs_group = paste(LFS,CTX, sep='\n')) 



am_groups <- am_groups %>% filter(gene_name %in% c("ASXL1", "BRINP3", "DNMT3A", "GATA2", "RAD21", "TET2", "TP53"))

mw_results <- am_groups %>%
  filter(group != "All possible") %>%        # same filtering as your plot
  filter(gene_name %in% c("ASXL1", "BRINP3", "DNMT3A", "GATA2", "RAD21", "TET2", "TP53")) %>%
  group_by(gene_name) %>%
  summarise(
    test = list(
      wilcox.test(
        am_pathogenicity.y ~ LFS, exact = FALSE
      )
    )
  ) %>%
  mutate(tidy = map(test, tidy)) %>%
  unnest(tidy) %>%
  dplyr::select(gene_name, statistic, p.value, method)
mw_results

obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = LFS, y = am_pathogenicity.y, color = LFS)) +
  facet_wrap(~ gene_name, nrow = 1) +

  geom_boxplot(outlier.shape = NA) +

  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +

  labs(y = "AlphaMissense\npathogenicity", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups

#ggsave("results/AlphaMissense_chip_supp.png", obsv_am_blood_groups, width = 6, height = 2, units = "in", dpi = 300)


pathogenic_proportion <- am_groups %>% mutate(pathogenic = if_else(am_pathogenicity.y > 0.564, "Pathogenic", "Non-Pathogenic")) %>%
  group_by(gene_name, LFS) %>%
  summarise(n_total = n(),
            n_path = sum(pathogenic == "Pathogenic"),
            n_nonpath = n_total-n_path,
            prop_path = n_path/n_total,
            .groups = "drop")
pathogenic_proportion

fishers_path <- pathogenic_proportion %>% 
  dplyr::select(gene_name, LFS, n_total, n_path, n_nonpath) %>%
  pivot_wider(
    names_from = LFS,
    values_from = c(n_total, n_path, n_nonpath),
    names_sep = "_"
  )
pathogenic_proportion
df_n <- pathogenic_proportion %>%
  distinct(gene_name, LFS, n_total) %>%
  mutate(n_lab = paste0("(n = ", n_total, ")"))

prop_plot <- ggplot(pathogenic_proportion, aes(x = LFS, y = prop_path, fill = LFS)) +
  facet_wrap(~ gene_name, nrow = 1) +
  geom_col(width = 0.7, color = "black") +
  geom_text(
    aes(label = sprintf("%.2f", prop_path)),
    vjust = -0.3,
    size = 2.6,
    color = "black"
  ) +
  geom_text(
    data = df_n,
    aes(x = LFS, y = -0.3, label = n_lab),
    inherit.aes = FALSE,
    vjust = 1.6,
    size = 2.6
  ) +
  scale_fill_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +
  scale_y_continuous(limits = c(-0.3, 1), expand = expansion(mult = c(-0.2, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(y = "Pathogenic\nproportion", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 10, 2, "pt")  # extra bottom room for the (n=...)
  )


#ggsave("results/AlphaMissense_proportion_supp.png", prop_plot, width = 6, height = 2, units = "in", dpi = 300)


combined <- obsv_am_blood_groups / prop_plot +
  plot_layout(heights = c(2, 1))

combined
ggsave("results/Manuscript_figures/Fig_S6/supp_tp53_alphamissense.png", combined, width = 7, height = 3, units = "in", dpi = 300)



####################################################
### Figure 3 LFS/CTX grouped alphamissense for tp53
####################################################

am_groups <- am_groups %>% filter(gene_name %in% c("TP53"))

unique(am_groups$lfs_group)
wilcox.test(am_pathogenicity.y ~ lfs_group, 
            data = am_groups %>% filter(lfs_group %in% c("LFS\nno-CTX", "non-LFS\nCTX")), exact = FALSE)

group_order <- c(
  "LFS\nno-CTX",
  "LFS\nCTX",
  "non-LFS\nno-CTX",
  "non-LFS\nCTX"
)
am_groups <- am_groups %>%
  mutate(lfs_group = factor(lfs_group, levels = group_order))


obsv_am_blood_groups <- ggplot(am_groups %>% filter(group != "All possible"), aes(x = lfs_group, y = am_pathogenicity.y)) +

  geom_boxplot(outlier.shape = NA, color = "black") +

  geom_jitter(width = 0.15, size = 1.5, stroke = 1.2, alpha = 1, aes(color = LFS, shape = CTX)) +
  scale_color_manual(values = c("LFS" = "#882255", "non-LFS" = "#44aa99")) +
  scale_shape_manual(values = c("CTX" = 16, "no-CTX" = 1)) +
  labs(y = "Pathogenicity score", x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "pt")
  )
obsv_am_blood_groups

#ggsave("results/tp53_alphamissense_blood.png", obsv_am_blood_groups, width = 3.5, height = 1.5, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_alphamissense_blood.png", obsv_am_blood_groups, width = 3.5, height = 1.5, units = "in", dpi = 300)

