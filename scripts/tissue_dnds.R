# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

#### Hunter Colegrove
#### 28 Jul 2025
#### dN/dS of tp53 in all tissues using classic method



tissue_order <- tissue_order_dnds
### input all possible mutations file from omega
### mutations in chip bed file regions

all_possible_muts_path <- "inputs/all_possible_sites_annotated.tsv.gz"
all_possible_muts <- read_delim(all_possible_muts_path, delim="\t")

## calculate expected mutations
expected_all <- all_possible_muts %>%
  filter(IMPACT != "intron_variant") %>%
  group_by(GENE) %>%
  summarise(
    synonymous = sum(IMPACT == "synonymous"),
    nonsynonymous = sum(IMPACT == "missense" | IMPACT == "nonsense" | IMPACT == "splice_region" | IMPACT == "essential_splice"),
    .groups = "drop"
  ) %>% 
  filter(synonymous >0 | nonsynonymous > 0)

############### group by LFS and distinguish syn/non-syn
observed_all_LFS <- maf_masked_coding %>%
  filter(
    Subject == "Patient", 
    Hugo_Symbol == "TP53",
    !is.na(exon_number) | Variant_Classification == "Splice_Site",
  ) %>%
  mutate(
    # Create new grouping variable
    Group = if_else(
      Subject %in% c("Patient", "Family member A", "Family member C"),
      "LFS",  # LFS group
      "non-LFS"   # Non-LFS group
    ),
    # Classify variants
    type = case_when(
      Variant_Classification %in% c("Silent") ~ "synonymous",
      Variant_Classification %in% c(
        "Missense_Mutation", "Nonsense_Mutation",
        "Nonstop_Mutation", "Splice_Site", "Splice_Region"
      ) ~ "nonsynonymous",
      TRUE ~ "exclude"
    )
  ) %>%
  group_by(Hugo_Symbol, Tissue, Group) %>%
  summarise(
    syn_obs = sum(type == "synonymous"),
    nonsyn_obs = sum(type == "nonsynonymous"),
    .groups = "drop"
  )

#### group samples
observed_all_LFS_grouped <- observed_all_LFS %>%
  left_join(tissue_categories, by="Tissue") %>%
  # replace any NA with "Other" if some samples not mapped
  mutate(Category = if_else(is.na(Category), "Other", Category)) %>%
  group_by(Hugo_Symbol, Group, Category) %>%
  summarise(
    syn_obs = sum(syn_obs, na.rm=TRUE),
    nonsyn_obs = sum(nonsyn_obs, na.rm=TRUE),
    .groups = "drop"
  )
observed_all_LFS_grouped
observed_all_LFS_summary <- observed_all_LFS_grouped %>%
  summarise(
    Hugo_Symbol = "TP53",
    Category = "All",
    Group = "LFS",
    syn_obs = sum(syn_obs, na.rm = TRUE),
    nonsyn_obs = sum(nonsyn_obs, na.rm = TRUE)
  )
observed_all_LFS_summary
observed_all_LFS <- bind_rows(observed_all_LFS_grouped, observed_all_LFS_summary)

# calculate dn/ds
dnds_all_LFS <- expected_all %>%
  inner_join(observed_all_LFS, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    pseudo_count_syn = (synonymous/(synonymous + nonsynonymous)),
    pseudo_count_non = (nonsynonymous/(synonymous + nonsynonymous)),
    dnds = ((nonsyn_obs + pseudo_count_non)/ nonsynonymous) / ((syn_obs + pseudo_count_syn)/ synonymous)
  )


## remove cases where no synonymous and no nonsynonymous
dnds_all_LFS <- dnds_all_LFS %>%
  filter(syn_obs >0 | nonsyn_obs > 0)

set.seed(1)
n_reps <- 200
df <- dnds_all_LFS %>%
  dplyr::select(GENE, Category, Group, synonymous, nonsynonymous, syn_obs, nonsyn_obs,
         pseudo_count_syn, pseudo_count_non)

boot_list <- vector("list", nrow(df))

for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  n <- r$syn_obs + r$nonsyn_obs
  
  labels <- c(rep("syn", r$syn_obs), rep("non", r$nonsyn_obs))
  
  draws <- replicate(n_reps, sample(labels, size = n, replace = TRUE))
  if (is.vector(draws)) draws <- matrix(draws, nrow = n)
  
  syn_b <- colSums(draws == "syn")
  non_b <- n - syn_b
  
  syn_adj <- syn_b + r$pseudo_count_syn
  non_syn_adj <- non_b + r$pseudo_count_non
  non_sites <- r$nonsynonymous
  syn_sites <- r$synonymous
  
  dnds <- ((non_syn_adj / non_sites) /
             (syn_adj / syn_sites))
  
  boot_list[[i]] <- tibble(
    GENE = r$GENE,
    Tissue = r$Category,
    Group = r$Group,
    draw = seq_len(n_reps),
    syn_b = syn_b,
    non_b = non_b,
    dnds = dnds
  )
}

boot_res <- bind_rows(boot_list)

boot_iqr <- boot_res %>%
  group_by(GENE, Tissue) %>%
  summarise(
    dnds_q25 = quantile(dnds, 0.25, na.rm = TRUE),
    dnds_med = median(dnds, na.rm = TRUE),
    dnds_q75 = quantile(dnds, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


dnds_all_LFS <- dnds_all_LFS %>%
  mutate(total_mut = syn_obs + nonsyn_obs) %>%
  mutate(label_fraction = paste0(nonsyn_obs, "\n―\n", syn_obs))

dnds_all_LFS <- dnds_all_LFS %>% left_join(boot_iqr, by=c("GENE" = "GENE", "Category" = "Tissue"))
group_labels <- c("All" = "All", "Blood" = "Blood", "Cancer" = "Cancer", "Solid" = "Solid\ntissues", "Sun-exposed skin" = "SkinSE")
pd <- position_dodge(width = 0.9)
breaks <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
dnds_classic <- ggplot(dnds_all_LFS, aes(x = Category, y = dnds, color = ifelse(Category == "All", "All", "Other"), fill = ifelse(Category == "All", "All", "Other"))) +
  scale_color_manual(values = c("All" = "#555555", "Other" = "#888888")) +
  scale_fill_manual(values = c("All" = "#555555", "Other" = "#DDDDDD")) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_errorbar(data = boot_iqr,
                aes(x = Tissue, ymin = dnds_q25, ymax = dnds_q75),
                position = pd, width = 0.2, inherit.aes = FALSE, na.rm = TRUE) +
  geom_point(shape = 21, size = 2.5, stroke = 1) +
  geom_text(
    aes(y = -0.1, label = label_fraction),
    position = position_dodge(width = 0.9),
    vjust = +1.15,
    hjust = 0.44,
    size = 3,
    lineheight = 0.5,
    color = 'black'
  ) +
  scale_y_continuous(limits = c(-1.3, 6.7), breaks = c(0,1,2,3,4,5,6), labels = c(0,1,2,3,4,5,6)) +
  scale_x_discrete(labels = group_labels) +
  labs(x = "Gene", y = expression(italic(d)[N]/italic(d)[S])) +
  theme_minimal() +
  theme(
    #axis.text.x = element_blank(),
    axis.text.x = element_markdown(angle = 0, hjust = 0.5, vjust = 1, size=8),
    legend.title = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=8),
    axis.title.y = element_text(size=8),
    legend.position = "none",
    plot.margin = margin(1,1,1,1)
  )
dnds_classic


#ggsave("results/dnds_naive_tissues_grouped_ms.png", dnds_classic, width = 3.75, height = 2, units = "in", dpi = 300)


## try smaller size in height
dnds_classic <- dnds_classic +
  scale_y_continuous(limits = c(-2, 6.7), breaks = c(0,1,2,3,4,5,6), labels = c(0,1,2,3,4,5,6)) +
  scale_color_manual(values = c("All" = "#888888", "Other" = "#888888")) +
  scale_fill_manual(values = c("All" = "#DDDDDD", "Other" = "#DDDDDD"))
dnds_classic


#ggsave("results/dnds_naive_tissues_grouped_ms.png", dnds_classic, width = 3.75, height = 1.5, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_4/tissue_dnds.png", dnds_classic, width = 3.75, height = 1.5, units = "in", dpi = 300)