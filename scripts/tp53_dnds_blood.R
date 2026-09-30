# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

#### Hunter Colegrove
#### 28 Jul 2025
#### dN/dS of chip genes using classic method


tissues = c("PBMC", "Buffy coat")
subjects <- c("Family member A", "Family member B", "Family member C", 
                "UW volunteer 1",  "UW volunteer 2",  "UW volunteer 3",  
                "UW volunteer 4" , "UW volunteer 5" , "UW volunteer 6", 
                "UW volunteer 7",  "Patient")

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

list_sub <- unique(maf_masked_coding$Subject)

############### group by LFS and distinguish syn/non-syn
unique(maf_masked_coding$Subject)
observed_all_LFS <- maf_masked_coding %>%
  filter(!is.na(gene_name) | Variant_Classification == "Splice_Site") %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(
    Tumor_Sample_Barcode != "DevDNA1_S1.1",
    Tissue %in% tissues,
    Subject %in% subjects
  ) %>%
  mutate(
    # Create new grouping variable
    Group = case_when(
      Subject %in% c("Patient") ~ "LFS\nCTx",
      Subject %in% c("Family member A", "Family member C") ~ "LFS\nno-CTx",
      Subject %in% c("UW volunteer 7") ~ "non-LFS\nCTx", 
      TRUE ~ "non-LFS\nno-CTx"
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
  group_by(Hugo_Symbol, Group) %>%
  summarise(
    syn_obs = sum(type == "synonymous"),
    nonsyn_obs = sum(type == "nonsynonymous"),
    .groups = "drop"
  )

# calculate dn/ds
dnds_all_LFS <- expected_all %>%
  inner_join(observed_all_LFS, by = c("GENE" = "Hugo_Symbol")) %>%
  mutate(
    pseudo_count_syn = (synonymous/(synonymous + nonsynonymous)),
    pseudo_count_non = (nonsynonymous/(synonymous + nonsynonymous)),
    dnds = ((nonsyn_obs + pseudo_count_non)/ nonsynonymous) / ((syn_obs + pseudo_count_syn)/ synonymous)
  )


set.seed(456)
n_reps <- 5000
df <- dnds_all_LFS %>%
  dplyr::select(GENE, Group, synonymous, nonsynonymous, syn_obs, nonsyn_obs,
         pseudo_count_syn, pseudo_count_non)
df <- df %>% filter(GENE == "TP53")

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
    Group = r$Group,
    draw = seq_len(n_reps),
    syn_b = syn_b,
    non_b = non_b,
    dnds = dnds
  )
}

boot_res <- bind_rows(boot_list)


boot_iqr <- boot_res %>%
  group_by(GENE, Group) %>%
  summarise(
    dnds_q25 = quantile(dnds, 0.25, na.rm = TRUE),
    dnds_q75 = quantile(dnds, 0.75, na.rm = TRUE),
    #dnds_q25 = quantile(dnds, 0.025, na.rm = TRUE),
    #dnds_q75 = quantile(dnds, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


dnds_all_LFS <- dnds_all_LFS %>%
  mutate(total_mut = syn_obs + nonsyn_obs) %>%
  mutate(label_fraction = paste0(nonsyn_obs, "\n―\n", syn_obs))

dnds_all_LFS <- dnds_all_LFS %>% left_join(boot_iqr, by=c("GENE", "Group"))

pd <- position_dodge(width = 0.9)
breaks <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
breaks <- c(0.25, 1, 4, 16)


lower_lim = -3.5
upper_lim = 30
breaks = c(0,10,20,30)

dnds_all_LFS <- dnds_all_LFS %>% filter(GENE == "TP53")
dnds_all_LFS
group_order <- c(
  "LFS\nno-CTx",
  "non-LFS\nno-CTx",
  "LFS\nCTx",
  "non-LFS\nCTx"
)

dnds_all_LFS <- dnds_all_LFS %>%
  mutate(Group = factor(Group, levels = group_order))

boot_iqr
boot_iqr <- boot_iqr %>%
  mutate(Group = factor(Group, levels = group_order))


dnds_all_LFS <- dnds_all_LFS %>%
  mutate(
    fill_group = case_when(
      Group == "LFS\nCTx"   ~ "fill_x",
      Group == "LFS\nno-CTx"  ~ "outline_x",
      Group == "non-LFS\nCTx"   ~ "fill_y",
      Group == "non-LFS\nno-CTx"  ~ "outline_y",
      TRUE                       ~ "other"
    )
  )

boot_iqr <- boot_iqr %>%
  mutate(
    fill_group = case_when(
      Group == "LFS\nCTx"   ~ "fill_x",
      Group == "LFS\nno-CTx"  ~ "outline_x",
      Group == "non-LFS\nCTx"   ~ "fill_y",
      Group == "non-LFS\nno-CTx"  ~ "outline_y",
      TRUE                       ~ "other"
    )
  )

fill_scale <- scale_fill_manual(
  values = c(
    fill_x    = "#882255",
    outline_x = NA,
    fill_y    = "#44aa99",
    outline_y = NA
  )
)

color_scale <- scale_color_manual(
  values = c(
    fill_x    = "#882255",
    outline_x = "#882255",
    fill_y    = "#44aa99",
    outline_y = "#44aa99"
  )
)

dnds_all_LFS

dnds_classic <- ggplot(dnds_all_LFS, aes(x = Group, y = dnds)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_errorbar(data = boot_iqr,
                aes(x = Group, ymin = dnds_q25, ymax = dnds_q75, color = fill_group),
                position = pd, width = 0.2, inherit.aes = FALSE, na.rm = TRUE) +
  geom_point(aes(color = fill_group, fill = fill_group), shape = 21, size = 2.5, stroke = 1) +
  geom_text(
    aes(y = -1, label = label_fraction),
    position = position_dodge(width = 1),
    vjust = 1,
    hjust = 0.5,
    size = 8*25.4/72.27,
    lineheight = 0.5,
    color = 'black'
  ) +
  fill_scale +
  color_scale +
  scale_y_continuous(limits = c(lower_lim, upper_lim), breaks = breaks) +
  labs(x = "Gene", y = expression(italic(d)[N]/italic(d)[S])) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size=8, margin=margin(0,0,0,0)),
    #axis.text.x = element_blank(),
    axis.text.y = element_text(size=8, margin=margin(0,-2,0,-2)),
    axis.title.y = element_text(size=8),
    legend.title = element_blank(), 
    axis.title.x = element_blank(),
    #legend.position = c(0.45,0.8),
    legend.position = "none",
    legend.key.size = unit(8, "pt"),
    plot.margin = margin(l=2,r=0,b=0,t=0),
    legend.text = element_text(size=8)
  )

dnds_classic


#ggsave("results/tp53_dnds_blood.png", dnds_classic, width = 2, height = 2.2, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_dnds_blood.png", dnds_classic, width = 2, height = 2.2, units = "in", dpi = 300)


boot_iqr %>% filter( Group %in% c("LFS\nno-CTx", "non-LFS\nno-CTx"))
dnds_classic_no_chemo <- ggplot(dnds_all_LFS %>% filter(Group %in% c("LFS\nno-CTx", "non-LFS\nno-CTx")), aes(x = Group, y = dnds)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_errorbar(data = boot_iqr %>% filter(Group %in% c("LFS\nno-CTx", "non-LFS\nno-CTx")),
                aes(x = Group, ymin = dnds_q25, ymax = dnds_q75, color = fill_group),
                position = pd, width = 0.2, inherit.aes = FALSE, na.rm = TRUE) +
  geom_point(aes(color = fill_group, fill = fill_group), shape = 21, size = 2.5, stroke = 1) +
  fill_scale +
  color_scale +
  # geom_text(
  #   aes(y = -.4, label = label_fraction),
  #   position = position_dodge(width = 1),
  #   vjust = 1,
  #   hjust = 0.5,
  #   size = 8*25.4/72.27, 
  #   lineheight = 0.5
  # ) +
  scale_y_continuous(limits = c(0, 2), breaks = c(0,1, 2)) +
  labs(x = "Gene", y = expression(italic(d)[N]/italic(d)[S])) +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(hjust = 0.5, size=8, margin=margin(0,0,0,0)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8, margin=margin(0,-2,0,-2)),
    #axis.title.y = element_text(size=8),
    axis.title.y = element_blank(),
    legend.title = element_blank(), 
    axis.title.x = element_blank(),
    #legend.position = c(0.45,0.8),
    legend.position = "none",
    legend.key.size = unit(8, "pt"),
    plot.margin = margin(l=1,r=0,b=2,t=3),
    legend.text = element_text(size=8)
  )

show(dnds_classic_no_chemo)
#ggsave("results/tp53_dnds_blood_nochemo_inset.png", dnds_classic_no_chemo, width = 0.75, height = 0.75, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_dnds_blood_nochemo_inset.png", dnds_classic_no_chemo, width = 0.75, height = 0.75, units = "in", dpi = 300)
