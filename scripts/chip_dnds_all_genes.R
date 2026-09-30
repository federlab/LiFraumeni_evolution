# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

#### Hunter Colegrove
#### 28 Jul 2025
#### dN/dS of chip genes using a naive method

### input all possible mutations file chip_dnds_target_annotate.R

all_possible_muts_path <- file.path("inputs/all_muts_annotated.tsv.gz")

all_possible_muts <- read_delim(all_possible_muts_path, delim="\t")
all_possible_muts
## calculate expected mutations
expected_all <- all_possible_muts %>%
  group_by(SYMBOL) %>%
  summarise(
    synonymous = sum(impact == "synonymous"),
    nonsynonymous = sum(impact == "missense" | impact == "nonsense" | impact == "splice"),
    .groups = "drop"
  ) %>%
  filter(synonymous >0 | nonsynonymous > 0)

############### group by LFS and distinguish syn/non-syn
observed_all_LFS <- maf_masked_coding %>%
  filter(!is.na(gene_name) | Variant_Classification == "Splice_Site") %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(
    Tumor_Sample_Barcode != "DevDNA1_S1.1",
    Hugo_Symbol %in% CHIP_genes,
    Tissue %in% c("PBMC", "Buffy coat")
  ) %>% 
  mutate(
    # Create new grouping variable
    Group = if_else(
      Subject %in% c("Patient", "Family member A", "Family member C"),
      "LFS",  
      "non-LFS"
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
  inner_join(observed_all_LFS, by = c("SYMBOL" = "Hugo_Symbol")) %>%
  mutate(
    pseudo_count_syn = (synonymous/(synonymous + nonsynonymous)),
    pseudo_count_non = (nonsynonymous/(synonymous + nonsynonymous)),
    dnds = ((nonsyn_obs + pseudo_count_non)/ nonsynonymous) / ((syn_obs + pseudo_count_syn)/ synonymous)
  )

## remove cases where no synonymous and no nonsynonymous
## synonymous or non-synonymous needs at least plot_n_min mutations
## synonymous + non-synonymous needs at least plot_n_total mutations
plot_n_min = 2
plot_n_total = 10
dnds_all_LFS <- dnds_all_LFS %>%
  group_by(SYMBOL) %>%
  filter(
    all(syn_obs >= plot_n_min & nonsyn_obs >= plot_n_min),
    all(all((syn_obs + nonsyn_obs) >= plot_n_total)),
    n_distinct(Group) == 2 
  ) %>%
  ungroup()

## bootstrap error bars
set.seed(456)
n_reps <- 5000
df <- dnds_all_LFS %>%
  dplyr::select(SYMBOL, Group, synonymous, nonsynonymous, syn_obs, nonsyn_obs,
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
    SYMBOL = r$SYMBOL,
    Group = r$Group,
    draw = seq_len(n_reps),
    syn_b = syn_b,
    non_b = non_b,
    dnds = dnds
  )
}

boot_res <- bind_rows(boot_list)

boot_iqr <- boot_res %>%
  group_by(SYMBOL, Group) %>%
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

dnds_all_LFS <- dnds_all_LFS %>% left_join(boot_iqr, by=c("SYMBOL", "Group"))
dnds_all_LFS
### plot 
patient_labels <- c("Family member A" = "LFS02", "Family member C" = "LFS03", "Patient" = "LFS01")
pd <- position_dodge(width = 0.9)
breaks <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
breaks <- c(0.25, 1, 4, 16)
dnds_classic <- ggplot(dnds_all_LFS, aes(x = SYMBOL, y = dnds, color = Group)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_errorbar(data = boot_iqr,
                aes(x = SYMBOL, color = Group, ymin = dnds_q25, ymax = dnds_q75),
                position = pd, width = 0.2, inherit.aes = FALSE, na.rm = TRUE) +
  geom_point(position = pd, size = 1.5) +
  geom_text(
    aes(y = -0.5, label = label_fraction, group = Group),
    position = position_dodge(width = 0.9),
    vjust = 1, hjust = 0.5,
    size = 8*25.4/72.27,
    lineheight = 0.5,
    color = 'black'
  ) +
  scale_y_continuous(limits = c(-6, 76), breaks = c(0,1,2,3,4,5), labels = c(0,1,2,3,4,5)) +
  scale_color_manual(values = group_colors) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=8, margin=margin(0,0,0,0)),
    axis.text.y = element_text(size=8, margin=margin(0,-2,0,-2)),
    axis.title.y = element_text(size=8),
    legend.title = element_blank(), 
    axis.title.x = element_blank(),
    legend.position = c(0.56,0.85),
    legend.key.size = unit(8, "pt"),
    plot.margin = margin(l=2,r=0,b=0,t=0),
    legend.text = element_text(size=8)
  ) + 
  coord_cartesian(ylim = c(-1.6, 5.5), clip = "off") +
  annotate(
    "text",
    x = -Inf,             # left side
    y = -0.5,
    label = "nonsyn\n—\nsyn",
    hjust = 0.8, vjust = 1,
    size = 8*25.4/72.27,
    lineheight = 0.5
  )
dnds_classic


#ggsave("results/dnds_chip_ms.png", dnds_classic, width = 3.5, height = 1.5, units = "in", dpi = 300)
file_out = "chip_dnds_all_genes.png"
ggsave(paste0("results/Manuscript_figures/Fig_2/", file_out), dnds_classic, width = 3.5, height = 1.5, units = "in", dpi = 300)

