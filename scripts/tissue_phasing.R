# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

### Hunter Colegrove
## First run phasing_tp53_181_indels.py


file <- "results/phasing_181_indels.csv"
if (!file.exists(file)) {
  stop(
    paste0(
      "Required input file not found:\n  ", file, "\n\n",
      "After running scripts/tissue_phasing_prepare.R,
      run phasing_tp53_181_indels.py in python to generate input file."
    ),
    call. = FALSE
  )
}

phased <- read_delim(file)
phased <- phased %>%
  arrange(bam_file, mut_pos, base_at_check) %>%
  distinct(ref_base, mut_base, mut_pos, bam_file, .keep_all = TRUE) %>%
  mutate(Tumor_Sample_Barcode = str_remove(bam_file, "\\.consensus\\.bam$")) %>%
  dplyr::rename(
    Chromosome = chrom,
    Start_Position = mut_pos,
    Tumor_Seq_Allele2 = mut_base,
    Phased = covers_check,
    Phased_base = base_at_check
  )

phased_short <- phased %>% 
  dplyr::select(Chromosome, Start_Position, Tumor_Seq_Allele2, Phased, Phased_base, Tumor_Sample_Barcode)

phased_df <- left_join(filt_maf, phased_short) %>%
  filter(!is.na(Phased))

total_possible <- nrow(phased_df)
phased_successful <- phased_df %>% filter(Phased == TRUE)
total_phased <- nrow(phased_successful)
base_counts <- phased_successful %>%
  count(Phased_base)

## how many muts were phased and which allele
phased_breakdown <- phased_df %>%
  mutate(Phased_allele = case_when(
    Phased_base == "C" ~ "WT",
    Phased_base == "T" ~ "LFS",
    TRUE ~ "Not phased"
  )) %>%
  count(Phased_allele) %>%
  dplyr::rename(Count = n) %>%
  mutate(Phased_allele = factor(Phased_allele, levels = c("LFS", "WT", "Not phased"))) %>%
  arrange(desc(Phased_allele)) %>%
  mutate(
    y_bottom = dplyr::lag(cumsum(Count), default = 0),
    y_top    = cumsum(Count),
    y_label  = (y_bottom + y_top) / 2,
    label_col = if_else(Phased_allele == "LFS", "white", "black")
  )

wt_n  <- phased_breakdown$Count[phased_breakdown$Phased_allele == "WT"]
lfs_n <- phased_breakdown$Count[phased_breakdown$Phased_allele == "LFS"]

bt <- binom.test(x = lfs_n, n = wt_n + lfs_n, p = 0.5, alternative = "two.sided")
p_lab <- sprintf("italic(p) == %s", sprintf("%.4f", bt$p.value))

phased_muts_181 <- ggplot(phased_breakdown, aes(x = "", fill = Phased_allele)) +
  geom_col(aes(y = Count), width = 1, color = "black", linewidth = 0.15) +
  geom_text(
    aes(y = y_label, label = Count, color = label_col),
    alpha = 1,
    size = 8 * 25.4 / 72.27
  ) +
  scale_fill_manual(
    values = c("Not phased" = "white", "WT" = "#DDDDDD", "LFS" = "black"),
    labels = c("Not phased" = "Not phased", "WT" = "Non-carrier copy", "LFS" = "Carrier copy")
  ) +
  scale_color_identity() +
  annotate(
    "text",
    x = 1.65, y = 50,
    label = p_lab,
    parse = TRUE,
    size = 8 * 25.4 / 72.27,
    hjust = 0
  ) +
  labs(y = "Count", x = "Mutations\n(+/- 150bp)", fill = "Allele") +
  coord_cartesian(ylim = c(0, 62), clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  #coord_cartesian(clip = "off") +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, margin = margin(2, 0, 0, 0), hjust = 0.1),
    axis.text.y  = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text  = element_text(size = 8, margin = margin(l = 2, b = 0)),
    legend.key.spacing.y = unit(1, "pt"),
    legend.key.size = unit(8, "pt"),
    legend.title = element_blank(),
    legend.margin = margin(0, 2, 0, -10),
    plot.margin = margin(0, 1, 1, 1),
    legend.position = "right",
    legend.justification = c(0.9,0.57)
  )

phased_muts_181

#ggsave("results/phased_mutations_181.png", phased_muts_181, width = 1.75, height = 1.5, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_4/phased_mutations_181.png", phased_muts_181, width = 1.75, height = 1.5, units = "in", dpi = 300)

