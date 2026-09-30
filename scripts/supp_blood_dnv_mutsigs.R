# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Supplemental Figure S5: CC>TT dinucleotide variants and SBS signatures in blood
## Top:    DNV counts per blood sample (CC>TT vs all other DNVs)
## Middle: SBS signature proportions per blood sample
## Bottom: total SBS mutations assigned per sample
##
## Requires maf_masked_coding (annotate_coding_depth.R, cross_contamination_filter.R)
## and 00.1_SharedConstants.R (subj_abbr, tissue_abbreviations).
## SBS activities are computed outside this pipeline (SigProfilerAssignment).

blood_sig_path <- "inputs/COSMIC/external/blood/blood_assignment_activities_26Dec.txt"
sample_id_mapping_path <- "inputs/sampleID_mapping.txt"

blood_tissues <- c("PBMC", "Bone marrow", "Buffy coat", "Whole blood", "Plasma")

## x-axis order shared by all panels (ordered by subject age)
blood_sample_order <- c(
  "PBMC (CON01)",
  "PBMC (CON04)",
  "PBMC (CON03)",
  "PBMC (CON02)",
  "Buffy (LFS01)",
  "Plasma (LFS01)",
  "WB (LFS01)",
  "BM (LFS01)",
  "PBMC (CON05)",
  "PBMC (LFS02)",
  "PBMC (CON06)",
  "PBMC (REL01)",
  "PBMC (LFS03)",
  "PBMC (CON07)"
)

## bold the subject code in axis labels, e.g. "PBMC (<b>CON01</b>)"
bold_subject_label <- function(x) gsub("\\(([^)]+)\\)", "(<b>\\1</b>)", x)

blood_sig_colors <- c("SBS1"   = "#E9D09E",
                      "SBS17a" = "#1AC753",
                      "SBS17b" = "#10551A",
                      "SBS40a" = "#FD9C2F",
                      "SBS5"   = "#A2D9E0",
                      "SBS7a"  = "#E898D6",
                      "SBS7b"  = "#BF077D",
                      "SBS7d"  = "#A335C2",
                      "SBSG"   = "black")

dnv_colors_blood <- c("CC>TT" = "#E69F00",
                      "Other" = "#888888")
dnv_group_order <- names(dnv_colors_blood)


###############################################################################
### DNV counts in blood
###############################################################################

dnv_blood <- maf_masked_coding %>%
  filter(Variant_Type == "DNP") %>%
  filter(Tissue %in% blood_tissues) %>%
  mutate(
    DNV_type = paste0(Reference_Allele, ">", Tumor_Seq_Allele2),
    ## CC>TT and its reverse complement GG>AA
    DNV_group = if_else(DNV_type %in% c("CC>TT", "GG>AA"), "CC>TT", "Other")
  ) %>%
  left_join(subj_abbr, by = "Subject") %>%
  left_join(tissue_abbreviations, by = "Tissue") %>%
  mutate(
    SampleLabel = factor(paste0(Tissue_abbr, " (", Subject_abbr, ")"), levels = blood_sample_order),
    DNV_group = factor(DNV_group, levels = dnv_group_order)
  ) %>%
  count(SampleLabel, DNV_group, name = "Count") %>%
  complete(
    SampleLabel = factor(blood_sample_order, levels = blood_sample_order),
    DNV_group = factor(dnv_group_order, levels = dnv_group_order),
    fill = list(Count = 0)
  )

dnv_counts_top <- ggplot(dnv_blood, aes(x = SampleLabel, y = Count, fill = DNV_group)) +
  geom_col() +
  scale_fill_manual(values = dnv_colors_blood,
                    labels = c("CC>TT" = "CC>TT", "Other" = "Other DNVs")) +
  scale_x_discrete(drop = FALSE, labels = bold_subject_label) +
  scale_y_continuous(breaks = c(0, 4, 8), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 8), clip = "off") +
  labs(y = "DNV\ncount") +
  theme_minimal() +
  theme(
    axis.text.x.bottom = element_blank(),
    axis.title.x       = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.y       = element_text(size = 8),
    axis.text.y        = element_text(size = 8, margin = margin(r = 0)),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(8, "pt"),
    legend.key.height  = unit(8, "pt"),
    legend.key.width   = unit(8, "pt"),
    legend.position    = "right",
    legend.margin      = margin(0, 0, 0, 0),
    plot.margin        = margin(1, 1, 8, 1)
  )


###############################################################################
### SBS signature proportions in blood
###############################################################################

sample_map_blood <- read_delim(sample_id_mapping_path, delim = "\t", quote = "\"") %>%
  mutate(Tissue  = str_trim(str_replace_all(tissue, '"', '')),
         Subject = str_remove(str_trim(subject), ":$")) %>%
  dplyr::select(Sample = sample, Subject, Tissue)

blood_sigs_long <- read_delim(blood_sig_path) %>%
  dplyr::rename(Sample = Samples) %>%
  pivot_longer(cols = starts_with("SBS"), names_to = "Signature", values_to = "Count") %>%
  left_join(sample_map_blood, by = "Sample") %>%
  left_join(subj_abbr, by = "Subject") %>%
  left_join(tissue_abbreviations, by = "Tissue") %>%
  mutate(SampleLabel = factor(paste0(Tissue_abbr, " (", Subject_abbr, ")"), levels = blood_sample_order))

## every signature sample must map onto the shared x-axis
stopifnot(!any(is.na(blood_sigs_long$SampleLabel)))

blood_sigs_bottom <- ggplot(blood_sigs_long, aes(x = SampleLabel, y = Count, fill = Signature)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "SBS\nproportion") +
  scale_fill_manual(values = blood_sig_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE, labels = bold_subject_label) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x.bottom = element_markdown(angle = 90, size = 8, vjust = 0.5, hjust = 1,
                                          margin = margin(t = 0)),
    axis.title.y      = element_text(size = 8),
    axis.title.x      = element_blank(),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(8, "pt"),
    legend.key.height = unit(8, "pt"),
    legend.key.width  = unit(8, "pt"),
    axis.ticks.x      = element_blank(),
    legend.position   = "right",
    axis.text.y       = element_text(size = 8, margin = margin(r = -20)),
    plot.margin       = margin(0, 1, 0, 1),
    legend.margin     = margin(0, 0, 0, 0)
  )

## total SBS mutations per sample
blood_sig_totals <- blood_sigs_long %>%
  group_by(SampleLabel) %>%
  summarise(total_mutations = sum(Count, na.rm = TRUE), .groups = "drop")

n_panel <- ggplot(blood_sig_totals, aes(x = SampleLabel, y = 1)) +
  geom_text(aes(label = total_mutations), size = 2.5, lineheight = 0.9) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(y = "n=") +
  theme_void() +
  theme(
    axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5, margin = margin(r = -30)),
    plot.margin  = margin(-6, 1, 1, 1)
  )


###############################################################################
### Combine
###############################################################################

blood_dnv_mutsigs <- dnv_counts_top / blood_sigs_bottom / n_panel +
  plot_layout(heights = c(0.35, 0.8, 0.2), guides = "collect") &
  theme(
    legend.position      = "right",
    legend.box           = "vertical",
    legend.spacing.y     = unit(25, "pt"),
    legend.justification = c(0.5, 0.9),
    legend.box.margin    = margin(t = 0, r = 0, b = 0, l = -10)
  )

blood_dnv_mutsigs
ggsave("results/Manuscript_figures/Fig_S5/blood_dnv_mutsigs.png", blood_dnv_mutsigs, width = 5, height = 2.5, units = "in", dpi = 300)
