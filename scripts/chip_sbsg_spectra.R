# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

##### plot the spectra of SBSG for figure 2G

spectra <- read_delim("inputs/sbs_spectra.txt")
spectra
sbsG_data <- spectra %>%
  dplyr::select(Type, SBSG) %>%
  mutate(
    Substitution = str_extract(Type, "(?<=\\[)[ACGT]>[ACGT](?=\\])"),
    Left  = str_sub(Type, 1, 1),
    Right = str_sub(Type, -1, -1),
    Context = paste0(Left, str_sub(Substitution, 1, 1), Right)  # e.g. ACA
  )

substitution_colors <- c(
  "C>A" = "#3DB5E8", 
  "C>G" = "#000000", 
  "C>T" = "#E32929",  
  "T>A" = "#CBC9CF", 
  "T>C" = "#64BC47", 
  "T>G" = "#ED9ACB"  
)

mutation_order <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

sbsG_data <- sbsG_data %>%
  mutate(Substitution = factor(Substitution, levels = mutation_order)) %>%
  group_by(Substitution) %>%
  mutate(order = row_number()) %>%
  ungroup()

sbsG <- ggplot(sbsG_data, aes(x = Context, y = SBSG, fill = Substitution)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.05) +
  scale_fill_manual(values = substitution_colors) +
  facet_grid(~ Substitution, scales = "free_x", space = "free_x") +
  labs(y = "SBSG mutation\nprobabilities") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "pt"),
    axis.text.x = element_text(angle = 90, size = 4, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    strip.text = element_text(size = 8, color = "black"),
    plot.margin = margin(0, 1, 0, 1, "pt")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)))



sbsG <- ggplot(sbsG_data, aes(x = Context, y = SBSG, fill = Substitution)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.05) +
  scale_fill_manual(values = substitution_colors) +
  facet_grid(~ Substitution, scales = "free_x", space = "free_x", switch = "x") +
  labs(y = "SBSG mutation probabilities") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "pt"),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, hjust = 0.7),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    strip.text = element_text(size = 8, hjust = 0.5, vjust = 1),
    plot.margin = margin(4, 1, -1, 1, "pt")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)))

ggsave("results/Manuscript_figures/Fig_2/sbsg_spectra.png", sbsG, height = 1.5, width = 1.9, units = "in", dpi = 300)
