# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Lollipop plot

# TP53 domains
filteredDomains <- read_delim(inputs$domains_file)

subjects_blood <- c("UW volunteer 1","UW volunteer 2","UW volunteer 3","UW volunteer 4",
                    "UW volunteer 5","UW volunteer 6","UW volunteer 7",
                    "Patient","Family member A","Family member C","Family member B")
## input cosmic data
cosmic_data <- read_delim("inputs/COSMIC/TP53_mutations_07OCT25_23_08_31.csv")
cosmic_counts <- cosmic_data %>% 
  group_by(Position) %>% 
  summarise(count = sum(Count))

blood_filtering <- maf_masked_coding %>% 
  filter(Subject %in% subjects_blood) %>%
  filter(Tissue %in% c("PBMC", "Buffy coat")) %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(prot.pos)) %>% 
  filter(!is.na(exon_number)) %>%
  filter(Variant_Type=="SNP")

gene <- "TP53"
my_counts <- blood_filtering %>%
  filter(Hugo_Symbol == gene) %>%
  filter(!is.na(prot.pos)) %>%
  group_by(prot.pos) %>%
  summarize(Count = n(), .groups = "drop")


#  Prepare cosmic counts 
cosmic_counts <- cosmic_data %>%
  group_by(Position) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  dplyr::rename(prot.pos = Position)

prot_len <- max(filteredDomains$aa.length)

test_plot <- ggplot(my_counts) +
  geom_linerange(aes(x = prot.pos, ymin = 0, ymax = Count), size = 0.25) +
  geom_point(aes(x = prot.pos, y = Count), size = 0.2) +
  scale_x_continuous(limits = c(0, prot_len)) +
  scale_y_continuous(name = "Mutation\nCount", 
                     limits = c(0,4),
                     breaks = c(0, 1, 2, 3, 4)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )


color_scale <- c("black", "#56B4DF", "#E69D00", "#009E74", "#D55C00", "#0071B2")
color_scale <- c("black", "#999999", "white")

filteredDomains <- filteredDomains %>%
  mutate(Label2 = factor(Label2,
                         levels = c("TAD",
                                    "DNA-binding domain",
                                    "Tetramerisation motif")))
test_domains_plot <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = prot_len, ymin = -0.25, ymax = 0.25),
            fill = "#cccccc") +
  geom_rect(
    data = filteredDomains,
    aes(xmin = Start, xmax = End, ymin = -0.35, ymax = 0.35, fill = Label2),
    color = "black",      # <-- adds border outline
  ) +
  scale_fill_manual(values = color_scale,
                    labels = c("TAD"= "TAD",
                               "DNA-binding domain"= "DBD",
                               "Tetramerisation motif" = "Tetramerisation motif"
  )) +
  scale_x_continuous(limits = c(0, prot_len)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  )

print(test_domains_plot)

# Legend for domains
domain_legend <- get_legend(
  test_domains_plot +
    theme(
      legend.key.size = unit(8, "pt"),
      legend.text = element_text(size = 8, margin = margin(0,0,0,1)), 
      legend.box.margin = margin(0, 0, 0, 0),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.margin = margin(2,0,5,45),
      legend.key.spacing.x = unit(3, "pt")
    )
)

test_plot_2 <- ggplot(cosmic_counts) +
  geom_linerange(aes(x = prot.pos, ymin = -Count, ymax = 0), size = 0.25) +
  scale_x_continuous(limits = c(0, prot_len)) +
  scale_y_continuous(labels = abs, name = "COSMIC\nCount") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(2, 2, -2, 2), 
    axis.title.y = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8)
  )

lollipop_slide <- cowplot::plot_grid(
  test_plot,
  test_domains_plot + theme(legend.position = "none"),
  test_plot_2,
  domain_legend,
  align = 'v', axis = "lr",
  rel_heights = c(4, 0.5, 2, 1.1), ncol = 1
) +
  theme(plot.margin = margin(b = -5))

# Show final plot
show(lollipop_slide)


#ggsave("results/tp53_lollipop_blood.png", lollipop_slide, width = 3.75, height = 1.5, units = "in", dpi = 300)
ggsave("results/Manuscript_figures/Fig_3/tp53_lollipop_blood.png", lollipop_slide, width = 3.75, height = 1.5, units = "in", dpi = 300)

