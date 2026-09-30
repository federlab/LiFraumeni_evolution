# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

sampleID_map <- "inputs/sampleID_mapping.txt"
path <- "inputs/COSMIC/external/tissues/"
sbs_counts <- "tissue_assignment_activities_26Dec.txt"


mutagenesis_tissues <- c("Liver",
                         "Colon",
                         "Skin",
                         "SkinNS",
                         "Esoph1",
                         "Esoph Ca1", 
                         "Esoph Ca2",
                         "Med Met",
                         "Lung Met",
                         "Liver Met1",
                         "Liver Met2")

tissue_abbreviations
cancer_samples_abbr <- c("Med Met", "Esoph Ca1", "Esoph Ca2", "Liver Met1", "Liver Met2", "Lung Met")
custom_label <- function(x) {
  sapply(x, function(t) {
    #abbrev <- tissue_labels[[t]]
    if (t %in% cancer_samples_abbr) {
      paste0("<span style='color:red;'>", t, "</span>")
    } else {
      t
    }
  })
}

sig_colors <- c("SBS1" = "#E9D09E", 
                "SBS17a"  = "#1AC753", 
                "SBS17b" = "#10551A", 
                "SBS40a" = "#FD9C2F", 
                "SBS5" = "#A2D9E0", 	
                "SBS7a" = "#E898D6", 
                "SBS7b" = "#BF077D", 
                "SBS7d" = "#A335C2")

tissue_sigs <- read_delim(paste0(path,sbs_counts))
tissue_sigs_long <- tissue_sigs %>% pivot_longer(
  cols = starts_with("SBS"), names_to = "Signature", values_to = "Count"
) %>%
  dplyr::rename(Sample = Samples)

sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)
sampleID <- sampleID %>%
  dplyr::rename(Sample = sample)

tissue_sigs_long <- left_join(tissue_sigs_long, sampleID, by = "Sample")
tissue_sigs_long <- tissue_sigs_long %>% left_join(tissue_abbreviations, by = c("tissue" = "Tissue"))

tissue_sigs_long <- tissue_sigs_long %>%
  filter(Tissue_abbr %in% mutagenesis_tissues) %>%
  mutate(Tissue_abbr = factor(Tissue_abbr, levels = mutagenesis_tissues))
  

tissue_sigs_long

tissue_sigs <- ggplot(tissue_sigs_long, aes(x=Tissue_abbr, y=Count, fill = Signature)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion") + 
  scale_x_discrete(labels = custom_label, expand = c(0,0)) +
  scale_fill_manual(values = sig_colors) + 
  theme_minimal() +
  theme(axis.text.x.bottom = element_markdown(angle = 90, size = 8, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=8, margin=margin(0,0,0,0)),
        legend.key.size = unit(8, "pt"),
        plot.margin = margin(2,5,2,5), 
        legend.margin = margin(0,1,0,-10),
        axis.ticks.x = element_blank())

ggsave("results/Manuscript_figures/Fig_4/tissue_mutsigs.png",tissue_sigs, width = 3.75, height = 1.5, units = "in", dpi = 300)
