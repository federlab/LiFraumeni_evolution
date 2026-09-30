# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

sampleID_map <- "inputs/sampleID_mapping.txt"
path <- "inputs/COSMIC/external/blood/"
sbs_counts <- "blood_assignment_activities_26Dec.txt"

ann_wide <- tibble::tibble(
  subject    = c("UW volunteer 1:","UW volunteer 2:","UW volunteer 3:","UW volunteer 4:",
                 "UW volunteer 5:","UW volunteer 6:","UW volunteer 7:",
                 "Patient:","Family member A:","Family member C:","Family member B:"),
  SampleCode = c("CON01","CON02","CON03","CON04","CON05","CON06","CON07",
                 "LFS01","LFS02","LFS03","REL01"))

#sample_order <- c("CON01", "CON02", "CON03", "CON04", "LFS01", "LFS02", "REL01", "CON05", "LFS03", "CON06", "CON07")
age_map <- c("CON01" = 25,
             "CON02" = 30,
             "CON03" = 27,
             "CON04" = 25,
             "LFS01" = 34,
             "LFS02" = 39,
             "REL01" = 61,
             "CON05" = 37,
             "LFS03" = 69,
             "CON06" = 60,
             "CON07" = 76)

sig_colors <- c("SBS1" = "#E9D09E", 
                "SBS17a"  = "#1AC753", 
                "SBS17b" = "#10551A", 
                "SBS40a" = "#FD9C2F", 
                "SBS5" = "#A2D9E0", 	
                "SBS7a" = "#E898D6", 
                "SBS7b" = "#BF077D", 
                "SBS7d" = "#A335C2",
                "SBSG" = "black")


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
tissue_sigs_long <- left_join(tissue_sigs_long, ann_wide, by = "subject")

tissue_sigs_long
tissue_sigs_label <- tissue_sigs_long %>%
  left_join(metadata, by = c("SampleCode" = "Subject_abbr", "tissue" = "Tissue")) %>%
  mutate(SampleLabel = paste0(Tissue_abbr, " (<b>", SampleCode, "</b>)")) %>%
  mutate(age = age_map[SampleCode]) %>%     
  arrange(age) %>%                          
  mutate(SampleLabel = factor(
    SampleLabel,
    levels = unique(SampleLabel)))            



blood_sigs <- ggplot(tissue_sigs_label, aes(x=SampleLabel, y=Count, fill = Signature)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion") + 
  scale_fill_manual(values = sig_colors) + 
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_minimal() +
  theme(axis.text.x.bottom = element_markdown(angle = 90, size = 8, vjust = 0.5, hjust = 1, margin=margin(t=0)),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=8, margin=margin(0,0,0,0)),
        legend.key.size = unit(8, "pt"),
        plot.margin = margin(1,1,1,1), 
        axis.ticks.x = element_blank(),
        # legend top 
        legend.position = "top",
        legend.justification = c(1,0.5),
        legend.margin = margin(b=-8),
        axis.text.y = element_text(size = 8, margin = margin(r=-20)))
        # legend right
        #legend.margin = margin(0,1,0,-10))
blood_sigs
ggsave("results/Manuscript_figures/Fig_2/chip_mutsigs.png", blood_sigs, width = 2, height = 1.5, units = "in", dpi = 300)
