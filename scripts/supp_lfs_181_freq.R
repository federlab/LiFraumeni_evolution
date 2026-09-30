# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Frequency of LFS mutation by tissue type



LFS_frequency <- MAF_table %>% 
  filter(prot.pos == "181") %>% 
  filter(Subject == "Patient") %>% 
  filter(Hugo_Symbol == "TP53") %>%
  filter(Tissue %in% tissue_order) %>%
  dplyr::select(Hugo_Symbol, Start_Position, Variant_Classification, Reference_Allele, Tumor_Seq_Allele2, t_depth, protein_variant, Tissue, VAF) %>%
  arrange(desc(VAF)) %>%
  mutate(Tissue = reorder(Tissue, VAF, FUN = max)) %>%
  mutate(Tumor = ifelse(Tissue %in% cancer_samples, "Tumor", "Non-Tumor"))
custom_label <- function(x) {
  sapply(x, function(t) {
    abbr <- tissue_abbreviations$Tissue_abbr[match(t, tissue_abbreviations$Tissue)]
    if (is.na(abbr)) abbr <- t
    if (t %in% cancer_samples) {
      paste0("<span style='color:red;'>", abbr, "</span>")
    } else {
      abbr
    }
  })
}

LFS_frequency
plot <- ggplot(LFS_frequency, aes(x=Tissue, y = VAF, color = Tumor)) + 
  geom_point() + 
  labs(y = "Germline p.R181H VAF") +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = "black") + 
  theme_minimal() + 
  scale_x_discrete(labels = custom_label) +
  scale_y_continuous(limits = c(0.4,1), breaks = seq(0,1,by = 0.1)) +
  theme(axis.text.x.bottom = element_markdown(angle = 90, size = 8, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =8),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = c(0.2,0.8)) + 
  scale_color_manual(name = NULL, values = c("Tumor" = "red", "Non-Tumor" = "black")) 
  
plot
ggsave("results/Manuscript_figures/Fig_S7/supp_lfs_181_freq.png", plot, width = 4, height = 2.5, units = "in", dpi = 300)
