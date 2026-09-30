# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

### Visualize sequencing depths across subjects and genes in CHIP panel
### Supplementary figure - sequencing depths by CHIP and MUT panel

sample_id_mapping_path <- "inputs/sampleID_mapping.txt"
sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))


#### all blood and tissues (exclude the following urine, buccal, and devDNA)
samples_all_exclude <- c("DNA03980_S17.1", "DNA03965_S30_S31.1", "DNA03966_S32.1", "DevDNA1_S1.1")

tissue_order <- tissue_order_qc



region_map <- tribble(
  ~chr,      ~full_region,
  "chr1",    "region_208",
  "chr2",    "region_2896",
  "chr4",    "region_4173",
  "chr6",    "region_5020",
  "chr7",    "region_5144",
  "chr8",    "region_5520",
  "chr9",    "region_5739",
  "chr10",   "region_784",
  "chr11",   "region_1111",
  "chr12",   "region_1355",
  "chr13",   "region_1501",
  "chr14",   "region_1725",
  "chr15",   "region_1904",
  "chr16",   "region_2115",
  "chr17",   "region_2378",
  "chr18",   "region_2457",
  "chr19",   "region_2739",
  "chr20",   "region_3388",
  "chr21",   "region_3515",
  "chr22",   "region_3703"
) %>%
  mutate(region = as.integer(sub("region_", "", full_region))) %>%
  dplyr::select(chr, region) %>%
  mutate(label = paste0(chr," (",region,")"))
region_list <- region_map$region


chip_depths <- final_masked_depth %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  filter(!is.na(gene_name)) %>%
  group_by(gene_name, subject, tissue, Samp, Chr) %>%
  summarise(
    meanDP_subject = mean(DP, na.rm = TRUE),
  ) %>% 
  left_join(subj_abbr, by=c("subject" = "Subject")) %>%
  left_join(tissue_abbreviations, by = c("tissue" = "Tissue")) %>%
  mutate(label = paste0(Tissue_abbr," (",Subject_abbr,")")) 


################################################################################
######## Average depth across samples by gene
################################################################################

### plot with sample labels

chip_depths_no_tp53 <- chip_depths %>% filter(gene_name != "TP53")
depth_by_gene_color <- ggplot(chip_depths_no_tp53, aes(x = gene_name, y = meanDP_subject, color = label)) +
  geom_point(size = 1) +
  theme_minimal() +
  scale_y_continuous(limits = c(0,30000)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size=8),
    legend.title = element_blank(), 
    legend.key.spacing.y = unit(-8, "pt"),
    legend.key.spacing.x = unit(-5, "pt"),
    legend.text = element_text(size = 8, margin=margin(l=-2)), 
    legend.margin = margin(l=-6),
    legend.box.spacing = unit(2, 'pt')
  ) +
  labs(y = "Mean depth (CHIP panel)")

depth_by_gene_color


######## Average depth across samples by mutagenesis region

chip_mut_depths <- final_masked_depth %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  mutate(region_split = sub("_.*", "", Gene),
         region = ifelse(region_split == "region",
                                sub(".*_", "", Gene), NA)
  ) %>%
  filter(region_split == "region") %>%
  group_by(region, Samp, tissue, subject) %>%
  summarise(
    meanDP_subject = mean(DP, na.rm = TRUE),
  ) %>% 
  left_join(subj_abbr, by=c("subject" = "Subject")) %>%
  left_join(tissue_abbreviations, by = c("tissue" = "Tissue")) %>%
  mutate(label_sample = paste0(Tissue_abbr," (",Subject_abbr,")"))  


chip_mut_depths <- chip_mut_depths %>% 
  mutate(region = as.integer(region)) %>% 
  left_join(region_map) 

chr_order <- paste0("chr", c(1:22, "X", "Y"))
chip_mut_depths <- chip_mut_depths %>%
  mutate(label = factor(
    label,
    levels = region_map %>%
      filter(chr %in% chr_order) %>%
      arrange(factor(chr, levels = chr_order)) %>%
      transmute(label = paste0(chr, " (", region, ")")) %>%
      pull(label)
  ))


depth_by_region <- ggplot(chip_mut_depths, aes(x = label, y = meanDP_subject, color = label_sample)) +
  geom_point(size = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size=8),
    legend.title = element_blank(), 
    legend.key.spacing.y = unit(-8, "pt"),
    legend.key.spacing.x = unit(-5, "pt"),
    legend.text = element_text(size = 8, margin=margin(l=-2)), 
    legend.margin = margin(l=-6),
    legend.box.spacing = unit(2, 'pt')
  ) +
  guides(color = guide_legend(ncol = 1)) +
  labs(y = "Mean depth (MUT panel)")

depth_by_region

depth_by_gene_color <- depth_by_gene_color + theme(legend.position = "none")
plot_combine <- (depth_by_gene_color / depth_by_region) +
  plot_layout(guides = "collect")

plot_combine
ggsave("results/Manuscript_figures/Fig_S1/qc_depth_by_panel.png", plot_combine, width = 6, height = 4)

################ depth by individual (tp53 only)

chip_depths_subject <- final_masked_depth %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  mutate(region_split = sub("_.*", "", Gene),
         region = ifelse(region_split == "region",
                         sub(".*_", "", Gene), NA)
  ) %>%
  filter((region_split == "region") | (!is.na(gene_name))) %>%
  mutate(plot_region = if_else(is.na(gene_name), region, gene_name)) %>%
  group_by(plot_region, Samp) %>%
  summarise(
    meanDP_subject = mean(DP, na.rm = TRUE),
    
    total_DP = sum(DP)
    
  ) %>%   
  left_join(sample_map, by=c(Samp = "sample")) %>% 
  left_join(subj_abbr, by = c(subject = "Subject")) %>% 
  left_join(tissue_abbreviations, by=c(tissue = "Tissue")) %>%
  mutate(tissue = factor(tissue, levels = tissue_order)) %>%
  mutate(label = paste0(Tissue_abbr," (",Subject_abbr,")")) 
  
tissue_order_df <- tibble(tissue = tissue_order, 
                          OrderID = seq_along(tissue_order))

chip_depths_subject <- chip_depths_subject %>%
  left_join(tissue_order_df, by="tissue")
blood_samps <- c("PBMC","Buffy coat")
chip_depths_subject <- chip_depths_subject %>%
  mutate(
    TissueGroup = if_else(tissue %in% blood_samps | grepl("PBMC", tissue, ignore.case = TRUE),
                          "Blood",
                          "Solid tissues"),
    SampleLabel = paste0(Tissue_abbr, " (<b>", Subject_abbr, "</b>)"),
    SampleLabel = fct_reorder(SampleLabel, OrderID, .fun = min) 
  )

chip_depths_subject_tp53 <- chip_depths_subject %>% filter(plot_region == "TP53") 

depth_per_subject <- ggplot(chip_depths_subject_tp53, aes(x = SampleLabel, y = meanDP_subject)) +
  geom_point(size=1) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000)) +
  theme(
    axis.text.x.bottom = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size =8),
    axis.title.x = element_blank(),
    legend.text = element_markdown(size =8), 
    legend.position = "none"
  ) +
  labs(y = "TP53 mean depth")

depth_per_subject
ggsave("results/Manuscript_figures/Fig_S1/depth_by_subject_chip.png", depth_per_subject, width = 6, height = 2)




#############################################################################
####### Save depth information as supplemental data table
#############################################################################


gene_map <- final_masked_depth %>% filter(!is.na(gene_name)) %>% distinct(gene_name, Chr)
chip_panel <- c("Whole blood", "Buffy coat", "PBMC", "Plasma")


chip_depths_subject_panelID <- chip_depths_subject %>%
  mutate(Panel = case_when(
    TissueGroup == "Blood" & (!(plot_region %in% region_list)) ~ "CHIP",
    TissueGroup == "Solid tissues" & (!(plot_region %in% region_list)) ~ "TP53",
    plot_region %in% region_list ~ "MUT"))
chip_depths_subject_panelID

#gene_chr <- chip_depths %>% ungroup() %>% distinct(gene_name, Chr)
depth_supp_table <- chip_depths_subject_panelID %>% 
  ungroup() %>%
  #dplyr::select(Subject_abbr, Tissue_abbr, Panel, chr, target, meanDP_subject) %>%
  dplyr::select(Subject_abbr, Tissue_abbr, Panel, plot_region, meanDP_subject) %>%
  dplyr::rename(Subject = Subject_abbr,
         Tissue = Tissue_abbr,
         Target = plot_region,
         meanDP = meanDP_subject) %>% 
  left_join(region_map %>% dplyr::select(chr, region) %>% mutate(region = as.character(region)), by = c(Target = "region")) %>%
  dplyr::rename(Chr = 'chr') %>%
  dplyr::select(Subject, Tissue, Panel, Chr, Target, meanDP) %>%
  left_join(gene_map, by = c(Target = "gene_name")) %>%
  mutate(Chr = coalesce(Chr.x, Chr.y)) %>%
  dplyr::select(-Chr.x, -Chr.y) %>%
  dplyr::select(Subject, Tissue, Panel, Chr, Target, meanDP) %>%
    mutate(Target = if_else(grepl("^[0-9]", Target),
                            paste0("region_", Target),
                            Target)) %>% 
  dplyr::rename(MeanDepth = "meanDP")

out_path <- "results/Manuscript_tables/"
out_file <- "supplemental_depth_panels.csv"
if (! dir.exists(out_path) ) {
  dir.create(out_path, recursive = TRUE)
}
write_csv(depth_supp_table, paste0(out_path,out_file))


###########
#### All mutations - supplement
###########

out_file_tp53 <- "supplemental_tp53_mutations.csv"
out_file_chip <- "supplemental_chip_mutations.csv"
out_file_mut <- "supplemental_mut_mutations.csv"

columns_select <- c("Hugo_Symbol", "NCBI_Build", "Chromosome", 
                    "Start_Position", "End_Position", "Variant_Classification",
                    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", 
                    "HGVSc", "HGVSp_Short", "prot.pos", "Exon_Number", 
                    "t_depth", "t_ref_count", "t_alt_count", "t_NC", 
                    "VAF", "NF",  "IMPACT", "FILTER", "Transcript_ID", 
                    "Consequence", "Subject_abbr", "Tissue_abbr")
columns_order <- c("Subject", "Tissue", "Hugo_Symbol", "NCBI_Build", 
                   "Chromosome", "Start_Position", "End_Position", 
                   "Variant_Classification", "Reference_Allele", 
                   "Tumor_Seq_Allele2", "HGVSc", "HGVSp_Short", "prot.pos", 
                   "Exon_Number", "t_depth", "t_ref_count", "t_alt_count", 
                   "t_NC", "VAF", "NF")
samples_all_exclude
supp_mutations_list <- maf_masked_coding %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  left_join(subj_abbr) %>%
  left_join(tissue_abbreviations) %>% 
  dplyr::select(all_of(columns_select)) %>%
  dplyr::rename(Subject = Subject_abbr, 
                Tissue = Tissue_abbr) %>%
  dplyr::select(all_of(columns_order)) 

tp53_supp_mutations <- supp_mutations_list %>%
  filter(Hugo_Symbol == "TP53") 

CHIP_supp_mutations <- supp_mutations_list %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>% 
  filter(Hugo_Symbol != "TP53")

MUT_supp_mutations <- maf_masked_coding %>% 
  filter(!(Samp %in% samples_all_exclude)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  left_join(subj_abbr) %>%
  left_join(tissue_abbreviations) %>% 
  left_join(testBed_MUT %>% dplyr::select(Chr, Pos, Gene) %>% dplyr::rename(Region = Gene), by = c(Chromosome = "Chr", Start_Position = "Pos")) %>% 
  dplyr::select(all_of(c(columns_select, "Region"))) %>%
  dplyr::rename(Subject = Subject_abbr, 
                Tissue = Tissue_abbr) %>%
  dplyr::select(all_of(c(columns_order, "Region"))) %>%
  filter(grepl("^region_", Region)) 
  
write_csv(tp53_supp_mutations, paste0(out_path,out_file_tp53))
write_csv(CHIP_supp_mutations, paste0(out_path,out_file_chip))
write_csv(MUT_supp_mutations, paste0(out_path,out_file_mut))


  