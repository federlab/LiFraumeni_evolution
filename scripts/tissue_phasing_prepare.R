# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

#### Hunter Colegrove
#### Filter mutations in LFS carriers that are nearby the germline LFS mutation.
#### Export as a file for use in phasing_tp53_181.py to find which allele 2nd hit mutations arise on. 


close_muts_181 <- maf_masked_coding %>% 
  filter(Hugo_Symbol == "TP53") %>% 
  filter(Start_Position <= (7675069 + 150) & Start_Position >= (7675069 - 150)) %>%
  filter(Subject %in% c("Patient", "Family member A", "Family member C")) %>%
  filter(Tissue != "Urine cells") #%>%
  #filter(!(Variant_Type == "DEL" | Variant_Type == "INS")) %>%
  #print(width = Inf, n = Inf)

write.table(close_muts_181, file = "results/close_muts_181.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

## phasing of SNP
# close_muts_72 <- maf_masked_coding %>%
#   filter(Hugo_Symbol == "TP53") %>%
#   filter(Start_Position <= (7676153 + 150) & Start_Position >= (7676153 - 150)) %>%
#   filter(Subject %in% c("Patient", "Family member A", "Family member C")) %>%
#   filter(Tissue != "Urine cells") %>%
#   #filter(!(Variant_Type == "DEL" | Variant_Type == "INS")) %>%
#   print(width = Inf, n = Inf)
# 
# write.table(close_muts_72, file = "results/close_muts_72.tsv",
#             sep = "\t", row.names = FALSE, quote = FALSE)
