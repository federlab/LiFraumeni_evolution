# Copyright (c) Brendan Kohrn and University of Washington

DP_table_full <- 
  read_delim(inputs$depth_file, delim="\t", col_types = cols(.default = "c" )) %>% 
  type_convert() %>% 
  mutate(NF = Ns/(DP))

# Mutagenesis bed
if (! file.exists("inputs/bed_file_MUT.csv") ) {
  testBed_MUT <- DP_table_full %>%
    dplyr::select(Chr, Pos) %>%
    unique()
  ## Check for in-panel
  testBed_MUT$InBed <- check_bed_overlap(targ_bed6_MUT,
                                          testBed_MUT$Chr,
                                          testBed_MUT$Pos)
  ## Check for region
  testBed_MUT$Gene <- check_bed_overlap(targ_bed12_MUT,
                                         testBed_MUT$Chr,
                                         testBed_MUT$Pos,
                                         returnGenes = TRUE)
  
  # Save testBed
  testBed_MUT %>% write_csv("inputs/bed_file_MUT.csv")
} else {
  testBed_MUT <- read_csv("inputs/bed_file_MUT.csv")
}

DP_table <- DP_table_full %>%
    left_join(testBed_MUT)
