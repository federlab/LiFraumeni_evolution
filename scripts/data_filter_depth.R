# Copyright (c) Brendan Kohrn and University of Washington

DP_filt <- 
  DP_table %>% 
  filter(NF <= inputs$maxNs) %>% 
  filter(DP >= inputs$min_depth)

## CHIP
# DP_filt_CHIP <- 
#   DP_table_CHIP %>% 
#   filter(NF <= inputs$maxNs) %>% 
#   filter(DP >= inputs$min_depth)

## CHIP + MUT
# DP_filt_CHIP_MUT <- 
#   DP_table_CHIP_MUT %>% 
#   filter(NF <= inputs$maxNs) %>% 
#   filter(DP >= inputs$min_depth)


