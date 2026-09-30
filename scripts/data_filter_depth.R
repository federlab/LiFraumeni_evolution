# Copyright (c) Brendan Kohrn and University of Washington

DP_filt <- 
  DP_table %>% 
  filter(NF <= inputs$maxNs) %>% 
  filter(DP >= inputs$min_depth)

