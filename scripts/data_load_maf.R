# Copyright (c) Brendan Kohrn and University of Washington

MAF_table <-
  loadMaf(inputs$concat_maf_file,
          alphamissense,
          intOGen
  )  %>% 
  mutate(VAF = t_alt_count / (t_ref_count + t_alt_count),
         t_NC = t_depth - t_ref_count - t_alt_count,
         NF = t_NC / (t_depth) ) %>%
  mutate(SampID1 = Tumor_Sample_Barcode,
         SampID = Tumor_Sample_Barcode,    
         PatientID = Tumor_Sample_Barcode  
  ) %>%
  left_join(sample_data, by=c("PatientID")) 
