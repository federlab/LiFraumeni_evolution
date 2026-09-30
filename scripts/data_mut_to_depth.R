# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## read genome.mut files and combine to a single file

file_list <- c(
  "inputs/mafs/DevDNA1_S1.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03962_S2.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03963_S3.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03964_S4.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03965_S30_S31.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03966_S32.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03968_S5.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03969_S6.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03970_S7.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03971_S8.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03972_S9.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03973_S10.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03974_S11.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03975_S12.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03976_S13.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03977_S14.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03978_S15.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03979_S16.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03980_S17.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03981_S18.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03982_S33.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03983_S34.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03984_S19.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03985_S35.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03986_S36.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03987_S20.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03988_S21.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03989_S37.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03990_S38.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03991_S39.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03992_S40.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03993_S22.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03994_S41.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03995_S23.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03996_S24.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03997_S25.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03998_S26.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA03999_S42.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA04000_S43.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA04001_S27.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA04002_S28.1.consensus.variant-calls.genome.mut",
  "inputs/mafs/DNA04003_S29.1.consensus.variant-calls.genome.mut"
)






processed_tables <- list()
if (! file.exists("inputs/mafs/combined_genome_pipeline.tsv") ) {
  for (file_path in file_list){
    mut_genome_table <- 
      read_delim(file_path, delim="\t", col_types = cols(.default = "c" )) %>% 
      type_convert()  %>%
      mutate(ERFAP = as.character(ERFAP))
    
    # rename columns to fit with depth file 
    DP_column_table <- mut_genome_table %>%
      dplyr::rename(Chr = contig) %>%
      dplyr::rename(Pos = start) %>%
      dplyr::rename(DP = depth) %>%
      dplyr::rename(Ns = no_calls) %>%
      dplyr::rename(Samp = sample)
    
    # write to new file
    output_file_path <- gsub("\\.genome\\.mut$", ".pipeline.tsv", file_path)
    write_delim(DP_column_table, output_file_path, delim="\t")
    
    # add processed file to list
    processed_tables[[file_path]] <- DP_column_table
  }
  
  combined_table <- bind_rows(processed_tables, .id = "SourceFile")
  write_delim(combined_table, 'inputs/mafs/combined_genome_pipeline.tsv', delim = "\t")
}
