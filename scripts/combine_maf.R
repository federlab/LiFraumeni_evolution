# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

### Hunter Colegrove
### 18 Sep 2024
### Li-Fraumeni
###

### combine all maf files with removed INFO contents and SV mutations 
### output a maf file containing all samples: allSamples_noINFO.maf


## path of files
dat_path <- "inputs/mafs/"

if (! file.exists("inputs/mafs/allSamples_noINFO.maf") ) {
  ## find all relevant mafs
  maf_NoINFO_list <- list.files(dat_path, pattern = "\\.variant-calls\\.noINFO\\.maf$", recursive = TRUE, full.names = TRUE)

  if (length(maf_NoINFO_list) == 0) {
    stop(
      paste0(
        "No *.variant-calls.noINFO.maf files found in ", dat_path, "\n\n",
        "Per-sample MAF files are not included in this repository. ",
        "Download them from dbGaP (phs004484.v1.p1) and place them in ", dat_path
      ),
      call. = FALSE
    )
  }

  ## combine the mafs - PHENO column has interfile discrepancy call string to not remove values
  combined_data <- maf_NoINFO_list %>%
    map_dfr(~ read_delim(.x, delim = "\t", skip = 1, col_types = cols(PHENO = col_character())))

  ## write with the leading "#version" sentinel so the output is a valid MAF
  out_filepath <- file.path(dat_path, "allSamples_noINFO.maf")
  writeLines("#version 2.4", out_filepath)
  write_delim(combined_data, out_filepath, delim="\t", append = TRUE, col_names = TRUE)
}
