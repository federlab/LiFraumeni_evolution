# Copyright (c) Brendan Kohrn and University of Washington

loadMaf <- function(inMafFile, 
                    AM_table, 
                    IntOGen_table) {
  cols_to_get = c("Hugo_Symbol",
                  "NCBI_Build",
                  "Chromosome",
                  "Start_Position",
                  "End_Position",
                  "Variant_Classification",
                  "Variant_Type",
                  "Reference_Allele",
                  "Tumor_Seq_Allele2",
                  "Tumor_Sample_Barcode",
                  "Samp",
                  "HGVSc",
                  "HGVSp",
                  "HGVSp_Short",
                  "Exon_Number",
                  "t_depth",
                  "t_ref_count",
                  "t_alt_count",
                  "Transcript_ID",
                  "Consequence",
                  "Existing_variation",
                  "IMPACT",
                  "FILTER",
                  "t_NC",
                  "fName")

  ## MAF files carry a leading "#version" sentinel line - skip it only when present
  maf_skip <- if (startsWith(readLines(inMafFile, n = 1), "#")) 1 else 0

  outData <- 
    read_delim(inMafFile, delim="\t", skip = maf_skip) %>%
    { print(.); . } %>%
    type_convert() %>% 
    mutate(Samp = Tumor_Sample_Barcode) %>%
    dplyr::select(any_of(cols_to_get)
    ) %>% 
    mutate(
      mutPosition = paste(
        Chromosome,
        Start_Position, 
        Reference_Allele, 
        Tumor_Seq_Allele2,
        sep=":")
    ) %>% 
    # Add variant clasifications
    left_join(variant_clasification_table, by=c("Variant_Classification")) %>%
    
    # Join with AlphaMissense data
    left_join(
      AM_table, 
      by=c(
        "Chromosome"="CHROM",
        "Start_Position" = "POS", 
        "Reference_Allele"="REF", 
        "Tumor_Seq_Allele2"="ALT"
      )
    ) %>% 
    # Finish classifying AM pathogenicity
    mutate(
      am_pathogenicity = case_when(
        !is.na(am_pathogenicity) ~ am_pathogenicity,
        Mutation_type %in% c("Nonsense_Mutation","Missense_Mutation","Indel") ~ 1,
        Variant_Classification %in% c("Splice_Region") & !is.na(Exon_Number) ~ 1,
        Variant_Classification %in% c("Splice_Region") & is.na(Exon_Number) ~ 0,
        Mutation_type == "Silent" ~ 0,
        # Mutation_type %in% c("Silent",) ~ 0,
        Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins",
                                      "Start_Loss","Splice_Site",
                                      "Nonstop_Mutation") ~ 1,
        Variant_Type %in% c("DNP","TNP","ONP") ~ 1,
        T ~ am_pathogenicity
      ),
      am_class = case_when(
        !is.na(am_class) ~ am_class,
        Mutation_type %in% c("Nonsense_Mutation","Missense_Mutation","Indel") ~ "likely_pathogenic",
        Variant_Classification %in% c("Splice_Region") & !is.na(Exon_Number) ~ "likely_pathogenic",
        Variant_Classification %in% c("Splice_Region") & is.na(Exon_Number) ~ "likely_benign",
        Mutation_type == "Silent" ~ "likely_benign",
        Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins",
                                      "Start_Loss","Splice_Site",
                                      "Nonstop_Mutation") ~ "likely_pathogenic",
        Variant_Type %in% c("DNP","TNP","ONP") ~ "likely_pathogenic",
        T ~ am_class
      )
    ) %>% 
    # Join with IntOGen data
    left_join(
      IntOGen_table, 
      by=c(
        "Chromosome" = "chr",
        "Start_Position"="pos",
        "Tumor_Seq_Allele2" = "alt"
      )
    ) %>% 
    # Finish classifying IntOGen values
    mutate(boostDM_score = case_when(
      Mutation_type %in% c("Indel") ~ 1,
      !is.na(boostDM_score) ~ boostDM_score,
      Mutation_type %in% c("Missense_Mutation","Nonsense_Mutation") ~ 1,
      Variant_Classification %in% c("Splice_Region") & !is.na(Exon_Number) ~ 1,
      Variant_Classification %in% c("Splice_Region") & is.na(Exon_Number) ~ 0,
      Mutation_type %in% c("Silent") ~ 0,
      Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins",
                                    "Start_Loss","Splice_Site",
                                    "Nonstop_Mutation") ~ 1,
      Variant_Type %in% c("DNP","TNP","ONP") ~ 1,
      T ~ boostDM_score
    ),
    boostDM_class = case_when(
      Mutation_type %in% c("Indel") ~ T,
      !is.na(boostDM_class) ~ boostDM_class,
      Mutation_type %in% c("Missense_Mutation","Nonsense_Mutation") ~ T,
      Variant_Classification %in% c("Splice_Region") & !is.na(Exon_Number) ~ T,
      Variant_Classification %in% c("Splice_Region") & is.na(Exon_Number) ~ F,
      Mutation_type %in% c("Silent") ~ F,
      Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins",
                                    "Start_Loss","Splice_Site",
                                    "Nonstop_Mutation") ~ T,
      Variant_Type %in% c("DNP","TNP","ONP") ~ T,
      T ~ boostDM_class
    )) %>% 
    
    # Join with depth file
    # left_join(bedTest %>%   
    #             mutate(coding_from_bed = if_else(InBed, "coding", "non-coding")), 
    #             by=c("Chromosome" = "Chr",
    #                "Start_Position" = "Pos")) %>%
    # left_join(bedTest %>%   
    #             mutate(coding_from_bed = if_else(InBed, "coding", "non-coding")),
    #           by=c("Chromosome" = "Chr", 
    #                "End_Position" = "Pos"), 
    #           suffix=c("_StartPosition","_EndPosition")) %>%
    
    # mutate(
    #   InBed = case_when(
    #     InBed_StartPosition == TRUE ~ TRUE,
    #     InBed_EndPosition == TRUE ~ TRUE, 
    #     TRUE ~ FALSE
    #   ),
    #   inMask = case_when(
    #     inMask_StartPosition == TRUE ~ TRUE,
    #     inMask_EndPosition == TRUE ~ TRUE, 
    #     TRUE ~ FALSE
    #   ),
    #   Gene = Gene_StartPosition,
    #   coding_from_bed = case_when(
    #     coding_from_bed_StartPosition == "coding" ~ "coding",
    #     coding_from_bed_EndPosition == "coding" ~ "coding", 
    #     TRUE ~ "non-coding"
    #   )
    # ) %>%

    # Extract protein positions
    extract(HGVSp_Short, 
            c("prot.ref","prot.pos","prot.alt"),
            "^p.([A-Z])([0-9]+)([A-Z=*]|_splice)$", 
            remove=FALSE, convert = TRUE) %>% 
    
    # # Reclassify coding
    # mutate(coding_from_maf = case_when(
    #   Variant_Classification %in% c("Splice_Region", "Splice_Site") & !is.na(Exon_Number) & !is.na(prot.pos) ~ "coding",
    #   Variant_Classification %in% c("Splice_Region", "Splice_Site") &  is.na(Exon_Number) ~ "non-coding",
    #   TRUE ~ coding
    # )) %>% 
    # dplyr::select(!c(coding)) %>% 
    # mutate(coding = coding_from_maf) %>%
    
    # mutate(
    #   isPathogenic_AM = case_when(
    #     Mutation_type == "Missense_Mutation" & am_class == "likely_pathogenic" ~ "Pathogenic",
    #     Mutation_type == "Splice" ~ "Pathogenic",
    #     Mutation_type %in% c("Indel","Nonsense_Mutation") ~ "Pathogenic",
    #     T ~ "Non-Pathogenic"),
    # 
    #   Mutation_Class = case_when(
    #     Reference_Allele == "-" ~ "Indel",
    #     Tumor_Seq_Allele2 == "-" ~ "Indel", 
    #     nchar(Reference_Allele) == 1 & nchar(Tumor_Seq_Allele2) == 1 ~ "SNP",
    #     nchar(Reference_Allele) > 1 & nchar(Reference_Allele) == nchar(Tumor_Seq_Allele2) ~ "MNP",
    #     TRUE ~ "Indel"
    #   )
    # ) 
  return(outData)
}

