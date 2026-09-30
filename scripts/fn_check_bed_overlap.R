# Copyright (c) Brendan Kohrn and University of Washington

check_bed_overlap <- function(inBed, inMutChrom, inMutPos, returnGenes = FALSE) {
  out_vect =c()
  genes_vect = c()
  for (mutIter in seq_along(inMutChrom)) {
    mutInBed = FALSE
    tmp_genes_vect = c()
    for (bedIter in seq_along(inBed$Name)) {
      bedRow <- inBed[bedIter,]
      if (inMutChrom[mutIter] == bedRow$Chrom && 
          inMutPos[mutIter] >= bedRow$Start && 
          inMutPos[mutIter] < bedRow$End) {
        mutInBed = TRUE
        if (returnGenes) {
          tmp_genes_vect = c(tmp_genes_vect, bedRow$Name)
        }
      }
    }
    num_bed_lines = length(tmp_genes_vect)
    if (num_bed_lines == 0) {
      genes_vect <- c(genes_vect, "")
    } else if (num_bed_lines == 1) {
      genes_vect <- c(genes_vect, tmp_genes_vect)
    } else {
      genes_vect <- c(genes_vect, "Multiple_Genes")
    }
    out_vect <- c(out_vect, mutInBed)
  }
  if (returnGenes) {
    return(genes_vect)
  } else {
    return(out_vect)
  }
}
