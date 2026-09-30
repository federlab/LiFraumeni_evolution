# Copyright (c) Brendan Kohrn and University of Washington

loadBed <- function(inBedFile) {
  outBed <- read_table(inBedFile, col_names = FALSE)
  #bedColNames = c("Chrom","Start","End", "Name","Score","Strand","thickStart","thickEnd","itemRgb","blockCounts","blockSizes","blockStarts")
  bedColNames = c("Chrom","Start","End", "Name")
  numCols = length(colnames(outBed))
  colnames(outBed) <- bedColNames[1:numCols]
  outBed$Start <- outBed$Start + 1
  outBed$End <- outBed$End + 1
  outBed$Target = factor(outBed$Name, levels = c(unique(outBed$Name),"Off_Target"))
  return(outBed)
}
