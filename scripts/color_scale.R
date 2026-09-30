# Copyright (c) Brendan Kohrn and University of Washington

Col.amClass.breaks <- c("likely_benign","likely_benign_LC",
                        "ambiguous","ambiguous_LC",
                        "likely_pathogenic","likely_pathogenic_LC")
Col.amClass.values <- c(
  "#DAE1F1",
  "#92A8D6",
  "#FAEEC1",
  "#F4D870",
  "#EBB28A",
  "#E08143"
)

Col.amClas.labels <- c(
  "Likely benign",
  "Likely benign large clone",
  "Ambiguous",
  "Ambiguous large clone",
  "Likely pathogenic",
  "Likely pathogenic large clone"
)

Col.amClass.fill <- scale_fill_manual(
  breaks = Col.amClass.breaks,
  values = Col.amClass.values,
  labels = Col.amClas.labels,
  na.value = NA,
  name = "AlphaMissense pathogenicity"
)
