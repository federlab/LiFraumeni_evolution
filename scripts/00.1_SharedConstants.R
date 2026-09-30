# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Shared constants used across multiple analysis scripts.

## CHIP panel genes ####
CHIP_genes <- c("NRAS", "BRINP3", "DNMT3A", "IDH1", "GATA2", "KIT", "TET2", "NPM1",
                "EZH2", "RAD21", "HNRNPK", "PTEN", "SMC3", "WT1", "KMT2A", "CBL", "KRAS",
                "PTPN11", "FLT3", "IDH2", "MYH11", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "SMC1A",
                "STAG2", "PHF6", "TP53")

## Subject groupings ####
mstp              <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4",
                       "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")
family            <- c("Family member A", "Family member B", "Family member C")
lfs_subjects      <- c("Patient", "Family member A", "Family member C")
ctx_subjects      <- c("Patient", "UW volunteer 7")
family_patient_blood         <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat")

## Subject ages ####
age_map <- c("UW volunteer 1" = 25,
             "UW volunteer 2" = 30,
             "UW volunteer 3" = 27,
             "UW volunteer 4" = 25,
             "Patient"        = 34,
             "Family member A"= 39,
             "Family member B"= 61,
             "UW volunteer 5" = 37,
             "Family member C"= 69,
             "UW volunteer 6" = 60,
             "UW volunteer 7" = 76)

## Subject abbreviations ####
subj_abbr <- tribble(
  ~Subject_abbr, ~Subject,
  "LFS01", "Patient",
  "LFS02", "Family member A",
  "REL01", "Family member B",
  "LFS03", "Family member C",
  "CON01", "UW volunteer 1",
  "CON02", "UW volunteer 2",
  "CON03", "UW volunteer 3",
  "CON04", "UW volunteer 4",
  "CON05", "UW volunteer 5",
  "CON06", "UW volunteer 6",
  "CON07", "UW volunteer 7"
)

## Tissue ordering and groupings ####
tissue_order <- c("Whole blood",
                  "Buffy coat",
                  "Plasma",
                  "Bone marrow",
                  "Thyroid",
                  "Mainstem bronchus",
                  "Lung",
                  "Esophagus 1",
                  "Esophagus 2",
                  "Gastric 1",
                  "Gastric 2",
                  "Cardiac muscle",
                  "Spleen",
                  "Liver",
                  "Colon",
                  "Omentum",
                  "Peritoneum",
                  "Renal",
                  "Testis",
                  "Skeletal muscle",
                  "Skin",
                  "Skin, non-sun-exposed",
                  "Mediastinal metastasis",
                  "Lung metastasis",
                  "Esophageal cancer 1",
                  "Esophageal cancer 2",
                  "Liver metastasis 1",
                  "Liver metastasis 2")

non_cancer_samples <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow",
                        "Thyroid", "Mainstem bronchus", "Lung", "Esophagus 1", "Esophagus 2",
                        "Gastric 1", "Gastric 2", "Cardiac muscle", "Spleen", "Liver", "Colon",
                        "Omentum", "Peritoneum", "Renal", "Testis", "Skeletal muscle",
                        "Skin", "Skin, non-sun-exposed")

cancer_samples <- c("All tissue-types", "Mediastinal metastasis", "Lung metastasis",
                    "Esophageal cancer 1", "Esophageal cancer 2",
                    "Liver metastasis 1", "Liver metastasis 2")

blood_samples <- c("Whole blood", "Bone marrow", "Plasma", "Buffy coat")

## Tissue ordering variants ####
tissue_order_qc <- c(
  "Bone marrow", "Whole blood", "Plasma", "Buffy coat", "PBMC",
  "Thyroid", "Mainstem bronchus", "Lung",
  "Esophagus 1", "Esophagus 2", "Gastric 1", "Gastric 2",
  "Cardiac muscle", "Spleen", "Liver", "Colon",
  "Omentum", "Peritoneum", "Renal", "Testis",
  "Skeletal muscle", "Skin", "Skin, non-sun-exposed",
  "Mediastinal metastasis", "Lung metastasis",
  "Esophageal cancer 1", "Esophageal cancer 2",
  "Liver metastasis 1", "Liver metastasis 2"
)

tissue_order_qc_panel <- c(
  "Bone marrow", "Whole blood", "Plasma", "Buffy coat", "PBMC",
  "Esophageal cancer 1", "Esophageal cancer 2",
  "Liver metastasis 1", "Liver metastasis 2",
  "Lung metastasis", "Mediastinal metastasis",
  "Liver", "Colon", "Skin", "Skin, non-sun-exposed",
  "Esophagus 1", "Esophagus 2", "Gastric 1", "Gastric 2",
  "Thyroid", "Mainstem bronchus", "Lung",
  "Cardiac muscle", "Spleen",
  "Omentum", "Peritoneum", "Renal", "Testis",
  "Skeletal muscle"
)

tissue_order_dnds <- c(tissue_order, "All tissues")

## Tissue category groupings (used in MF and dN/dS tissue analyses) ####
tissue_categories <- tribble(
  ~Tissue,                  ~Category,
  "Whole blood",            "Blood",
  "Buffy coat",             "Blood",
  "Plasma",                 "Blood",
  "Bone marrow",            "Blood",
  "Thyroid",                "Solid",
  "Mainstem bronchus",      "Solid",
  "Lung",                   "Solid",
  "Esophagus 1",            "Solid",
  "Esophagus 2",            "Solid",
  "Gastric 1",              "Solid",
  "Gastric 2",              "Solid",
  "Cardiac muscle",         "Solid",
  "Spleen",                 "Solid",
  "Liver",                  "Solid",
  "Colon",                  "Solid",
  "Omentum",                "Solid",
  "Peritoneum",             "Solid",
  "Renal",                  "Solid",
  "Testis",                 "Solid",
  "Skeletal muscle",        "Solid",
  "Skin, non-sun-exposed",  "Solid",
  "Skin",                   "Sun-exposed skin",
  "Mediastinal metastasis", "Cancer",
  "Lung metastasis",        "Cancer",
  "Esophageal cancer 1",    "Cancer",
  "Esophageal cancer 2",    "Cancer",
  "Liver metastasis 1",     "Cancer",
  "Liver metastasis 2",     "Cancer",
  "All tissues",            "All"
)

## Tissue abbreviations (used in overlap and heatmap figures) ####
tissue_abbreviations <- tribble(
  ~Tissue,                  ~Tissue_abbr,
  "Whole blood",            "WB",
  "Buffy coat",             "Buffy",
  "Plasma",                 "Plasma",
  "Bone marrow",            "BM",
  "Thyroid",                "Thyroid",
  "Mainstem bronchus",      "Bronchus",
  "Lung",                   "Lung",
  "Esophagus 1",            "Esoph1",
  "Esophagus 2",            "Esoph2",
  "Gastric 1",              "Gast1",
  "Gastric 2",              "Gast2",
  "Cardiac muscle",         "Cardiac",
  "Spleen",                 "Spleen",
  "Liver",                  "Liver",
  "Colon",                  "Colon",
  "Omentum",                "Omentum",
  "Peritoneum",             "Peritoneum",
  "Renal",                  "Renal",
  "Testis",                 "Testis",
  "Skeletal muscle",        "Skeletal",
  "Skin",                   "Skin",
  "Skin, non-sun-exposed",  "SkinNS",
  "Mediastinal metastasis", "Med Met",
  "Lung metastasis",        "Lung Met",
  "Esophageal cancer 1",    "Esoph Ca1",
  "Esophageal cancer 2",    "Esoph Ca2",
  "Liver metastasis 1",     "Liver Met1",
  "Liver metastasis 2",     "Liver Met2",
  "PBMC",                   "PBMC",
  "Urine cells",            "Urine"
)


## Color palettes ####
group_colors <- c("non-LFS" = "#44AA99", "LFS" = "#882255")

shape_group_colors <- c(
  "non-LFS/no-CTx" = "#44AA99",
  "non-LFS/CTx"    = "#44AA99",
  "LFS/no-CTx"     = "#882255",
  "LFS/CTx"        = "#882255"
)
