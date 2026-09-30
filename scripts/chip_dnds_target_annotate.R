# Copyright (c) 2026 Hunter Colegrove, Alison Feder, and University of Washington

## Reformats the VEP annotated all mutations file and counts syn/non-syn
## Run chip_dnds_target_sites.R
## Use Ensembl VEP webtool to annotate consequences


# rank order consequence types from:
# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
consequence_types <- c(
  "splice_acceptor_variant",
  "splice_donor_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "missense_variant",
  "splice_donor_5th_base_variant",
  "splice_region_variant",
  "splice_donor_region_variant",
  "splice_polypyrimidine_tract_variant",
  "stop_retained_variant",
  "synonymous_variant",
  "mature_miRNA_variant",
  "5_prime_UTR_variant",
  "3_prime_UTR_variant",
  "non_coding_transcript_exon_variant",
  "intron_variant",
  "NMD_transcript_variant",
  "non_coding_transcript_variant",
  "upstream_gene_variant",
  "downstream_gene_variant",
  "regulatory_region_variant",
  "intergenic_variant"
)

# function to pick highest ranking consequence
rank_consequence <- function(consequence_string){
  cons <- unlist(strsplit(consequence_string, ","))
  cons <- trimws(cons)
  idx <- match(cons, consequence_types)
  final <- cons[which.min(idx)]
  return(final[1])
}

all_muts_vep_path <- file.path("inputs/all_muts_vep.txt.gz")

all_muts_vep <- read_delim(all_muts_vep_path)
all_muts_vep_clean <- all_muts_vep %>%
  separate(Location, into = c("Chromosome", "Position"), sep = ":", remove = FALSE) %>%
  mutate(Mutation_ID = `#Uploaded_variation`) %>%
  dplyr::select(Chromosome, Position, Mutation_ID, Consequence, SYMBOL)

rank_muts_vep <- all_muts_vep_clean %>%
  mutate(rank_consequence = sapply(Consequence, rank_consequence)) 

unique(rank_muts_vep$rank_consequence)


## group mutations into nonsense, missense, synonymous, and splice
## remove other mutations
mutation_groups <- c(
                     #nonsense
                     "stop_gained" = "nonsense", 
                     "start_lost" = "nonsense",
                     "stop_lost" = "nonsense",
                     
                     #missense
                     "missense_variant" = "missense",
                     
                     #synonymous
                     "synonymous_variant" = "synonymous",
                     "stop_retained_variant" = "synonymous",
                     
                     #splice
                     "splice_donor_5th_base_variant" = "splice",
                     "splice_region_variant" = "splice",
                     "splice_acceptor_variant" = "splice",
                     "splice_donor_variant" = "splice",
                     "splice_polypyrimidine_tract_variant" = "splice", 
                     "splice_donor_region_variant" = "splice")

rank_muts_vep_annotate <- rank_muts_vep %>%
  mutate(impact = mutation_groups[rank_consequence]) %>%
  filter(!is.na(impact))

## save file
all_possible_muts_annotated_path <- file.path("inputs/all_muts_annotated.tsv.gz")
write_delim(rank_muts_vep_annotate, all_possible_muts_annotated_path, delim="\t")

