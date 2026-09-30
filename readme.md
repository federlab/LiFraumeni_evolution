# Ultra-deep duplex sequencing reveals unique features of somatic evolution in the normal tissues of a family with Li-Fraumeni syndrome

This repository contains the analysis pipeline used to generate the analysis and figures in the pre-print, available at [doi.org/10.64898/2026.01.12.699071](https://doi.org/10.64898/2026.01.12.699071).

Data to run pipeline and recreate figures can be accessed through dbGaP accession phs004484.v1.p1.

Sample .BAM, .MAF, and .MUT files will be required to run analysis in entirety.

## Pipeline overview
This repository is organized as an RStudio project and is intended to be run in R with the exception of phasing analysis run as a python script.

Analyses are run using a main script (scripts/00.0_MainScript.R).  
Note: The order of script execution within the main script must be preserved, as later steps depend on objects, annotations, and intermediate results generated throughout the workflow. 

### Using the pipeline
1. Clone this repository: `git clone https://github.com/federlab/LiFraumeni_evolution.git`
2. Open `LiFraumeni_evolution.Rproj` in RStudio
3. Install required packages:
- Run `setup.R` (at the repo root) once in R to install all dependencies automatically.
- For exact package versions used in the published analysis, see [pre-print supplemental materials](https://www.biorxiv.org/content/10.64898/2026.01.12.699071v1.supplementary-material).
4. Supply required reference files:
- See following section for list of necessary files and recommended locations.
5. Configure inputs:
- Open `processing_config.txt` and update the file paths to point to your local copies of the required input files. Filtering parameters are pre-set to values used in the published analysis and can be left as-is to reproduce results.
6. Run scripts/00.0_MainScript.R
- Note: To generate the final figure, phasing_tp53_181_indels.py will need to be run in python between tissue_phasing_prepare.R and tissue_phasing.R. This script requires pandas and pysam, and consensus BAMs placed at `inputs/BAMs/<sample>/<sample>.consensus.bam`.

## Required input files
Large input files are excluded from version control and must be provided prior to running pipeline.
Below are suggested filepath locations for input files. Edit processing_config.txt if input files are located in different locations.

1. Human reference genome; inputs/ref_genomes/hg38.fa
2. Gene annotation file, inputs/gencode.v38.annotation.gff3.gz:
- https://www.gencodegenes.org/human/release_38.html
3. Alphamissense annotations file, inputs/alphamissense/AlphaMissense_hg38.tsv.gz:
- https://zenodo.org/records/8208688
- A TP53-only subset of this file must also be created and saved as inputs/alphamissense/AlphaMissense_hg38.TP53.tsv (path set by the AlphaMissense key in processing_config.txt). It contains the four header lines and all rows for the canonical TP53 transcript ENST00000445888.6:
```
cd inputs/alphamissense
gunzip -c AlphaMissense_hg38.tsv.gz | head -4 > AlphaMissense_hg38.TP53.tsv
gunzip -c AlphaMissense_hg38.tsv.gz | grep ENST00000445888 >> AlphaMissense_hg38.TP53.tsv
```
4. VEP annotations for CHIP panel target sites, inputs/all_muts_vep.txt.gz:
- Generated part-way through the pipeline: scripts/chip_dnds_target_sites.R writes results/allmutations.tsv, which must be annotated with VEP and saved to the path above before scripts/chip_dnds_target_annotate.R is run
- https://www.ensembl.org/vep
