
# AMR Public Data Collection Pipeline Guide

The pipeline includes data from three main sources:

### 1. NDARO Database (National Database of Antibiotic Resistant Organisms)
- **Data Retrieved**: 28-12-2023
- **Antibiotics Count**: 107
- **Format**: NCBI Assemblies, NCBI BioSample Antibiograms
- **AMR Table Source**: [NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRtable/)

**Notes**: The NDARO AMR table can be downloaded in CSV format and requires no preprocessing before pipeline input.

### 2. BV BRC Database (Bacterial and Viral Bioinformatics Resource Center)
- **Data Retrieved**: 29-12-2023
- **Antibiotics Count**: 83
- **Format**: NCBI Assembly or SRA IDs, Patric Antibiograms
- **AMR Table Source**: [BV BRC](https://www.bv-brc.org/view/Bacteria/2#view_tab=amr)

**Preprocessing**:
  - Use the `R/data_preprocess/bv_brc_analysis.R` script.
  - This script filters to include only ESKAPEE pathogenic bacteria, assigns missing Assembly IDs, and fills in SRA IDs if Assembly IDs are missing.
  - Output: Filtered AMR table `R/final_datasets/cdc_narms_data.csv`, used as an input to the Snakemake pipeline.

### 3. CDC NARMS Database (National Antimicrobial Resistance Monitoring System)
- **Data Retrieved**: 02-01-2024
- **Format**: NCBI Assemblies, CDC NARMS Antibiograms
- **AMR Table Source**: [CDC NARMS](https://wwwn.cdc.gov/narmsnow/)

**Preprocessing**:
  - Use the `R/data_preprocess/narms_split_antibiotics.R` script.
  - This script filters data, keeps only Assembly IDs, assigns Biosample IDs if Assembly IDs are unavailable, and retrieves missing SRA IDs.
  - Output: Filtered NARMS AMR table `R/final_datasets/BV_BRC_antibiotics.csv`, used as input for data downloading in the Snakemake pipeline.


## Running the Snakemake Pipeline

Once data is downloaded and preprocessed, execute the Snakemake pipeline to download assemblyand antibiograms data. To run the Snakemake pipeline, use:

```bash
snakemake --cores <ncores> --snakefile Snakefile_milestone_1
# <ncores>: Number of cores to use
# Snakefile_milestone_1: Snakemake file to run
```

Ensure all preprocessed AMR tables are in place as specified before starting the pipeline.

## Additional Scripts
You can run the scripts in the `R/data_postprocess/` (`combine_antibiograms.R` and `combine_strains_metadata.R`) directory to create the concatenated QCed strains and antibiogram tables. These scripts are used to generate the final data tables for analysis. Also, they produce some pie plots about the distribution of genus and about the distribution of used data sources.

From the 3 folders `Assemblies/ndaro`, `Assemblies/cdc_narms` & `Assemblies/bv_brc`, copy all the files to `Assemblies/all_assemblies`. Then we will keep only the assemblies that are contained in `Assemblies/strains_metadata.csv` (N50 values cutoff) to do this we have created the script `Assemblies/keep_filtered_assemblies.sh` to filter the assemblies that are not in the final data tables. This script will move to the `Assemblies/filtered_assemblies` directory the assemblies that are in the final data tables.
