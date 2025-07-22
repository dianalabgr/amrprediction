# AMR Public Data Collection Pipeline

Firstly, we have to download AMR data, including Whole Genome Sequencing (WGS) data and antibiograms. We prefer genome assemblies, because of their non-redundant information, so they take up less disk size. Antibiograms include antibiotics, Minimum Inhibitory Concentration (MIC) values, resistant phenotype, measure platforms and vendor.

Both of WGS and antibiograms data will be the input of a classification algorithm, for predicting MIC value for a given genome assembly, therefore also the output class (Resistant, Susceptible). In addition, the input of the algorithm will be the AMR sequences for each assembly, derived from an annotation pipeline. Below is the databases, from which we have collected the data.

### NDARO Database (National Database of Antibiotic Resistant Organisms)
---

Data Retrieved: 28-12-2023

Number of different antibiotics (including antibiotic combinations): **107**

Format:
-   NCBI Assemblies
-   NCBI BioSample Antibiograms

Filters:
- AST phenotypes: I, ND, R, S
- Enterococcus faecalis
- Enterococcus faecium
- Staphylococcus pseudintermedius
- Staphylococcus aureus
- Klebsiella pneumoniae
- Klebsiella oxytoca
- Acinetobacter baumannii
- Pseudomonas aeruginosa
- Pseudomonas putida
- Enterobacter hormaechei
- Enterobacter cloacae
- Enterobacter roggenkampii
- Enterobacter asburiae
- Enterobacter kobei
- Enterobacter ludwigii
- Enterobacter bugandensis
- E.coli and Shigella

We downloaded the NDARO AMR table as csv file from here: [NDARO](https://www.ncbi.nlm.nih.gov/pathogens/isolates/#AST_phenotypes:(*=I%20*=ND%20*=R%20*=S)%20AND%20taxgroup_name:(%22Enterococcus%20faecalis%22%20%22Enterococcus%20faecium%22%20%22Staphylococcus%20pseudintermedius%22%20%22Staphylococcus%20aureus%22%20%22Klebsiella%20pneumoniae%22%20%22Klebsiella%20oxytoca%22%20%22Acinetobacter%20baumannii%22%20%22Pseudomonas%20aeruginosa%22%20%22Pseudomonas%20putida%22%20%22Enterobacter%20hormaechei%22%20%22Enterobacter%20cloacae%22%20%22Enterobacter%20roggenkampii%22%20%22Enterobacter%20asburiae%22%20%22Enterobacter%20kobei%22%20%22Enterobacter%20ludwigii%22%20%22Enterobacter%20bugandensis%22%20%22E.coli%20and%20Shigella%22)). We have not to do any data preprocess for this file, because it is already in a good format.


### BV BRC Database (Bacterial and Viral Bioinformatics Resource Center)
---

Data Retrieved: 29-12-2023

Number of different antibiotics (including antibiotic combinations): **83**

Format:
-   NCBI Assembly or SRA IDs
-   Patric Antibiograms

Filters: None, because we downloaded the whole AMR table and we filtered it with the code of `data_acquisition/R/data preprocess/bv_brc_analysis.R` script.

We have downloaded the AMR table from here: [BV BRC](https://www.bv-brc.org/view/Bacteria/2#view_tab=amr). We have to do some data preprocess for this file, using he code of `data_acquisition/R/data preprocess/bv_brc_analysis.R` script. This script filters the data by keeping only ESKAPEE pathogenic bacteria, assigns the NCBI Assembly IDs based on the BV BRC IDs and gets the NCBI SRA IDs for the missing Assembly IDs. Then, it produces the filtered BV BRC AMR table (file `data_acquisition/R/final_datasets/bv_brc_antibiograms.csv`), which is the input of the snakemake pipeline for data downloading.


### CDC NARMS Database (National Antimicrobial Resistance Monitoring System)
---

Data Retrieved: 02-1-2024

Format:
-   NCBI Assemblies
-   CDC NARMS Antibiograms


We have downloaded the AMR table from here: [CDC NARMS](https://wwwn.cdc.gov/narmsnow/). We have to do some data transformations for this file, using the code of `data_acquisition/R/data preprocess/narms_split_antibiotics.R` script. This script filters the data by keeping only assembly IDs an if it is missing, it assigns the NCBI Biosample IDs, to get from them the SRA IDs. Then, it produces the filtered NARMS AMR table (file `data_acquisition/R/final_datasets/cdc_narms_data.csv`), which is the input of the snakemake pipeline for data downloading.

---

The next step of our analysis is to run the snakemake pipeline, which will download the data from the above databases and will convert them if needed, in order to be the input of the annotation pipeline. The pipeline is written in the `data_acquisition/Snakefile_milestone_1` file, and it is executed by the below command:
```bash
snakemake --cores ncores --snakefile Snakefile_milestone_1
        # ncores: number of cores to use
        # Snakefile_milestone_1: the snakemake file
```

---

After the data downloading and converting, we have to run the postprocessing scripts, located in the `R/data postprocess` directory. More analytically, the `data_acquisition/R/data postprocess/combine_antibiograms.R` script combines the antibiograms from the three databases and prepares the final antibiograms table. This table helps us to split it into separate tables for each antibiotic. The `data_acquisition/R/data postprocess/combine_strains_metadata.R` script creates a final table with assemblies' metadata, including the location of the isolated samples, the sequencing method, the IDs in the databases, etc.