
# 🔬 Demo: Predicting the Antibiogram of a New Genome Assembly

This demo walks you through a complete use case: predicting the **in silico antibiogram** of a publicly available genome using the pipelines provided in this repository.

---

## 🧩 Step-by-Step Instructions

### 1. Clone the Pipeline

Download or clone the `annotate_newAssemblies/` folder from this repository.

### 2. Download the Genome (Escherichia coli)

- Visit the [NCBI Genome Browser](https://www.ncbi.nlm.nih.gov/datasets/genome/?biosample=SAMN14843973)
- Download the genome in **FASTA format**
- Extract the file `GCA_016987635.1_PDT000967611.1_genomic.fna` into the `annotate_newAssemblies/` folder

### 3. Create a File List

In `annotate_newAssemblies/`, create a file named `Escherichia_files.txt` with the following content:

GCA_016987635.1_PDT000967611.1_genomic.fna


---

## 🛠 Set Up the Environment

### 4. Download Required Annotation Files

- Download [`ProcessNewAssemblies_nedded_files.zip`](https://zenodo.org/record/16213507) from Zenodo
- Unzip it to your working directory

### 5. Create rRNA Output Folders

Inside the following three directories:

~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/
├── results_16s/
├── results_23s/
└── results_5s/


Create these subfolders inside **each** of them:

Acinetobacter_baumannii/
Enterobacter_cloacae/
Enterococcus_faecium/
Escherichia_coli/
Klebsiella_pneumoniae/
Pseudomonas_aeruginosa/
Staphylococcus_aureus/


### 6. Update the Script Path

Edit the file:

~/ProcessNewAssemblies_nedded_files/scripts/kmer_find_newAssemblies.py


At line **172**, replace the default paths with:

```python
kmer_vocab_dna = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flatteneddna_stats_.npy', allow_pickle=True)
kmer_vocab_protein = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedprotein_stats_.npy', allow_pickle=True)
kmer_vocab_rrna = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedrrna_stats_.npy', allow_pickle=True)

🚀 Run the Annotation Pipeline

cd annotate_newAssemblies/
snakemake --snakefile Snakefile_annotation_newAssemblies --cores 4

This will annotate the genome and extract k-mers required for prediction.
🤖 Predict the Antibiogram
7. Load Pretrained Models

    Download models_and_shap_values.zip

    Unzip it inside the predict_newAssemblies/ folder

8. Run the Jupyter Notebook

jupyter notebook

    Open the notebook in predict_newAssemblies/

    Run all cells to:

        Load trained models

        Predict resistance profiles

        Output the in silico antibiogram
