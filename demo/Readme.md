# 🔬 Demo: Predicting the Antibiogram of a New Genome Assembly

This demo walks you through a complete example using a public **Escherichia coli** genome to predict its **in silico antibiogram**, using the pipelines and models included in this repository.

---

## 1. Download the Repository & Set Up Environment

### Clone the Repository

```bash
git clone https://github.com/dianalabgr/amrprediction.git
cd amrprediction/Predict_NewAssemblies
```
Create a Python Environment (Recommended)
```bash
conda create -n amr_prediction_env python=3.12.7 
conda activate amr_prediction_env 
pip install -r ../requirements/software_requirements.txt
```
---

## 2. Download and Prepare Inputs
Download Genome Assembly (FASTA)

- Go to the [NCBI genome browser](https://www.ncbi.nlm.nih.gov/datasets/genome/?biosample=SAMN14843973)
- Download the genome in **FASTA format**
- Move and unzip GCA_016987635.1_PDT000967611.1_genomic.fna into the annotate_newAssemblies/ directory.

Edit the file Escherichia_files.txt in annotate_newAssemblies/ to contain only:

GCA_016987635.1_PDT000967611.1_genomic.fna

or if you use a different genome, the name of your genome

    ⚠️ If using a genome from a different genus, you should also modify line 14 of Snakefile_annotation_newAssemblies accordingly.
---

## 3. Prepare Annotation Resources
#### Download Dependencies

-Download [`ProcessNewAssemblies_nedded_files.zip`](https://zenodo.org/record/16213507)

-Unzip it into: annotate_newAssemblies/

#### Create Required Subfolders

#### Navigate to: ~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/

Inside each of the following folders:

    results_16s/

    results_23s/

    results_5s/

#### Create the following subfolders:
```bash
mkdir Acinetobacter_baumannii
mkdir Enterobacter_cloacae
mkdir Enterococcus_faecium
mkdir Escherichia_coli
mkdir Klebsiella_pneumoniae
mkdir Pseudomonas_aeruginosa
mkdir Staphylococcus_aureus
```

✅ You must create all 7 species folders inside each results_* folder.

---

## 4. Modify Paths in Scripts
#### Update k-mer Script Paths

Open:
```bash
~/ProcessNewAssemblies_nedded_files/scripts/kmer_find_newAssemblies.py
```
At line 172, update the paths:
```python
kmer_vocab_dna = np.load('<your_path>/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flatteneddna_stats_.npy', allow_pickle=True)
kmer_vocab_protein = np.load('<your_path>/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedprotein_stats_.npy', allow_pickle=True)
kmer_vocab_rrna = np.load('<your_path>/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedrrna_stats_.npy', allow_pickle=True)
```

#### Update Snakemake File

-Open the snakemake file annotate_newAssemblies/Snakefile_annotation_newAssemblies

-Replace all instances of ProcessNewAssemblies_nedded_files with your full local path, for example:
```swift
/home/youruser/path/to/ProcessNewAssemblies_nedded_files
```
---

## 5. Run Annotation Pipeline
```bash
cd annotate_newAssemblies/
snakemake --snakefile Snakefile_annotation_newAssemblies --cores 4
```

This will annotate the genome and extract k-mer features for prediction.

---

## 6. Predict the Antibiogram
#### Download Trained Models

 -Download [`models_and_shap_values.zip](https://zenodo.org/record/16213507)
 -Unzip into: Predict_newAssemblies/predict_newAssemblies/

#### Run the Notebook
```bash
jupyter notebook
```

#### Then open and run:
```
predict_newAssemblies.ipynb
```

This will:

  -Load the trained models

  -Predict antibiotic resistance

  -Display the in silico antibiogram

---

## Output

You will get a resistance profile prediction for the input genome.
