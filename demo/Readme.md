# 🔬 Demo: Predicting the Antibiogram of a New Genome Assembly

This demo walks you through a complete example using a public **Escherichia coli** genome to predict its **in silico antibiogram** using the pipelines and models in this repository.

---

## 1. Download and Prepare Inputs
### Downlaod the github repo and specifically the folder Predict_NewAssemblies which will be your working directory 

### Download Genome Assembly (FASTA)

- Go to the [NCBI genome browser](https://www.ncbi.nlm.nih.gov/datasets/genome/?biosample=SAMN14843973)
- Download the genome in **FASTA format**
- Unzip `GCA_016987635.1_PDT000967611.1_genomic.fna` into the folder:  
  `annotate_newAssemblies/`

### Create Input File List

In `annotate_newAssemblies/`, change the file called `Escherichia_files.txt` to have only the following line (or the name of your assembly that you want to annotate):

GCA_016987635.1_PDT000967611.1_genomic.fna


---

## 2. Prepare Annotation Resources

### Download Dependencies

- Download [`ProcessNewAssemblies_nedded_files.zip`](https://zenodo.org/record/16213507)
- Unzip it inside the annotate_newAssemblies directory

### Create Required Subfolderss

Go to:
~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/


Inside **each of the following folders**:

results_16s/
results_23s/
results_5s/


Create these subfolders **(inside each folder)**:

-Acinetobacter_baumannii/

-Enterobacter_cloacae/

-Enterococcus_faecium/

-Escherichia_coli/

-Klebsiella_pneumoniae/

-Pseudomonas_aeruginosa/

-Staphylococcus_aureus/

bash 

mkdir Acinetobacter_baumannii
mkdir Enterobacter_cloacae
mkdir Enterococcus_faecium
mkdir Escherichia_coli
mkdir Klebsiella_pneumoniae
mkdir Pseudomonas_aeruginosa
mkdir Staphylococcus_aureus


✅ You must create **all 7 species folders** inside **each `results_*` folder**.

---

## 3. Modify the k-mer Script

Open:

~/ProcessNewAssemblies_nedded_files/scripts/kmer_find_newAssemblies.py


At **line 172**, replace the default hardcoded paths with your correct path that specifies where the folder ProcessNewAssemblies_nedded_files is located:

```python
kmer_vocab_dna = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flatteneddna_stats_.npy', allow_pickle=True)
kmer_vocab_protein = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedprotein_stats_.npy', allow_pickle=True)
kmer_vocab_rrna = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedrrna_stats_.npy', allow_pickle=True)
```

---

## 4. Modify the snakemake file 

Open:

~/annotate_newAssemblies/Snakefile_annotation_newAssemblies

Replace the path for ProcessNewAssemblies_nedded_files to the full path in your local machine (e.g. my full path is /home/argis/Desktop/pasteur/amr_prediction_project/paper/github/Demo/ProcessNewAssemblies_nedded_files)

---

## 4. Run Annotation Pipeline

Open a terminal and run:

cd annotate_newAssemblies/

snakemake --snakefile Snakefile_annotation_newAssemblies --cores 4

This will annotate the genome and extract k-mer features for prediction.

---

## 5. Predict the Antibiogram
### Download Models

    Download models_and_shap_values.zip

    Unzip it into the folder: predict_newAssemblies/

📓 Run the Notebook

jupyter notebook

    Open the notebook in predict_newAssemblies/

    Run all cells to:

        Load the trained models

        Predict antibiotic resistance

        View the in silico antibiogram

✅ Output

You will get:

    A resistance profile prediction for the input genome

 
