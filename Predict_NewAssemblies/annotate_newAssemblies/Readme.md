This repository includes a Snakemake pipeline that automates the annotation of a new bacterial genome assembly and creates the necessary k-mer encodings (DNA, protein, and rRNA) for downstream antimicrobial resistance prediction.
📦 Required Files

Download and extract the following from Zenodo:

    ProcessNewAssemblies_nedded_files.zip
    Contains all required data, scripts, and resources for annotation and k-mer generation.

📁 Folder Preparation

Before running the workflow, make sure the following folders exist:
Navigate to:

~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/

Then, create the following subdirectories under each of the results_* folders (results_23s, results_5s, results_16s):

Acinetobacter_baumannii/
Enterobacter_cloacae/
Enterococcus_faecium/
Escherichia_coli/
Klebsiella_pneumoniae/
Pseudomonas_aeruginosa/
Staphylococcus_aureus/

These folders are necessary for saving the rRNA annotation results for each species.
🔧 Script Adjustment

In the script:

~/ProcessNewAssemblies_nedded_files/scripts/kmer_find_newAssemblies.py

Locate line 172 and update the hardcoded path to reflect your local setup. Example:

# Load the k-mer vocabulary for DNA, proteins, and rRNA
kmer_vocab_dna = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flatteneddna_stats_.npy', allow_pickle=True)
kmer_vocab_protein = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedprotein_stats_.npy', allow_pickle=True)
kmer_vocab_rrna = np.load('~/ProcessNewAssemblies_nedded_files/Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedrrna_stats_.npy', allow_pickle=True)

▶️ Running the Pipeline

Once everything is set up, navigate to the main folder and run:

snakemake --snakefile Snakefile_annotation_newAssemblies --cores 4

This will annotate the new assembly and generate k-mer encodings for input into the machine learning models.
