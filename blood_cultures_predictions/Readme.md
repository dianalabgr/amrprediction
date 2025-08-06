This folder contains a Snakemake workflow and Jupyter notebook for processing blood culture assemblies, performing genomic annotation, extracting k-mers, and predicting in silico antibiograms using machine learning models.
📦 Required Files

Before running the pipeline, download the following datasets from Zenodo:

    bloodcultures_assemblies_antibiograms.zip – Contains the assemblies and antibiograms.

    ProcessNewAssemblies_nedded_files.zip – Contains required annotation files and scripts.

    models_and_ShapValues.zip – Contains trained models and SHAP values.

    ResFinder supplementary file for metrics comparison.

🛠 Setup Instructions

    Unzip the required files in your working directory.

    In the folder:

~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/

 inside each of the directories:

    results_23s/

    results_5s/

    results_16s/

Create these subdirectories:

    Acinetobacter_baumannii/
    Enterobacter_cloacae/
    Enterococcus_faecium/
    Escherichia_coli/
    Klebsiella_pneumoniae/
    Pseudomonas_aeruginosa/
    Staphylococcus_aureus/

    This is required so that the rRNA analysis can save the results for each species correctly.

    Update paths in the scripts according to your local setup.

🔁 Workflow

    The Snakemake pipeline performs annotation of the assemblies and extracts relevant k-mers.

    The Jupyter Notebook (.ipynb) loads trained models and predicts resistance phenotypes from the k-mers.

    The notebook also computes evaluation metrics comparing the ML-based predictions to ResFinder-based predictions.
