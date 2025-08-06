
In the snakemake file for the assemblies from blood cultures we make the annotation and we find the kmers of the assemblies of the blood cultures. These files are saved in Zenodo folder bloodcultures_assemblies_antibiograms.zip.
In order to run, ytou need to also download from zenodo folder ProcessNewAssemblies_nedded_files.zip and make this change 
In the folder ~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/results_23s In the folder ~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/results_5s and In the folder ~/ProcessNewAssemblies_nedded_files/annotation/annotation/rRNA_genes/results_16s
create the folders 
Acinetobacter_baumannii
Enterobacter_cloacae
Enterococcus_faecium
Escherichia_coli
Klebsiella_pneumoniae
Pseudomonas_aeruginosa
Staphylococcus_aureus
in order to be able for the results of the rRNA analysis to be able to be saved there.,                    
And also in the scripts change the paths accordinggly to your setup. 

In the ipynb file we predict the in silico antibiograms using the models that we have created and also we create the metrics from the resfinder prediction using the files models_and_ShapValues.zip and the resfinder supplementary file.  
