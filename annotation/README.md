# Steps for running the annotation pipeline

1. Afer we have downloaded the assemblies, we have to run the annotation pipeline to extract the appropriate sequences from them. First, we have to copy the assemblies that exist in the concatenated antibiogram file to the `assemblies` folder. We can do this by running the jupyter notebook named `scripts/copy_assemblies.ipynb`.

2. We have downloaded antimicrobial peptides from AMRFinderPlus ([link](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-05-02.2/ReferenceGeneCatalog.txt), file `annotation/amr_databases/ReferenceGeneCatalog_v3.12.txt`), CARD ([link](https://card.mcmaster.ca/download/0/broadstreet-v3.2.9.tar.bz2), file `annotation/amr_databases/aro_index_v3.8.9.tsv`) and Resfinder ([link](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/), file `annotation/amr_databases/resfinder_dna.fsa`) databases.

3. Then, we ran `annotation_CARD_AMRfinderplus.R` script to concatenate the two databases, to split core and plus proteins and to remove the duplicate ones.

4. Based on the above results, we downloaded the protein sequences using `download_proteins.sh` script.

5. For Resfinder, we downloaded from bitlocker the files of resfinder database (https://bitbucket.org/genomicepidemiology/resfinder_db.git/src) we ran `scripts/resfinder_translator.py` to translate DNA sequences to protein sequences, using the below command. We saved the results in `amr_databases/resfinder_protein.fasta` file.
```bash
python scripts/resfinder_translator.py annotation/amr_databases/resfinder_dna.fsa annotation/amr_databases/resfinder_protein.fasta
```

6. We then concatenate all proteins from protein_accessions_plus and protein_accessions_core folders, using the below bash command:
```bash
cat annotation/protein_accessions_*/* > annotation/concatenated_fasta.fasta 
```

7. We ran `scripts/append_gene_names.ipynb` to include gene names in the fasta headers.
python3 append_gene_name.py

8. We concatenated the `concatenated_fasta.fasta` file with `resfinder_protein.fasta` file, using the below bash command:
```bash
cat annotation/concatenated_fasta_appened.fasta annotation/amr_databases/resfinder_protein.fasta > annotation/all_proteins.fasta
```

9. We ran CD-HIT to remove redundant proteins and for clustering the proteins with 90% identity. We used the below command:
```bash
cd-hit -i annotation/all_proteins.fasta -o annotation/all_proteins_cdhit.fasta -c 0.9 -d 80 -M 0 -g 1 -G 0 -aL 0.9
        # -M 0     # memory limit
        # -c 0.9   # sequence identity threshold
        # -aL 0.9  # alignment coverage for the longer sequence
        # -G 0     # do not use global sequence identity
        # -g 1     # cluster sequences into the highest similarity cluster
        # -d 50    # word_length
```


10. We have ready the antimicrobial peptides. So, next step is to find matches between the antimicrobial peptides and the assemblies we downloaded, using BLASTx routine.

11. We used the results of BLASTx to find the antimicrobial peptide sequences in the assemblies and to remove noise and partially matched sequences, using specific alignment coverage and e-value thresholds.

12. We finally ran snakemake pipeline to extract antimicrobial peptides (AMPs), their upstream regions and 5s, 16s and 23s rRNA copies, using the below command:
```bash
snakemake --cores ncores --snakefile Snakefile_create_annotation_Neededfiles
        # ncores: number of cores to use
        # Snakefile_create_annotation_Neededfiles: the snakemake file
```
and
```bash
snakemake --cores ncores --snakefile Snakefile_create_annotation_Neededfiles_2
        # ncores: number of cores to use
        # Snakefile_create_annotation_Neededfiles_2: the snakemake file
```
To create all needed files.


13. We finally have the extracted antimicrobial peptides, their upstream regions and 5s, 16s and 23s rRNA copies for each assembly, by running the below command:
```bash
snakemake --cores 20 --snakefile Snakefile_annotation_assemblies
        # ncores: number of cores to use
        # Snakefile_annotation_assemblies: the snakemake file
```

Now we are ready to run the machine learning models, using the scripts and notebooks located on the `Machine_Learning` folder.
