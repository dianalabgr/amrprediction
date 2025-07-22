#!/bin/bash

# --------------------------------
# Download Plus Protein Sequences
# --------------------------------
# Specify the input file containing accession numbers
input_file="protein_accessions_plus.txt"

# Create a folder to store downloaded files
output_folder="protein_accessions_plus"
mkdir -p "$output_folder"

# Loop through each accession number in the file
while IFS= read -r accession; do
    # Construct the URL for the FASTA file
    fasta_url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$accession&rettype=fasta"
    # Use wget to download the FASTA file and save it in the output folder
    wget -nc -c "$fasta_url" -O "${output_folder}/${accession}.fasta"
done < "$input_file"



# --------------------------------
# Download Core Protein Sequences
# --------------------------------
input_file="protein_accessions_core.txt"

output_folder="protein_accessions_core"
mkdir -p "$output_folder"

while IFS= read -r accession; do
    fasta_url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$accession&rettype=fasta"

    wget -nc -c "$fasta_url" -O "${output_folder}/${accession}.fasta"
done < "$input_file"
