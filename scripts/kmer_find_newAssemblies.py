import os
import pickle
import os
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
import argparse

def build_dataset(kmer_size, assembly, path_dna, path_protein, path_5s, path_16s, path_23s):
    """
    Build the dataset including upstream sequences and matching them with k-mer vocabularies.

    Args:
        kmer_size (int): Size of the k-mers.
        assembly (str): Name of the assembly.
        directory (str): Directory containing upstream sequence files.

    Returns:
        tuple: Contains filtered data lists X, X_dna_chars, and headers_DNA.
    """
    # Get upstream sequences as an additional feature for the model
    upstream = {}
    
    # Get upstream sequences based on the assembly
    #path = os.path.abspath(os.path.join(directory, 'upstream_seqs_final', f'{assembly}.fna'))
    
    with open(path_dna, 'r') as file:
        # Read the fasta file and store a dictionary with the assembly as the key and the sequence as the value
        fasta = {}
        for line in file:
            if line.startswith('>'):
                header = line.strip()
                fasta[header] = ''
            else:
                fasta[header] += line.strip()

        upstream[assembly] = fasta
    print("Got upstream sequences")

        # Get proteins sequences as an additional feature for the model
    proteins = {}
    
    # Get upstream sequences based on the assembly
    #path = os.path.abspath(os.path.join(directory, 'aligned_proteins_final', f'{assembly}.fna'))
    
    with open(path_protein, 'r') as file:
        # Read the fasta file and store a dictionary with the assembly as the key and the sequence as the value
        fasta = {}
        for line in file:
            if line.startswith('>'):
                header = line.strip()
                fasta[header] = ''
            else:
                fasta[header] += line.strip()

        proteins[assembly] = fasta
    print("Got proteins sequences")

        # Get 5s sequences as an additional feature for the model
    rrna_5s = {}
    
    # Get upstream sequences based on the assembly
    #path = os.path.abspath(os.path.join(directory, 'rRNAs/final_5s_seqs', f'{assembly}.fna.txt'))
    
    with open(path_5s, 'r') as file:
        rrna_5s[assembly] = file.read().strip()
    print("Got 5s rRNA sequences")

        # Get 16s sequences as an additional feature for the model
    rrna_16s = {}
    
    # Get upstream sequences based on the assembly
    #path = os.path.abspath(os.path.join(directory, 'rRNAs/final_16s_seqs', f'{assembly}.fna.txt'))
    
    with open(path_16s, 'r') as file:
        rrna_16s[assembly] = file.read().strip()
    print("Got 16s rRNA sequences")

    # Get 23s sequences as an additional feature for the model
    rrna_23s = {}
    
    # Get upstream sequences based on the assembly
    #path = os.path.abspath(os.path.join(directory, 'rRNAs/final_23s_seqs', f'{assembly}.fna.txt'))
    
    with open(path_23s, 'r') as file:
        rrna_23s[assembly] = file.read().strip()
    print("Got 23s rRNA sequences")

    # Prepare the data
    X = []
    count = 0
    kmer_val = kmer_size
    X_dna_chars = []
    headers_DNA = []
    X_protein_chars = []
    headers_PROT = []
    X_rrna_chars = []
    headers_RRNA = []
    sep_features = []

    # Add each upstream sequence of fasta as a separate feature
    for header, sequence in upstream[assembly].items():
        # Replace sequences made entirely of 'X' with empty strings
        if sequence == 'X' * len(sequence):
            sep_features.append('X' * kmer_val)
        else:
            # Clean up sequences and handle unwanted characters
            cleaned_sequence = sequence.split('*')[0].replace('-', '').replace('K', 'G').replace('M', 'A').replace('R', 'G').replace('Y', 'C').replace('S', 'G')
            sep_features.append(cleaned_sequence)
        
        # Collect headers only once if count is 0
        if count == 0:
            headers_DNA.append(header)
        
    X_dna_chars.append(sep_features)
    #print(f"Collected headers: {headers_DNA}")

    sep_features = []
    # Add each protein sequence of fasta as a separate feature
    for header, sequence in proteins[assembly].items():
        if sequence == 'X' * len(sequence):
            sep_features.append('X'*kmer_val)
        else:
            sep_features.append(sequence.split('*')[0].replace('-', '').replace('B', 'D'))
        if count==0:
            headers_PROT.append(header)
    X_protein_chars.append(sep_features)
    # print(count)

    
    sep_features = []
    
    # Get the 5s rRNA sequence
    rna_5s = rrna_5s[assembly]
    # Get the 16s rRNA sequence
    rna_16s = rrna_16s[assembly]
    # Get the 23s rRNA sequence
    rna_23s = rrna_23s[assembly]
    if count==0:
        headers_RRNA.append(['rna_5s', 'rna_16s', 'rna_23s'])
    
    sep_features.append(rna_5s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
    sep_features.append(rna_16s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
    sep_features.append(rna_23s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
    
    X_rrna_chars.append(sep_features)
    
    # Return the prepared data
    return X, X_dna_chars, headers_DNA, X_protein_chars, headers_PROT, X_rrna_chars, headers_RRNA 


def process_kmers(kmer_size, assembly, headers_DNA, X_dna_chars, headers_PROT, X_protein_chars, headers_RRNA, X_rrna_chars, directory):
    """
    Processes k-mers for DNA, protein, and rRNA sequences and returns k-mer counts.

    Parameters:
    kmer_size (int): Size of the k-mers.
    dna_vocab_path (str): File path to the k-mer vocabulary for DNA.
    protein_vocab_path (str): File path to the k-mer vocabulary for protein.
    rrna_vocab_path (str): File path to the k-mer vocabulary for rRNA.
    headers_DNA (list): List of headers for DNA sequences.
    X_dna_chars (list): List of DNA sequences corresponding to the headers.
    headers_PROT (list): List of headers for protein sequences.
    X_protein_chars (list): List of protein sequences corresponding to the headers.
    headers_RRNA (list): List of headers for rRNA sequences.
    X_rrna_chars (list): List of rRNA sequences corresponding to the headers.

    Returns:
    tuple: Tuple containing lists of k-mer counts for DNA, protein, and rRNA.
    """

    # Load the k-mer vocabulary for DNA, proteins, and rRNA
    kmer_vocab_dna = np.load('../Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flatteneddna_stats_.npy', allow_pickle=True)
    kmer_vocab_protein = np.load('../Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedprotein_stats_.npy', allow_pickle=True)
    kmer_vocab_rrna = np.load('../Machine_Learning/build_dataset/data/kmer_'+str(kmer_size)+'_flattenedrrna_stats_.npy', allow_pickle=True)


    # Convert the loaded data into DataFrames
    dna_kmer_df = pd.DataFrame(kmer_vocab_dna, columns=['Index', 'Kmer', 'DNA'])
    protein_kmer_df = pd.DataFrame(kmer_vocab_protein, columns=['Index', 'Kmer', 'Protein'])
    rrna_kmer_df = pd.DataFrame(kmer_vocab_rrna, columns=['Index', 'Kmer', 'rRNA'])
    
    # Clean the 'rRNA' column to remove parentheses and commas
    rrna_kmer_df['rRNA'] = rrna_kmer_df['rRNA'].apply(lambda x: x[0])

    # Group k-mers by their associated DNA, protein, and rRNA
    dna_kmer_dict = dna_kmer_df.groupby('DNA')['Kmer'].apply(list).to_dict()
    protein_kmer_dict = protein_kmer_df.groupby('Protein')['Kmer'].apply(list).to_dict()
    rrna_kmer_dict = rrna_kmer_df.groupby('rRNA')['Kmer'].apply(list).to_dict()

    # Initialize lists to hold k-mer counts for DNA, proteins, and rRNA
    all_kmer_counts_dna = []
    all_kmer_counts_protein = []
    all_kmer_counts_rrna = []

    # Process k-mers for DNA sequences
    for dna, kmer_vocab in dna_kmer_dict.items():
        dna_key = dna.split('_')[0]
        matching_index = None
        for i, header in enumerate(headers_DNA):
            header_key = header.split('_')[0]
            if header_key == dna_key:
                matching_index = i
                break

        if matching_index is not None:
            sequence = X_dna_chars[0][matching_index]
            vectorizer = CountVectorizer(analyzer='char', ngram_range=(len(kmer_vocab[0]), len(kmer_vocab[0])), vocabulary=kmer_vocab)
            kmer_matrix = vectorizer.transform([sequence])
            kmer_counts = kmer_matrix.toarray().flatten().tolist()
            all_kmer_counts_dna.extend(kmer_counts)
        else:
            print(f"DNA {dna} not found in headers_DNA. Filling with 0s.")
            all_kmer_counts_dna.extend([0] * len(kmer_vocab))

    # Process k-mers for protein sequences
    for protein, kmer_vocab in protein_kmer_dict.items():
        protein_key = protein.split('_')[0]
        matching_index = None
        for i, header in enumerate(headers_PROT):
            header_key = header.split('_')[0]
            if header_key == protein_key:
                matching_index = i
                break

        if matching_index is not None:
            sequence = X_protein_chars[0][matching_index]
            vectorizer = CountVectorizer(analyzer='char', ngram_range=(len(kmer_vocab[0]), len(kmer_vocab[0])), vocabulary=kmer_vocab)
            kmer_matrix = vectorizer.transform([sequence])
            kmer_counts = kmer_matrix.toarray().flatten().tolist()
            all_kmer_counts_protein.extend(kmer_counts)
        else:
            print(f"Protein {protein} not found in headers_PROT. Filling with 0s.")
            all_kmer_counts_protein.extend([0] * len(kmer_vocab))

    # Process k-mers for rRNA sequences
    for rrna, kmer_vocab in rrna_kmer_dict.items():
        matching_index = None
        for i, header in enumerate(headers_RRNA[0]):
            if header == rrna:
                matching_index = i
                break

        if matching_index is not None:
            sequence = X_rrna_chars[0][matching_index]
            vectorizer = CountVectorizer(analyzer='char', ngram_range=(len(kmer_vocab[0]), len(kmer_vocab[0])), vocabulary=kmer_vocab)
            kmer_matrix = vectorizer.transform([sequence])
            kmer_counts = kmer_matrix.toarray().flatten().tolist()
            all_kmer_counts_rrna.extend(kmer_counts)
        else:
            print(f"rRNA {rrna} not found in headers_RRNA. Filling with 0s.")
            all_kmer_counts_rrna.extend([0] * len(kmer_vocab))

    #np.save(directory+'kmer_' + str(kmer_size) + '_flattened_' + "rrna_" + assembly + '_.npy', all_kmer_counts_rrna)
    #np.save(directory+'kmer_' + str(kmer_size) + '_flattened_' + "dna_" + assembly + '_.npy', all_kmer_counts_dna)
    #np.save(directory+'kmer_' + str(kmer_size) + '_flattened_' + "protein_" + assembly + '_.npy', all_kmer_counts_protein)
	
	# Define the single file name with the compressed extension
    filename = os.path.basename(assembly)
    file_path = directory + 'kmer_' + str(kmer_size) + '_flattened_' + filename + '_.npz'

    # Save all three numpy arrays into one compressed file
    np.savez_compressed(file_path,
                    all_kmer_counts_rrna=all_kmer_counts_rrna,
                    all_kmer_counts_dna=all_kmer_counts_dna,
                    all_kmer_counts_protein=all_kmer_counts_protein)

# Return the k-mer counts
    return all_kmer_counts_dna, all_kmer_counts_protein, all_kmer_counts_rrna


# Set up argument parser
parser = argparse.ArgumentParser(description='Process k-mers for DNA, proteins, and rRNA sequences.')
parser.add_argument('--assembly', type=str, required=True, help='Name of assembly.')
parser.add_argument('--kmer_size', type=int, required=True, help='Size of the k-mers.')
parser.add_argument('--path_dna', type=str, required=True, help='File path to DNA.')
parser.add_argument('--path_protein', type=str, required=True, help='File path to  proteins.')
parser.add_argument('--path_5s', type=str, required=True, help='File path to 5s.')
parser.add_argument('--path_16s', type=str, required=True, help='File path to  16s.')
parser.add_argument('--path_23s', type=str, required=True, help='File path to 23s.')
parser.add_argument('--directory', type=str, required=True, help='File path to directory of kmers saved')

# Parse arguments
args = parser.parse_args()

X, X_dna_chars, headers_DNA, X_protein_chars, headers_PROT, X_rrna_chars, headers_RRNA = build_dataset(args.kmer_size, args.assembly, args.path_dna, args.path_protein, args.path_5s, args.path_16s, args.path_23s)

process_kmers(args.kmer_size, args.assembly, headers_DNA, X_dna_chars, headers_PROT, X_protein_chars, headers_RRNA, X_rrna_chars, args.directory)
