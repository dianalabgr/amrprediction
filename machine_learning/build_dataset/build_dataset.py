import os
import pickle


def build_dataset(kmer_size):
    files = os.listdir("../../Data_Acquisition/Assemblies/filtered_assemblies/")
    files = [os.path.basename(file).split('.fna')[0] for file in files]

    # Remove newline characters from each line
    # files = [line.strip() for line in files]
    
    print(len(files))
    
    # Get upstream sequences as an additional feature for the model
    upstream = {}

    for assembly in files:
        # Get upstream sequences based on the assembly
        path = os.path.abspath(f'../../annotation/annotated_assemblies/upstream_seqs_final/{assembly}.fna')
        
        with open(path, 'r') as file:
            # read the fasta file and store a dictionary with the assembly as the key and the sequence as the value
            fasta = {}
            for line in file:
                if line.startswith('>'):
                    header = line.strip()
                    fasta[header] = ''
                else:
                    fasta[header] += line.strip()

            upstream[assembly] = fasta
    print("Got upstream sequences")

    # Get protein sequences as an additional feature for the model
    proteins = {}

    for assembly in files:
        # Get upstream sequences based on the assembly
        path = os.path.abspath(f'../../annotation/annotated_assemblies/aligned_proteins_final/{assembly}.fna')
        
        with open(path, 'r') as file:
            # read the fasta file and store a dictionary with the assembly as the key and the sequence as the value
            fasta = {}
            for line in file:
                if line.startswith('>'):
                    header = line.strip()
                    fasta[header] = ''
                else:
                    fasta[header] += line.strip()

            proteins[assembly] = fasta
    print("Got protein sequences")

    # Get 5s rRNA sequences as an additional feature for the model
    rrna_5s = {}

    for assembly in files:
        # Get upstream sequences based on the assembly
        path = os.path.abspath(f'../../annotation/annotated_assemblies/rRNAs/final_5s_seqs/{assembly}.fna.txt')
        
        with open(path, 'r') as file:
            rrna_5s[assembly] = file.read().strip()
    print("Got 5s rRNA sequences")

    # Get 16s rRNA sequences as an additional feature for the model
    rrna_16s = {}

    for assembly in files:
        # Get upstream sequences based on the assembly
        path = os.path.abspath(f'../../annotation/annotated_assemblies/rRNAs/final_16s_seqs/{assembly}.fna.txt')
        
        with open(path, 'r') as file:
            rrna_16s[assembly] = file.read().strip()
    print("Got 16s rRNA sequences")


    # Get 23s rRNA sequences as an additional feature for the model
    rrna_23s = {}

    for assembly in files:
        # Get upstream sequences based on the assembly
        path = os.path.abspath(f'../../annotation/annotated_assemblies/rRNAs/final_23s_seqs/{assembly}.fna.txt')
        
        with open(path, 'r') as file:
            rrna_23s[assembly] = file.read().strip()
    print("Got 23s rRNA sequences")

    print("data loaded")
    
    
    # Prepare the data
    X = []
	
    count=0 
	
    kmer_val = kmer_size
    #kmer_val = 5
    X_dna_chars = []
    headers_DNA=[]
    X_protein_chars = []
    headers_PROT=[]
    X_rrna_chars = []
    headers_RRNA=[]

    for assembly in files:
        if assembly in files:
            
            features = []
            sep_features = []
            # Add each upstream sequence of fasta as a separate feature
            for header, sequence in upstream[assembly].items():
                # we check if sequence is all with X and we change them with ''
                if sequence == 'X' * len(sequence):
                    features.append('X'*kmer_val)
                    sep_features.append('X'*kmer_val)
                else:
                    # remove all charachters after '*'
                    features.append(sequence.split('*')[0].replace('-', '').replace('K', 'G').replace('M', 'A').replace('R', 'G').replace('Y', 'C').replace('S', 'G'))
                    sep_features.append(sequence.split('*')[0].replace('-', '').replace('K', 'G').replace('M', 'A').replace('R', 'G').replace('Y', 'C').replace('S', 'G'))
                if count==0:
                    headers_DNA.append(header)
                
            X_dna_chars.append(sep_features)
            # print(count)
            
            sep_features = []
            # Add each protein sequence of fasta as a separate feature
            for header, sequence in proteins[assembly].items():
                if sequence == 'X' * len(sequence):
                    features.append('X'*kmer_val)
                    sep_features.append('X'*kmer_val)
                else:
                    features.append(sequence.split('*')[0].replace('-', '').replace('B', 'D'))
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
            
            features.append(rna_5s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
            features.append(rna_16s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
            features.append(rna_23s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
            
            sep_features.append(rna_5s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
            sep_features.append(rna_16s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
            sep_features.append(rna_23s.replace('-', '').replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('D', 'G'))
            
            X_rrna_chars.append(sep_features)
            
            # Append the features to the X list
            X.append(features)
            count=1


    pickle.dump(X, open('./data/X_no_kmers.pkl', 'wb'))
    pickle.dump(X_dna_chars, open('./data/X_dna_chars.pkl', 'wb'))
    pickle.dump(headers_DNA, open('./data/headers_dna.pkl', 'wb'))    
    pickle.dump(X_protein_chars, open('./data/X_protein_chars.pkl', 'wb'))
    pickle.dump(headers_PROT, open('./data/headers_protein.pkl', 'wb'))        
    pickle.dump(X_rrna_chars, open('./data/X_rrna_chars.pkl', 'wb'))
    pickle.dump(headers_RRNA, open('./data/headers_rrna.pkl', 'wb')) 
    
    return()


build_dataset(5)
