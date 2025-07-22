import pandas as pd
import os
import re
import sys

# Redirect output to log file
sys.stdout = open(snakemake.log[0], 'w')


# Initialize path variables
strain_metadata = pd.read_csv("../data_acquisition/Assemblies/strains_metadata.csv")

reference_genomes = {
    'Klebsiella': 'Klebsiella_pneumoniae',
    'Escherichia': 'Escherichia_coli',
    'Enterobacter': 'Enterobacter_cloacae',
    'Pseudomonas': 'Pseudomonas_aeruginosa',
    'Staphylococcus': 'Staphylococcus_aureus',
    'Acinetobacter': 'Acinetobacter_baumannii',
    'Enterococcus': 'Enterococcus_faecium'
}

assembly = os.path.basename(snakemake.input[0])

# get genus from genus column of strain_metadata based on assembly variable
# if assembly is equal to sra, assembly, or biosample, return genus
mask = (strain_metadata["sra"] == assembly.replace(".fna", "")) | \
        (strain_metadata["assembly"] == assembly.replace(".fna", "")) | \
        (strain_metadata["biosample"] == assembly.replace(".fna", ""))
genus = strain_metadata.loc[mask, "genus"].to_string(index=False).strip() if not strain_metadata.loc[mask, "genus"].empty else None

ref_genome = reference_genomes[genus]
dir_name = f"{snakemake.params[0]}{ref_genome}/"

# Get rRNA copies as fasta files
copies_length = []
protein_length = {}
for fasta in os.listdir(dir_name):
    match = re.match(r'^(\d+)\.fasta$', fasta)
    if match:
        copies_length.append(f"{dir_name}{match.group(1)}.fasta")


for copy in copies_length:
    # Find length of each protein/dna sequence on db_file
    with open(copy, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line.split(' ')[0][1:].rstrip()
                protein_length[protein_id] = 0
            else:
                protein_length[protein_id] += len(line.strip())


# Run blastn on the assembly file against all databases in dir_name folder and filter the results
for db in os.listdir(dir_name):
    match = re.match(r'^(db\d+)\.ndb$', db)
    if match:
        os.system(f"../tools/ncbi-blast-2.15.0+/bin/blastn -query {snakemake.input[0]} \
                  -db {dir_name}{match.group(1)} -outfmt \"6 std qseq sseq\" -out {dir_name}{assembly}_{match.group(1)}.out -num_threads {snakemake.threads}")
        
        # Create a dataframe from the BLAST file
        if os.path.getsize(f"{dir_name}{assembly}_{match.group(1)}.out") == 0:
            df = pd.DataFrame()
            df.to_csv(f"{dir_name}{assembly}_{match.group(1)}_filtered.out", sep='\t', header=False, index=False)
        else:
            df = pd.read_csv(f"{dir_name}{assembly}_{match.group(1)}.out", sep='\t', header=None)
        
            df.columns = ['query_id', 'subject_id', 'pct_identity', 'aln_length', 'n_of_mismatches', 'gap_openings', 
                'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score', 'q_seq', 's_seq']
            
            
            # Create a new column in the dataframe with the total percent identity (aln_length / protein length) protein name is based on second column
            df['pct_coverage'] = df.apply(lambda row: (row['aln_length']) / protein_length[row['subject_id']], axis=1) * 100

            # filter the dataframe to keep total_pct_identity > 50 and pct_identity > 70
            df = df[(df['pct_coverage'] > 80) & (df['pct_identity'] > 80)]
            

            df.to_csv(f"{dir_name}{assembly}_{match.group(1)}_filtered.out", sep='\t', header=False, index=False)
            os.system(f"rm {dir_name}{assembly}_{match.group(1)}.out ")

# make regex to match the file name f'.db{number}'
files = []
for db in os.listdir(dir_name):
    if re.match(assembly +  r'_db(\d+)\_filtered.out$', db):
            file = os.path.join(dir_name, db)
            files.append(file)

blastn_data = {}

for file in files:
    if os.path.getsize(file) > 0:
        df = pd.read_csv(file, delimiter = "\t" , header=None)
        blastn_data[df[1].iloc[0]] = df

bitscores = []

for key, values in blastn_data.items():
    max_index = max(values[11])
    bitscores.append(max_index)


data = pd.DataFrame({'contig': [df.iloc[0, 0] for df in blastn_data.values()],
        'gene': blastn_data.keys(),
        'bitscore': bitscores,
        'length': None,
        'sequence': None,
        'start': None,
        'end': None})
data = data.sort_values(by='bitscore', ascending=False)

length = []
for key in data['gene']:
    length.append(protein_length[key])

data['length'] = length

unique_coordinates = set()

if len(blastn_data) == 0:
    for i in range(len(files)):
        # read each fasta file (i+1.fasta) and count the sequence length
        with open(f"{dir_name}{i+1}.fasta", 'r') as f:
            length = 0
            for line in f:
                if not line.startswith('>'):
                    length += len(line.strip())
        
        data.at[i, 'sequence'] = 'X'*length
        data.at[i, 'start'] = 0
        data.at[i, 'end'] = 0
        data.at[i, 'length'] = length
        data.at[i, 'bitscore'] = 0
        data.at[i, 'gene'] = None
        data.at[i, 'contig'] = None
else:
    for index, row in data.iterrows():
        subject_id = row['gene']

        if blastn_data[subject_id].empty:
            data.at[index, 'sequence'] = 'X'*data['length'].iloc[index]
            data.at[index, 'start'] = 0
            data.at[index, 'end'] = 0
        else:
            max_bitscore = blastn_data[subject_id].sort_values(by = 11, ascending=False)
            max_bitscore = max_bitscore.iloc[0]

            if max_bitscore[8] > max_bitscore[9]:
                # reverse complement the sequence
                sequence = max_bitscore[12]
                sequence = sequence[::-1]
                sequence = sequence.translate(str.maketrans('ATCG', 'TAGC'))
            else:
                sequence = max_bitscore[12]
            
            data.at[index, 'sequence'] = sequence
            data.at[index, 'start'] = max_bitscore[6]
            data.at[index, 'end'] = max_bitscore[7]

            unique_coordinates.add((max_bitscore[0], max_bitscore[6], max_bitscore[7]))
            
            # Remove the coordinates from the blastn_data, and also when start and end coords are within the range of the previous coordinates
            for subject_id, df in blastn_data.items():
                blastn_data[subject_id] = df[~df[[0, 6, 7]].apply(
                    lambda x: (x[0], x[6], x[7]) in unique_coordinates or (x[6] >= max_bitscore[6] and x[7] <= max_bitscore[7]), axis=1
                )]
                # add empty dataframes to the dictionary to avoid key errors
                if blastn_data[subject_id].empty:
                    blastn_data[subject_id] = pd.DataFrame(columns=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

print(data)

data = data.sort_values(by='gene', ascending=True)

final_seqs = ''
for index, row in data.iterrows():
    final_seqs = str(final_seqs) + str(row['sequence'])

with open(snakemake.output[0], 'w') as f:
    f.write(final_seqs)
    

    
