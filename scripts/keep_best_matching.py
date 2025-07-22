import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Keep best mathcing sequences')
parser.add_argument('--file', help='Path to BLAST file (tabular format)')
parser.add_argument('--output', help='Output file to write the best matching sequences to')
parser.add_argument('--db', help='Database file (cd-hit output fasta)')

args = parser.parse_args()
blast_file = args.file
output_file = args.output
db_file = args.db


# Create a dataframe from the BLAST file
df = pd.read_csv(blast_file, sep='\t', header=None)

df.columns = ['query_id', 'subject_id', 'pct_identity', 'aln_length', 'n_of_mismatches', 'gap_openings', 
              'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score', 'q_seq', 's_seq']


# Find length of each protein/dna sequence on db_file
protein_length = {}
with open(db_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            protein_id = line.split(' ')[0][1:].rstrip()
            protein_length[protein_id] = 0
        else:
            protein_length[protein_id] += len(line.strip())


# Create a new column in the dataframe with the total percent identity (aln_length / protein length) protein name is based on second column
df['pct_coverage'] = df.apply(lambda row: (row['aln_length'] - row['gap_openings']) / protein_length[row['subject_id']], axis=1) * 100

# Filter s_start > 5
df = df[df['s_start'] <= 5]


# filter the dataframe to keep total_pct_identity > 70 and pct_identity > 70
df = df[(df['pct_coverage'] > 70) & (df['pct_identity'] > 70)]

# Group dataframe by query_id & subject_id and get the max value of pct_identity
df2 = df.groupby('subject_id').agg({'e_value': 'min', 'pct_identity': 'max'})
result_df = pd.merge(df, df2, on=['subject_id', 'pct_identity', 'e_value'], how='inner')
result_df = result_df.reindex(columns=df.columns.tolist())

# Based on subject_id, remove duplicates and keep the best matching sequence
result_df = result_df.drop_duplicates(subset='subject_id', keep='first')


result_df.to_csv(output_file, sep='\t', index=False)

