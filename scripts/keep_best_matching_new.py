import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Keep best mathcing sequences')
parser.add_argument('--file', help='Path to BLAST file (tabular format)')
parser.add_argument('--fasta', help='Path to fasta file')
parser.add_argument('--output', help='Output file to write the best matching sequences to')
parser.add_argument('--db', help='Database file (cd-hit output fasta)')

args = parser.parse_args()
blast_file = args.file
fasta_file = args.fasta
output_file = args.output
db_file = args.db

#read original fasta files 
def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        entry_name = None
        sequence_length = 0

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if entry_name is not None:
                    sequences[entry_name] = sequence_length
                entry_name = line[1:]
                entry_name = entry_name.split()[0]
                sequence_length = 0
            else:
                sequence_length += len(line)

        if entry_name is not None:
            sequences[entry_name] = sequence_length

    return sequences

#Read the length of the different contigs
sequences_dict = read_fasta(fasta_file)


# Function to get length from dictionary
def get_length(name):
    return sequences_dict.get(name, 0)  # Default to 0 if name not found


# Create a dataframe from the BLAST file
df = pd.read_csv(blast_file, sep='\t', header=None)

df.columns = ['query_id', 'subject_id', 'pct_identity', 'aln_length', 'n_of_mismatches', 'gap_openings', 
              'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score', 'q_seq', 's_seq', 'slen']


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
df['pct_coverage'] = df.apply(lambda row: row['aln_length'] / row['slen'], axis=1) * 100

# Get the contig length for each contig
df['contig_length'] = df["query_id"].apply(get_length)

# Filter s_start > 20
df = df[((df['s_start'] <= 20) | (df['q_start']==1) | (df['q_start']==df["contig_length"]) |
		(df['q_end']==1) | (df['q_end'] == df["contig_length"]))]


# filter the dataframe to keep total_pct_identity > 70 and pct_identity > 70
df = df[(df['pct_coverage'] > 70) & (df['pct_identity'] > 70) | 
		((df['q_start']==1) & (df['pct_identity'] > 70)) | 
		((df['q_end'] == df["contig_length"]) & (df['pct_identity'] > 70)) |
		((df['q_end']==1) & (df['pct_identity'] > 70)) | 
		((df['q_start'] == df["contig_length"]) & (df['pct_identity'] > 70))]


# Keep only the proteins with different coordinates of the same contig 
# Assume filtered_df2 is your original DataFrame with the necessary columns
df1 = pd.DataFrame()  # Empty DataFrame to store the filtered rows

# Iterate through the grouped DataFrame by `query_id`
for index, row in df.iterrows():
    inner_break = False
    node=row.iloc[0]
    min_value = min(row.iloc[7], row.iloc[6])
    max_value = max(row.iloc[7], row.iloc[6])
    #print(node)
    for index_2, row_2 in df[df["query_id"] == node].iterrows():
        #print("second",row_2.iloc[0],row_2.iloc[1])
        if row.equals(row_2):
            continue
        else: 
            min_value_2 = min(row_2.iloc[7], row_2.iloc[6])
            max_value_2 = max(row_2.iloc[7], row_2.iloc[6])
    #TO check if there is overlapping 
            if (min_value>=max_value_2):
                continue
            elif (max_value<=min_value_2):
                continue
            #This if else is refering if only a proportion of the row2 is inside row1, here I will set a threshold of overlapping of
    #20% of overlapping to be a problem 
            elif (min_value>=min_value_2) & (min_value<=max_value_2) & (max_value>=max_value_2):
                ratio=float((max_value_2-min_value)/(max_value_2-min_value_2))
                #print("ratio1")
                #print(ratio)
                if (ratio>=0.5) : 
                    if (row.iloc[10]<row_2.iloc[10]) and (row.iloc[11]>row_2.iloc[11] ):
                        #    print("Yeah")
                            continue
                    elif (row.iloc[10]<=row_2.iloc[10]) and (row.iloc[11]>row_2.iloc[11]):
                    #    print("Yeah")
                        continue
                    elif (row.iloc[10]<row_2.iloc[10]) and (row.iloc[11]>=row_2.iloc[11]):
                    #    print("Yeah")
                        continue
                    elif (row.iloc[10]==row_2.iloc[10]) and (row.iloc[11]==row_2.iloc[11]) \
                    and (row.loc["pct_coverage"]==row_2.loc["pct_coverage"]):
                     #   print("Yeah")
                        continue
                    elif (row.iloc[10]==row_2.iloc[10]) and (row.iloc[11]==row_2.iloc[11]) \
                    and (row.loc["pct_coverage"]>row_2.loc["pct_coverage"]):
                      #  print("No1111")
                        inner_break = True
                        break  # Break out of the inner loop
                    else:
                     #   print("No22222")
                        inner_break = True
                        break  # Break out of the inner loop
                else:
                    continue
            elif (min_value<=min_value_2)&(max_value<=max_value_2)&(max_value>=min_value_2):
                ratio=float((max_value-min_value_2)/(max_value_2-min_value_2))
               # print(ratio)
                if (ratio>=0.5) : 
                    if (row.iloc[10]<row_2.iloc[10]) and (row.iloc[11]>row_2.iloc[11]):
                        #    print("yes")
                            continue
                    elif (row.iloc[10]<=row_2.iloc[10]) and (row.iloc[11]>row_2.iloc[11]):
                    #    print("yes")
                        continue
                    elif (row.iloc[10]<row_2.iloc[10]) and (row.iloc[11]>=row_2.iloc[11]):
                   #     print("yes")
                        continue
                    elif (row.iloc[10]==row_2.iloc[10]) and (row.iloc[11]==row_2.iloc[11])\
                    and (row.loc["pct_coverage"]==row_2.loc["pct_coverage"]):
                    #    print("Yeah")
                        continue
                    elif (row.iloc[10]==row_2.iloc[10]) and (row.iloc[11]==row_2.iloc[11]) \
                    and (row.loc["pct_coverage"]>row_2.loc["pct_coverage"]):
                   #     print("No2")
                        inner_break = True
                        break  # Break out of the inner loop
                    else:
                    #    print("No2")
                        inner_break = True
                        break  # Break out of the inner loop
                else: 
                    continue
            elif (min_value>=min_value_2) & (min_value<=max_value_2) & (max_value<=max_value_2):
                ratio=float((max_value-min_value)/(max_value_2-min_value_2))
                if ((ratio>=0.5) and (row.loc["pct_coverage"]>=80)):
                    if (row.iloc[10]<=row_2.iloc[10]) and (row.iloc[11]>=row_2.iloc[11]):
                     #   print("Yeah")
                        continue
                    
                    else:
                     #   print("No")
                        inner_break = True
                        break  # Break out of the inner loop
                else:
                    continue
            else:
                continue
    
    if (inner_break==False):
        # Break out of the outer loop if a break condition was met in the inner loop
        #df1 = df1.append(row, ignore_index=True)
        row_df = pd.DataFrame([row])

        # Concatenate df1 with the new row DataFrame
        df1 = pd.concat([df1, row_df], ignore_index=True)
    
df2=df1.drop_duplicates(subset=["bit_score","q_seq"])                    
result_df = df2.reindex(columns=df2.columns.tolist())

# Based on subject_id, remove duplicates and keep the best matching sequence
result_df = result_df.drop_duplicates(subset='subject_id', keep='first')


result_df.to_csv(output_file, sep='\t', index=False)

