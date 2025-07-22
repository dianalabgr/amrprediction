from Bio import SeqIO
import pandas as pd

gene_names = []
input = "../annotation/concatenated_fasta.fasta"
output = "../annotation/concatenated_fasta_appened.fasta"

card = "../annotation/amr_databases/aro_index_v3.8.9.tsv"
amrfinder = "../annotation/amr_databases/ReferenceGeneCatalog_v3.12.txt"


with open(input, "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_names.append(record.name)
        
card = pd.read_csv(card, sep = "\t")
amrfinder = pd.read_csv(amrfinder, sep = "\t")

def concatenate_gene_names(element):
    updated_header = {}
    try:
        query = amrfinder.loc[(amrfinder['genbank_protein_accession'] == element) | (amrfinder['refseq_protein_accession'] == element), 'gene_family']
        if len(query) > 0:
            updated_header[element] = element + "_" + query.iloc[0]
        else:
            if element == "sp|P24734.3|AMPR_PSEAE":
                element = 'P24734.3'
            updated_header[element] = element + "_" + card.loc[card['Protein Accession'] == element, 'CARD Short Name'].iloc[0]
    except IndexError:
        print('Index error', element)
    except Exception as e:
        print(e)
    
    return updated_header[element]


updated_headers = {element: concatenate_gene_names(element) for element in gene_names}

updated_records = []
with open(input, "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        record.id = updated_headers[record.id]
        record.description = record.description
        updated_records.append(record)


with open(output, "w") as output_file:
    SeqIO.write(updated_records, output_file, "fasta")
