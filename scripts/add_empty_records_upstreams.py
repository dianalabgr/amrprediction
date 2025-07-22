import re


# Get Fasta headers
headers = []
current_header = None

# Read Cd-hit fasta file
with open(snakemake.input[1], 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('>'):
            current_header = line.split(' ')[0][1:]
            headers.append(current_header)


# Read upstream/protein sequences
sequences = {}
current_header = None

with open(snakemake.input[0], 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('>'):
            current_header = line.split(' ')[0][1:]
            sequences[current_header] = ''
        else:
            sequences[current_header] += line


# Find with regex if headers are in the sequences.keys() and append this record if does not exist.

for header in headers:
    pattern = re.compile(r"^{}".format(re.escape(header)))
    
    if not any(pattern.search(key) for key in sequences.keys()):
        sequences[header] = 'X'*300

# Sort sequences by headers
sequences = dict(sorted(sequences.items()))

# Write the new fasta file
with open(snakemake.output[0], 'w') as file:
    for header, sequence in sequences.items():
        file.write(f'>{header}\n{sequence}\n')