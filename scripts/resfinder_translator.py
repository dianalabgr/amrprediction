import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_fasta(input_file, output_file):
    with open(input_file, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    translated_records = []
    for record in records:
        seq = record.seq.translate()
        # Remove stop codons ('*')
        seq = seq.rstrip('*')
        translated_record = SeqRecord(seq, id=record.id, description=record.description)
        translated_records.append(translated_record)

    with open(output_file, "w") as f:
        SeqIO.write(translated_records, f, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    translate_fasta(input_file, output_file)

