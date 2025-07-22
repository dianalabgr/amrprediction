#BiocManager::install("Biostrings")
library(Biostrings)

args <- commandArgs(TRUE)

file   <- args[1]
coords <- args[2]
output_dna  <- args[3]
output_prot <- args[4]

assembly <- Biostrings::readDNAStringSet(file)

coords <- data.table::fread(coords, sep = "\t", header = T)

names(assembly) <- sapply(names(assembly), function(x) {
  strsplit(x, " ")[[1]][1]
}, USE.NAMES = F)


upstream_length <- 300


# Get Upstream sequences
up_sequences <- data.frame()
protein <- data.frame()

for (i in 1:nrow(coords)) {
  # Get protein aligned sequence
  protein <- rbind(protein, data.frame(seqs = coords$q_seq[i], names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
  
  if (coords$q_start[i] < coords$q_end[i]) {
    # Forward strand
    if ((coords$q_start[i] - upstream_length) >=1 ) {
      up_sequences <- rbind(up_sequences, data.frame(seqs = substr(toString(assembly[[coords$query_id[i]]]), 
                                                          coords$q_start[i] - upstream_length, coords$q_start[i] - 1), 
                                                   names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
    }
    else if (coords$q_start[i] > 1 ) {
      up_sequences <- rbind(up_sequences, data.frame(seqs = substr(toString(assembly[[coords$query_id[i]]]), 
                                                                   1, coords$q_start[i] - 1), 
                                                     names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
    }
    else {
      up_sequences <- rbind(up_sequences, data.frame(seqs = paste(rep("N", 300), collapse = ""), 
                                                     names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
    }
  }
                       
  else {
    # Reverse strand
    if ((coords$q_end[i] - upstream_length) >=1 ) {
    up_sequences <- rbind(up_sequences, data.frame(seqs = toString(Biostrings::reverseComplement(Biostrings::DNAString(substr(assembly[[coords$query_id[i]]], 
                                                          coords$q_end[i] - upstream_length, coords$q_end[i] - 1)))),
                                                   names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
    }
    else if (coords$q_end[i] > 1 ) {
      up_sequences <- rbind(up_sequences, data.frame(seqs = substr(toString(assembly[[coords$query_id[i]]]), 
                                                                   1, coords$q_end[i] - 1), 
                                                     names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
    }
    else {
      up_sequences <- rbind(up_sequences, data.frame(seqs = paste(rep("N", 300), collapse = ""), 
                                                     names = paste(coords$subject_id[i], coords$query_id[i], sep = "_")))
    }
  }
}

up_sequences$seqs <- sapply(up_sequences$seqs, function(seq) {
  seq_length <- nchar(seq)
  if (seq_length < 300) {
    padding <- strrep("N", 300 - seq_length)
    paste0(padding, seq)
  } else {
    seq
  }
})

# Write upstream sequences to file with FASTA format
dna_set <- Biostrings::DNAStringSet(up_sequences$seqs)
names(dna_set) <- up_sequences$names


# Write aligned protein sequence in FASTA format
protein_set <- Biostrings::AAStringSet(protein$seqs)
names(protein_set) <- protein$names

Biostrings::writeXStringSet(dna_set, output_dna)
Biostrings::writeXStringSet(protein_set, output_prot)

