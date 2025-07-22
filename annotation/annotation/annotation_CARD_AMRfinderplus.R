# Packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("dplyr")
})


# CARD data
card = read.csv(
  file = '../annotation/annotation/amr_databases/aro_index_v3.8.9.tsv',
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

# AMRfinderplus data
amrfinderplus = read.csv(
  file = '../annotation/annotation/amr_databases/ReferenceGeneCatalog_v3.12.txt',
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

# Counting and verifying that I have 86 empty cells in CARD data in the column Protein.Accession
empty_string_count <- sum(card$Protein.Accession == "")
print(empty_string_count)

# No empty strings for nucleotides
empty_string_nuc <- sum(card$DNA.Accession == "")
print(empty_string_nuc)

# Counting and verifying that I have 363 empty cells in AMRfinderplus data in the column genbank_protein_accession
empty_string_count <- sum(amrfinderplus$genbank_protein_accession == "")
print(empty_string_count)

# No empty strings for nucleotides
empty_string_nuc <- sum(amrfinderplus$genbank_nucleotide_accession == "")
print(empty_string_nuc)


# convert empty cells to NAs
card$Protein.Accession[card$Protein.Accession == ""] <- NA
amrfinderplus$genbank_protein_accession[amrfinderplus$genbank_protein_accession == ""] <- NA
amrfinderplus$refseq_protein_accession[amrfinderplus$refseq_protein_accession == ""] <- NA


# Counting NAs
# CARD
protein_NA_card = sum(is.na(card$Protein.Accession))
nucleotide_NA_card = sum(is.na(card$DNA.Accession))

# AMRfinderplus
protein_NA_amr = sum(is.na(amrfinderplus$genbank_protein_accession))
nucleotide_NA_amr = sum(is.na(amrfinderplus$genbank_nucleotide_accession))
rotein_NA_amr_refseq = sum(is.na(amrfinderplus$refseq_protein_accession))
protein_NA_card
nucleotide_NA_card
protein_NA_amr
nucleotide_NA_amr


# Divide AMRfinderplus data to core (antimicrobial resistance genes)
# and plus(other than AMR, e.g virulence genes, biocide resistance genes etc)
# Note: for 3 entries scope was "non-reported"

amr_core <- subset(amrfinderplus, scope == 'core')

amr_plus <- subset(amrfinderplus, scope == 'plus')

head(amr_core)
dim(amr_core)
head(amr_plus)
dim(amr_plus)


#make separate dfs for protein and nucleotide accession numbers, for each database. Exclude NAs

# AMR
# protein
amr_core_protein <- na.omit(data.frame(genbank_protein_accession = ifelse(is.na(amr_core$genbank_protein_accession), 
                                                                          amr_core$refseq_protein_accession, 
                                                                          amr_core$genbank_protein_accession)))


# amr_core_protein = subset(amr_core, !is.na(genbank_protein_accession), select = 'genbank_protein_accession')

amr_plus_protein <- na.omit(data.frame(genbank_protein_accession = ifelse(is.na(amr_plus$genbank_protein_accession), 
                                                                          amr_plus$refseq_protein_accession, 
                                                                          amr_plus$genbank_protein_accession)))
# amr_plus_protein = subset(amr_plus, !is.na(genbank_protein_accession), select = 'genbank_protein_accession')
# nucleotide
amr_core_nucleotide = subset(amr_core,!is.na(genbank_nucleotide_accession), select = 'genbank_nucleotide_accession')
amr_plus_nucleotide = subset(amr_plus,!is.na(genbank_nucleotide_accession), select = 'genbank_nucleotide_accession')


# CARD
# protein
card_protein = subset(card, !is.na(Protein.Accession), select = 'Protein.Accession')
# nucleotide
card_nucleotide = subset(card, !is.na(DNA.Accession), select = 'DNA.Accession')

# Protein
print("Proteins")
dim(card_protein)[1]
dim(amr_core_protein)[1]
dim(amr_plus_protein)[1]


# nucleotide
print("Nucleotides")
dim(card_nucleotide)[1]
dim(amr_core_nucleotide)[1]
dim(amr_plus_nucleotide)[1]


# The number of overlapping accessions between CARD and AMRfinderplus core
common_protein_accessions <- card_protein$Protein.Accession %in% amr_core_protein$genbank_protein_accession
common_nucleotide_accessions <- card_nucleotide$DNA.Accession %in% amr_core_nucleotide$genbank_nucleotide_accession

print("Between CARD and AMRfinderplus core")
cat("Number of common protein accessions:",
    sum(common_protein_accessions),
    "\n")

cat("Number of common nucleotide accessions:",
    sum(common_nucleotide_accessions),
    "\n")


# The number of overlapping accessions between CARD and AMRfinderplus plus
common_protein_accessions <-
  card_protein$Protein.Accession %in% amr_plus_protein$genbank_protein_accession
common_nucleotide_accessions <-
  card_nucleotide$DNA.Accession %in% amr_plus_nucleotide$genbank_nucleotide_accession
print("Between CARD and AMRfinderplus plus")
cat("Number of common protein accessions:",
    sum(common_protein_accessions),
    "\n")
cat("Number of common nucleotide accessions:",
    sum(common_nucleotide_accessions),
    "\n")


# CARD
length(unique(card_protein$Protein.Accession))
length(unique(card_nucleotide$DNA.Accession))

# AMR core
length(unique(amr_core_protein$genbank_protein_accession))
length(unique(amr_core_nucleotide$genbank_nucleotide_accession))

# AMR plus
length(unique(amr_plus_protein$genbank_protein_accession))
length(unique(amr_plus_nucleotide$genbank_nucleotide_accession))


# Concatenating keeping only unique values
# Core proteins
print("core proteins")
protein_all_core_unique <-
  unique(c(
    card_protein$Protein.Accession,
    amr_core_protein$genbank_protein_accession
  ))
length(protein_all_core_unique)

# core nucleotides
print("core nucleotides")
nucleotide_all_core_unique <-
  unique(
    c(
      card_nucleotide$DNA.Accession,
      amr_core_nucleotide$genbank_nucleotide_accession
    )
  )
length(nucleotide_all_core_unique)

# plus proteins
print("plus proteins")
protein_plus_unique <-
  unique(amr_plus_protein$genbank_protein_accession)
length(protein_plus_unique)

# intersection between protein_plus_unique and protein_all_core_unique
intersect_protein <-
  intersect(protein_all_core_unique, protein_plus_unique)
length(intersect_protein)

# remove the accession numbers that can also be found in core
protein_plus_unique <-
  setdiff(protein_plus_unique, intersect_protein)
length(protein_plus_unique)

# plus nucleotides
print("plus nucleotides")
nucleotide_plus_unique <-
  unique(amr_plus_nucleotide$genbank_nucleotide_accession)
length(nucleotide_plus_unique)

# intersection between protein_plus_unique and protein_all_core_unique
intersect_nucleotide <-
  intersect(nucleotide_all_core_unique, nucleotide_plus_unique)
length(intersect_nucleotide)

# remove the accession numbers that can also be found in core
nucleotide_plus_unique <-
  setdiff(nucleotide_plus_unique, intersect_nucleotide)
length(nucleotide_plus_unique)



# save to text files where each line is an accession number

# core proteins
protein_core <- sapply(protein_all_core_unique, as.character)
writeLines(protein_core, "../annotation/annotation/protein_accessions_core.txt")

# # core nucleotides
# nucleotide_core <- sapply(nucleotide_all_core_unique, as.character)
# writeLines(nucleotide_core, "../annotation/nucleotide_accessions_core.txt")

# plus proteins
protein_plus <- sapply(protein_plus_unique, as.character)
writeLines(protein_plus, "../annotation/annotation/protein_accessions_plus.txt")

# # plus nucleotides
# nucleotide_plus <- sapply(nucleotide_plus_unique, as.character)
# writeLines(nucleotide_plus, "../annotation/nucleotide_accessions_plus.txt")


# For genes that have not Refseq protein id and Genbank protein id
non_translated_genes <- amrfinderplus$refseq_nucleotide_accession[which(is.na(amrfinderplus$refseq_protein_accession) & 
                                                                  is.na(amrfinderplus$genbank_protein_accession))]

sum(is.na(non_translated_genes))
writeLines(unique(non_translated_genes), "../annotation/annotation/non_translated_genes.txt")
