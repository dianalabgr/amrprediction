library(dplyr)
library(stringr)
library(readr)
library(data.table)

# Load the BV BRC AMR data
bv_brc <- fread("R/initial_datasets/BVBRC_genome_amr.csv", header = TRUE, sep = ',', colClasses = c("Genome ID" = "character"))

# We want to isolate only the below 7 ESKAPEE pathogens from the dataset
strains <- c("Klebsiella", "Enterococcus", "Acinetobacter", "Pseudomonas", "Staphylococcus", "Enterobacter", "Escherichia")


bv_brc_filtered <- data.frame()

# Keep only the strains we are interested in
for (i in strains) {
  bv_brc_filtered <- rbind(bv_brc_filtered, bv_brc |>
    filter(if_any(where(is.character), ~str_detect(.x, i)))
  )
}

# Remove the columns that are not needed
bv_brc_filtered <- bv_brc_filtered[,2:14] |>
  select(-c(Measurement, `Laboratory Typing Method`, `Laboratory Typing Method Version`, Vendor))


# Change cells ending with '/'
adjustments <- grep("/$", bv_brc_filtered$`Measurement Value`)
for (i in adjustments) {
  bv_brc_filtered$`Measurement Value`[i] <- as.numeric(sub("/$", "", bv_brc_filtered$`Measurement Value`[i]))
}

# bv_brc_filtered$`Measurement Value` <- lapply(bv_brc_filtered$`Measurement Value`, function(x) eval(parse(text = x)))
bv_brc_filtered$`Measurement Value` <- as.character(bv_brc_filtered$`Measurement Value`)

bv_brc_filtered$Antibiotic <- gsub("/", "-", bv_brc_filtered$Antibiotic)


# --------------------------------
# Combine Genome IDs with NCBI IDs
# --------------------------------
# Description: We will combine the Genome IDs with the NCBI IDs to get the final assembly files
# We will use the NCBI API to get the SRA Run accession numbers for the Genome IDs, through the Entrez API
# If the assembly accession number is not available, we will use the SRA accession number to get the assembly accession number

genome_data <- read_csv("R/initial_datasets/BVBRC_genome.csv", col_types = list(`Genome ID`= col_character()))

genome_data <- genome_data |>
  select(c(`Genome ID`, `SRA Accession`, `Assembly Accession`))

genome_data$`Genome ID` <- format(as.numeric(genome_data$`Genome ID`), nsmall = 6)
bv_brc_filtered$`Genome ID` <- format(as.numeric(bv_brc_filtered$`Genome ID`), nsmall = 6)

bv_brc_filtered_new <- bv_brc_filtered |>
  distinct(`Genome ID`)

assemblies <- vector("list", length=nrow(bv_brc_filtered_new))
for (i in 1:nrow(bv_brc_filtered_new)) {
  temp <- genome_data[match(bv_brc_filtered_new$`Genome ID`[i], genome_data$`Genome ID`),]
  
  if (is.na(temp$`Assembly Accession`)) {
    if (!is.na(temp$`SRA Accession`) & substring(temp$`SRA Accession`, 1, 3) == "ERS") {
      assemblies[[i]] <- system(paste('efetch -db sra -id', temp$`SRA Accession`,
                                      '-format docsum | xtract -pattern Runs -element Run@acc'), intern = T)
    } else {
      assemblies[[i]] <- temp$`SRA Accession`
    }
  } else {
    assemblies[[i]] <- temp$`Assembly Accession`
  }
}

bv_brc_filtered_new$Data <- do.call(rbind, assemblies)

bv_brc_filtered <- bv_brc_filtered |>
  inner_join(bv_brc_filtered_new, by = "Genome ID")


# Remove rows where SRA files were unpaired, so we do not have final assembly
bv_brc_assemblies <- list.files("../data_acquisition/Assemblies/bv_brc/", full.names = TRUE)

bv_brc_filtered$isValid <- sapply(bv_brc_filtered$Data, function(x) {
  ifelse(length(grep(x, bv_brc_assemblies)), TRUE, NA)
}, USE.NAMES = F)

# Get unique valid assemblies
length(unique(bv_brc_filtered$Data[which(!is.na(bv_brc_filtered$isValid))]))

# Remove non valid assemblies
bv_brc_filtered <- bv_brc_filtered |>
  filter(!is.na(isValid))


# Some columns has more than 1 IDs separated with commas. We keep the first of them
bv_brc_filtered$Data <- sapply(bv_brc_filtered$Data, function(x) {
  strsplit(x, ",")[[1]][1]
}, USE.NAMES = F)


bv_brc_filtered <- bv_brc_filtered |>
  filter(!is.na(Data))

# Save the file and copy it to antibiograms folder
fwrite(x = bv_brc_filtered, file = "R/final_datasets/bv_brc_antibiograms.csv", sep = ",")

file.copy("R/final_datasets/bv_brc_antibiograms.csv", "../data_acquisition/Antibiograms/")
