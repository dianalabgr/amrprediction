suppressPackageStartupMessages({
  library(readxl)
})


######################################
# Data Loading and Preprocessing
######################################
narms_data <- read_xlsx("R/initial_datasets/narms_IsolateData.xlsx")

narms_antibiogams <- cbind(Genus = narms_data$Genus,
                           Species = narms_data$Species, narms_data[, 17:length(colnames(narms_data))])


######################################
# Data transformation
######################################
# Calculate the total number of rows for the final dataframe
total_rows <- length(seq(3, length(colnames(narms_antibiogams)), by = 4)) * nrow(narms_antibiogams)

# Preallocate the dataframe
narms_antibiogams_new <- data.frame(
  Genus = character(total_rows),
  Species = character(total_rows),
  Antibiotic = character(total_rows),
  Equiv = character(total_rows),
  Result = character(total_rows),
  Phenotype = character(total_rows),
  Data = character(total_rows),
  stringsAsFactors = FALSE
)

# Initialize row index
row_index <- 1

# Function to get the appropriate database IDs
get_database_ids <- function(ncbi_accession, wgs_id) {
  id <- ""
  if (grepl("^SAMN", ncbi_accession) &
      grepl("^20", wgs_id, "^20")) {
    id <- ncbi_accession
  } else if (grepl("^SAMN", ncbi_accession) &
             grepl("^PNUSA", wgs_id)) {
    id <- wgs_id
  } else if (is.na(ncbi_accession) &
             grepl("^PNUSA", wgs_id)) {
    id <- wgs_id
  } else if (!grepl("^SAMN", ncbi_accession) &
             grepl("^PNUSA", wgs_id)) {
    id <- wgs_id
  } else {
    id <- NA
    # print(c(ncbi_accession, wgs_id))
  }
  
  return(id)
}

# Populate the pre-allocated dataframe
for (i in seq(3, length(colnames(narms_antibiogams)), by = 4)) {
  for (j in seq(1, nrow(narms_antibiogams))) {
    antibiotic <- substring(colnames(narms_antibiogams)[i], 1, 3)

    narms_antibiogams_new[row_index, ] <- data.frame(
      Genus = narms_antibiogams$Genus[j],
      Species = narms_antibiogams$Species[j],
      Antibiotic = antibiotic,
      Equiv = narms_antibiogams[j, i],
      Result = narms_antibiogams[j, i + 1],
      Phenotype = narms_antibiogams[j, i + 2],
      Data = get_database_ids(narms_data$`NCBI Accession Number`[j], narms_data$`WGS ID`[j]),
      stringsAsFactors = FALSE
    )

    # Move to the next row
    row_index <- row_index + 1
  }
}


# filter rows with NAs in Phenotype and Data
narms_antibiogams_new <- narms_antibiogams_new[!is.na(narms_antibiogams_new$Phenotype), ]
narms_antibiogams_new <- narms_antibiogams_new[!is.na(narms_antibiogams_new$Data), ]

# Save the final data frame
write.csv(x = narms_antibiogams_new, file = "../data_acquisition/Antibiograms/cdc_narms_antibiograms.csv")
write.csv(x = narms_antibiogams_new, file = "R/final_datasets/cdc_narms_data.csv")


######################################
# Split the data into antibiotics
# One file per antibiotic
######################################
antibiotics <- unique(narms_antibiogams_new$Antibiotic)

# Split the data into antibiotics
for (antibiotic in antibiotics) {
  antibiotic_data <- narms_antibiogams_new[narms_antibiogams_new$Antibiotic == antibiotic, ]
  write.csv(x = antibiotic_data, file = paste0("../data_acquisition/Antibiograms/cdc_narms/", antibiotic, ".csv"))
}
