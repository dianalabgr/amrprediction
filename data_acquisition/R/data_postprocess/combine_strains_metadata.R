library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(cowplot)


# ----------------------------
# 0. Load Dataframes
# ----------------------------
ndaro_final_antobiograms <- read.csv(file = "../data_acquisition/Antibiograms/ndaro_antibiograms.csv")
ndaro_data <- read.csv(file = "R/initial_datasets/NDARO_isolates.csv")

bv_brc_filtered <- read.csv(file = "./Antibiograms/bv_brc_antibiograms.csv")
bv_brc  <- read.csv(file = "R/initial_datasets/BVBRC_genome_amr.csv")

narms_data_filtered_new <- read.csv(file = "../data_acquisition/Antibiograms/cdc_narms_antibiograms.csv")
narms_data <- readxl::read_xlsx(path = "R/initial_datasets/narms_IsolateData.xlsx")


# ----------------------------
# 1. Append dataframes
# ----------------------------

# Add BV BRC strain info
bv_brc_strains <- data.frame(
  "name" = sapply((bv_brc_filtered |> distinct(Data, .keep_all = T))$Genome.ID, function(x) {
    unique(bv_brc$Genome.Name[as.numeric(bv_brc$Genome.ID) == as.numeric(x)][[1]])
  }, USE.NAMES = F),
  "source" = sapply((bv_brc_filtered |> distinct(Data, .keep_all = T))$Genome.ID, function(x) {
    unique(bv_brc$Source[as.numeric(bv_brc$Genome.ID) == as.numeric(x)][[1]])
  }, USE.NAMES = F),
  "host" = NA,
  "year" = sapply((bv_brc_filtered |> distinct(Data, .keep_all = T))$Genome.ID, function(x) {
    unique(bv_brc$Testing.Standard.Year[as.numeric(bv_brc$Genome.ID) == as.numeric(x)][[1]])
  }, USE.NAMES = F),
  "platform" = NA,
  "location" = NA,
  "contigs" = NA,
  "method" = NA,
  "sra" = sapply(unique(bv_brc_filtered$Data), function(x) {
    if (substr(x, 2, 3) == "RR")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "biosample" = sapply(unique(bv_brc_filtered$Data), function(x) {
    if (substr(x, 1, 4) == "SAMN")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "assembly" = sapply(unique(bv_brc_filtered$Data), function(x) {
    if (substr(x, 1, 2) == "GC")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "database" = "BV BRC"
)


# Append NDARO strain info
ndaro_strains <- data.frame(
  "name" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    index <- which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)
    paste(ndaro_data$X.Organism.group[index], ndaro_data$Strain[index])
  }, USE.NAMES = F),
  "source" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    ndaro_data$Isolation.source[which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)]
  }, USE.NAMES = F),
  "host" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    ndaro_data$Host[which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)]
  }, USE.NAMES = F),
  "year" = NA,
  "platform" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    ndaro_data$Platform[which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)]
  }, USE.NAMES = F),
  "location" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    ndaro_data$Location[which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)]
  }, USE.NAMES = F),
  "contigs" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    ndaro_data$Contigs[which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)]
  }, USE.NAMES = F),
  "method" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    ndaro_data$Method[which(ndaro_data$BioSample == x | ndaro_data$Assembly == x)]
  }, USE.NAMES = F),
  "sra" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    if (substr(x, 2, 3) == "RR")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "biosample" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    if (substr(x, 1, 4) == "SAMN")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "assembly" = sapply(unique(ndaro_final_antobiograms$Data), function(x) {
    if (substr(x, 1, 2) == "GC")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "database" = "NDARO"
)


# Remove rows with non-existent assembly ids (25 unique assemblies) due to technical error
list_of_assemblies <- list.files("../data_acquisition/Assemblies/cdc_narms", full.names = TRUE)

narms_data_filtered_new <- narms_data_filtered_new |>
  filter(Data %in% sapply(list_of_assemblies, function(x) {
    strsplit(basename(x), ".fna")[[1]]
  }, USE.NAMES = F))

# Append CDC NARMS strain info
cdc_narms_strains <- data.frame(
  "name" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    index <- which(narms_data$`NCBI Accession Number` == x | narms_data$`WGS ID` == x)
    
    paste(narms_data$Genus[index][[1]],
          narms_data$Species[index][[1]], 
          narms_data$Serotype[index][[1]])
  }, USE.NAMES = F),
    
  "source" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    narms_data$`Specimen Source`[which(narms_data$`NCBI Accession Number` == x | narms_data$`WGS ID` == x)][[1]]
  }, USE.NAMES = F),
  "host" = NA,
  "year" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    narms_data$`Data Year`[which(narms_data$`NCBI Accession Number` == x | narms_data$`WGS ID` == x)][[1]]
  }, USE.NAMES = F),
  "platform" = NA,
  "location" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    narms_data$`Region Name`[which(narms_data$`NCBI Accession Number` == x | narms_data$`WGS ID` == x)][[1]]
  }, USE.NAMES = F),
  "contigs" = NA,
  "method" = NA,
  "sra" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    if (substr(x, 2, 3) == "RR")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "biosample" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    if (substr(x, 1, 4) == "SAMN" | substr(x, 1, 4) == "PNUS")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "assembly" = sapply(unique(narms_data_filtered_new$Data), function(x) {
    if (substr(x, 1, 2) == "GC")
      x[[1]]
    else
      NA
  }, USE.NAMES = F),
  "database" = "CDC NARMS"
)


# ----------------------------
# 2. Remove duplicate ids (between NDARO & BV BRC)
# ----------------------------
# For the same id we remove the rows from BV BRC database, 
# because they contain more NAs in metadata

# Append data frames
strain_metadata <- rbind(bv_brc_strains, cdc_narms_strains, ndaro_strains)

# All duplicates are only assembly ids
length(na.omit(intersect(ndaro_strains$assembly, bv_brc_strains$assembly)))
length(na.omit(intersect(ndaro_strains$biosample, bv_brc_strains$biosample)))
length(na.omit(intersect(ndaro_strains$sra, bv_brc_strains$sra)))


strain_metadata <- strain_metadata |>
  filter(!(assembly %in% na.omit(intersect(ndaro_strains$assembly, bv_brc_strains$assembly)) & database == "BV BRC"))

bv_brc_strains <- bv_brc_strains |>
  filter(!(assembly %in% na.omit(intersect(ndaro_strains$assembly, bv_brc_strains$assembly)) & database == "BV BRC"))


# Export strain data frames
write.csv(ndaro_strains, "../data_acquisition/Assemblies/ndaro_strains.csv")
write.csv(bv_brc_strains, "../data_acquisition/Assemblies/bv_brc_strains.csv")
write.csv(cdc_narms_strains, "../data_acquisition/Assemblies/cdc_narms_strains.csv")

# ----------------------------
# 3. Add N50 values
# ----------------------------

# Add N50 values for NDARO data frame
ids <- apply(data.frame(strain_metadata$sra, strain_metadata$biosample, strain_metadata$assembly), 1, function(row) {
  paste(na.omit(row))
})

strain_metadata$N50 <- sapply(ids, FUN = function(x) {
  ifelse(length(ndaro_data$N50[which(ndaro_data$Assembly == x)]), 
         ndaro_data$N50[which(ndaro_data$Assembly == x)], 
         NA)
})


list_of_assemblies <- c(
  list.files("../data_acquisition/Assemblies/ndaro", full.names = TRUE),
  list.files("../data_acquisition/Assemblies/bv_brc", full.names = TRUE),
  list.files("../data_acquisition/Assemblies/cdc_narms", full.names = TRUE)
)

for (i in 1:nrow(strain_metadata)) {
  tryCatch({
    if (is.na(strain_metadata$N50[i])) {
      strain_metadata$N50[i] <- Biostrings::N50(Biostrings::readDNAStringSet(list_of_assemblies[
        ifelse(!is.na(grep(ids[i], list_of_assemblies)[1]),
               grep(ids[i],
               list_of_assemblies)[1],
        next)])@ranges@width)
    }
  },
  error = function(e) {
    print(e)
  },
  finally = {
    next
  })
}


# ----------------------------
# 4. Add Genus column
# ----------------------------

strain_metadata$genus <- sapply(strain_metadata$name, function(x) {
  strsplit(x = x, split = " ")[[1]][1]
}, USE.NAMES = F)

# Adjust "E.coli" and "Escherichia"
strain_metadata$genus <- ifelse(strain_metadata$genus == "E.coli", "Escherichia", strain_metadata$genus)

unique(strain_metadata$genus)
table(strain_metadata$genus)

# Summary for Klebsiella
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[1])])

# Summary for Enterococcus
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[2])])

# Summary for Acinetobacter
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[3])])

# Summary for Pseudomonas
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[4])])

# Summary for Staphylococcus
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[5])])

# Summary for Enterobacter
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[6])])

# Summary for Escherichia & Shigella!
summary(strain_metadata$N50[which(strain_metadata$genus == unique(strain_metadata$genus)[7])])

source("./R/ggplot_theme_Publication.R")

# Visualize histograms and boxplots
df <- data.frame(table(strain_metadata$genus), N50=1:7)
colnames(df) <- c("genus", "freq", "N50")

bar_unfiltered <- ggplot(strain_metadata) +
  aes(x = genus, y = N50) +
  geom_bar(stat = "summary", fun = "mean", color = "#023743FF", fill = "#476F84FF") +
  scale_y_continuous(labels = number_format(big.mark = ".", decimal.mark = ",")) +
  labs(title = "N50 Mean Values per Genus", subtitle = "Unfiltered Dataset", y = "N50 Values") +
  geom_text(data = df, aes(x = genus, label = paste("counts\n", freq)), color = "white", vjust = -0.3) +
  theme_Publication() +
  theme(
    plot.title = element_text(size = 14L, face = "bold"),
    plot.subtitle = element_text(size = 12L),
    axis.title.y = element_text(size = 14L),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11L, angle = 25, margin = margin(t = 15))
  )


box_unfiltered <- ggplot(strain_metadata) +
  aes(x = genus, y = N50) +
  geom_boxplot(color = "dodgerblue4", lwd = 0.4, outlier.size = 1.4, outlier.alpha = 0.6) +
  scale_y_continuous(labels = number_format(big.mark = ".", decimal.mark = ",")) +
  labs(title = "N50 Boxplots per Genus (Unfiltered)", x = "Genus") +
  theme_gray() +
  theme(
    plot.title = element_text(size = 14L, face = "bold"),
    plot.subtitle = element_text(size = 12L),
    axis.title.y = element_text(size = 14L),
    axis.title.x = element_text(size = 14L)
  )


density_plot <- strain_metadata |>
  # filter(genus == "Staphylococcus") |>
  ggplot() +
    aes(x = N50) +
    geom_density(adjust = 1L, color = "#023743FF", fill = "#476F84FF") +
    scale_y_sqrt() +
    scale_x_continuous(labels = number_format(big.mark = ".", decimal.mark = ","), breaks = seq(0, 7500000, by = 1000000)) +
    labs(title = "N50 Distribution", 
         subtitle = "Filtered Dataset",
         y = "Density", x = "N50 Value") +
    # geom_vline(xintercept = 100000, color = "red") +
    theme_Publication() +
    theme(
      plot.title = element_text(size = 14L, face = "bold"),
      plot.subtitle = element_text(size = 12L),
      axis.title.y = element_text(size = 14L),
      axis.title.x = element_text(size = 13L)
    )


# ----------------------------
# 5. QC assemblies (filter low quality assemblies)
# ----------------------------

# 1607 strains were removed
strain_metadata_filtered <- strain_metadata |>
  filter(!((genus == "Enterobacter" | genus == "Staphylococcus") & N50 < 20000)) |>
  filter(!(N50 < 50000 & genus != "Enterobacter" & genus != "Staphylococcus"))

df <- data.frame(table(strain_metadata_filtered$genus), N50=1:7)
colnames(df) <- c("genus", "freq", "N50")

# ----------------------------
# 6. Plot N50 Distribution
# ----------------------------
bar_filtered <- ggplot(strain_metadata_filtered) +
  aes(x = genus, y = N50) +
  geom_bar(stat = "summary", fun="mean", color = "#023743FF", fill = "#476F84FF") +
  scale_y_continuous(labels = number_format(big.mark = ".", decimal.mark = ",")) +
  labs(title = "N50 Mean Values per Genus", subtitle = "Filtered Dataset") +
  geom_text(data = df, aes(x = genus, label = paste("counts\n", freq)), color = "white", vjust = -0.3) +
  theme_Publication() +
  theme(
    plot.title = element_text(size = 14L, face = "bold"),
    plot.subtitle = element_text(size = 12L),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11L, angle = 25, margin = margin(t = 15))
  )


box_filtered <- ggplot(strain_metadata_filtered) +
  aes(x = genus, y = N50) +
  geom_boxplot(color = "dodgerblue4", lwd = 0.4, outlier.size = 1.4, outlier.alpha = 0.6) +
  scale_y_continuous(labels = number_format(big.mark = ".", decimal.mark = ",")) +
  labs(title = "N50 Boxplots per Genus (Filtered)", x = "Genus") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14L, face = "bold"),
    plot.subtitle = element_text(size = 12L),
    axis.title.y = element_text(size = 14L),
    axis.title.x = element_text(size = 14L)
  )


bottom_row <- plot_grid(density_plot, labels = "C", label_size = 18L, nrow = 1, ncol = 1, align = 'h')
top_row <- plot_grid(bar_unfiltered, bar_filtered, nrow = 1, labels = c("A", "B"), label_size = 18L) +
  draw_label("Bacterial Genus", x = 0.5, y = 0, vjust = -0.5)

plot_grid(top_row, bottom_row, ncol = 1)

grid.arrange(box_unfiltered, box_filtered)


# ----------------------------
# 7. Save file
# ----------------------------

fwrite(strain_metadata_filtered, file = "../data_acquisition/Assemblies/strains_metadata.csv", sep = ",")
