suppressPackageStartupMessages({
  library(ggplot2)
  library(shiny)
  library(plotly)
  library(dplyr)
  library(DT)
})


# ----------------------------
# 0. Load Dataframes
# ----------------------------
ndaro_final_antobiograms <- fread(input = "../data_acquisition/Antibiograms/ndaro_antibiograms.csv", header = T)
bv_brc_filtered <- fread(input = "../data_acquisition/Antibiograms/bv_brc_antibiograms.csv", header = T)
narms_data_filtered_new <- fread(input = "../data_acquisition/Antibiograms/cdc_narms_antibiograms.csv", header = T)

strain_metadata_filtered <- fread(input = "../data_acquisition/Assemblies/strains_metadata.csv", header = T)


# ----------------------------
# 1. Append dataframes
# ----------------------------
# Add BV BRC data frame from `bv_brc_analysis.R`
final_table <- data.frame(
  "antibiotic" = bv_brc_filtered$Antibiotic,
  "sign" = bv_brc_filtered$`Measurement Sign`,
  "value" = bv_brc_filtered$`Measurement Value`,
  "phenotype" = bv_brc_filtered$`Resistant Phenotype`,
  "measurement" = bv_brc_filtered$`Measurement Unit`,
  "organism" = bv_brc_filtered$`Genome Name`,
  "platform" = bv_brc_filtered$`Laboratory Typing Platform`,
  "standard" = bv_brc_filtered$`Testing Standard`,
  "database" = "BV BRC",
  "id" = bv_brc_filtered$Data
)

# Add NDARO data frame from `ndaro_split_antibiotics.R`
final_table <- rbind(final_table, data.frame(
  "antibiotic" = ndaro_final_antobiograms$Antibiotic,
  "sign" = ndaro_final_antobiograms$`Measurement sign`,
  "value" = ndaro_final_antobiograms$Measurement,
  "phenotype" = ndaro_final_antobiograms$`Resistance phenotype`,
  "measurement" = ndaro_final_antobiograms$`Measurement units`,
  "organism" = ndaro_final_antobiograms$Organisms,
  "platform" = ndaro_final_antobiograms$`Laboratory typing platform`,
  "standard" = ndaro_final_antobiograms$`Testing standard`,
  "database" = "NDARO",
  "id" = ndaro_final_antobiograms$Data
))

narms_antibiotics <- c("AMI" = "amikacin",
                       "AMP" = "ampicillin",
                       "ATM" = "aztreonam",
                       "AUG" = "amoxicillin-clavulanic acid",
                       "AXO" = "ceftriaxone",
                       "AZM" = "azithromycin",
                       "CAZ" = "ceftazidime",
                       "CCV" = "ceftazidime-clavulanic acid",
                       "CEP" = "cephalothin",
                       "CEQ" = "cefquinome",
                       "CHL" = "chloramphenicol",
                       "CIP" = "ciprofloxacin",
                       "CLI" = "clindamycin",
                       "COL" = "colistin",
                       "COT" = "trimethoprim-sulfamethoxazole",
                       "CTC" = "cefotaxime-clavulanic acid",
                       "CTX" = "cefotaxime",
                       "ERY" = "erythromycin",
                       "FEP" = "cefepime",
                       "FFN" = "florfenicol",
                       "FIS" = "sulfisoxazole",
                       "FOS" = "fosfomycin",
                       "FOX" = "cefoxitin",
                       "GEN" = "gentamicin",
                       "HYG" = "hygromycin",
                       "IMI" = "imipenem",
                       "KAN" = "kanamycin",
                       "LIN" = "lincomycin",
                       "LZD" = "linezolid",
                       "MER" = "meropenem",
                       "NAL" = "naladixic acid",
                       "PTZ" = "piperacillin-tazobactam",
                       "RIF" = "rifampin",
                       "SMX" = "sulfamethoxazole",
                       "STR" = "streptomycin",
                       "TEL" = "telithromycin",
                       "TET" = "tetracycline",
                       "TIO" = "ceftiofur",
                       "TMP" = "trimethoprim",
                       "TOB" = "tobramycin")

# Add CDC NARMS data frame from `narms_split_antibiotics.R`
for (i in seq(1, nrow(narms_data_filtered_new))) {
  antibiotic <- substring(narms_data_filtered_new$Antibiotic[i], 1, 3)
  final_table <- rbind(final_table, 
                              data.frame("organism" = paste(narms_data_filtered_new$Genus[i], narms_data_filtered_new$Species[i]),
                                         "antibiotic" = narms_antibiotics[[antibiotic]],
                                         "sign" = narms_data_filtered_new$Equiv[i],
                                         "value" = narms_data_filtered_new$Result[i],
                                         "phenotype" = narms_data_filtered_new$Phenotype[i],
                                         "measurement" = NA,
                                         "platform" = NA,
                                         "standard" = NA,
                                         "database" = "CDC NARMS",
                                         "id" = narms_data_filtered_new$Data[i]))
}


# ----------------------------
# 2. Plot antibiotic counts
# ----------------------------
ggplot(data.frame(sort(table(final_table$antibiotic)))) +
  aes(x = Var1, y = Freq) +
  geom_col(fill = "#228B22") +
  geom_text(aes(label = Freq), hjust = -0.7,  size = 2.6) +
  coord_flip() +
  labs(x = "Antibiotic") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# ----------------------------
# 3. Remove low counts antibiotics
# ----------------------------
cutoff <- 400
antibiotic_counts <- data.frame(sort(table(final_table$antibiotic))) |>
  filter(Freq > cutoff)

final_table <- final_table |>
  filter(antibiotic %in% antibiotic_counts$Var1)


# ----------------------------
# 4. Data Preprocess
# ----------------------------

# Remove discrepancies due to capital letters
final_table$antibiotic <- tolower(final_table$antibiotic)

# Change amipicillin-sulbactam to ampicillin-sulbactam
final_table$antibiotic[which(final_table$antibiotic == "amipicillin-sulbactam")] <- "ampicillin-sulbactam"

# Preprocess Typing standard values
final_table$standard[which(final_table$standard == "clsi")] <- "CLSI"
final_table$standard[which(final_table$standard == "missing")] <- NA
final_table$standard[which(final_table$standard == "")] <- NA

# Transform signs
final_table$sign <- ifelse(final_table$sign == "", "==", final_table$sign)
final_table$sign[final_table$sign == "="] <- "=="
final_table$sign[final_table$sign == "³ "] <- "=="


# ----------------------------
# 5. Change phenotype to S, R, I & N
# ----------------------------
final_table$phenotype <- sapply(final_table$phenotype, function(x){
  ifelse(grepl("resistant", tolower(x)), "R", 
         ifelse(grepl("nonsusceptible", tolower(x)), "NS",
                ifelse(grepl("susceptible", tolower(x)), "S", 
                       ifelse(grepl("intermediate", tolower(x)), "I",
                              ifelse(grepl("not defined", tolower(x)), "N", x)))))
})


# ----------------------------
# 6. Keep assemblies only from strains metadata
# ----------------------------
strain_metadata_filtered <- read.csv("../data_acquisition/Assemblies/strains_metadata.csv")

final_table <- final_table |>
  filter(id %in% unique(na.omit(c(strain_metadata_filtered$sra,
                                  strain_metadata_filtered$biosample,
                                  strain_metadata_filtered$assembly))))

write.csv(final_table, "../data_acquisition/Antibiograms/total_antibiograms_argis.csv")



# ----------------------------
# 7. Shiny App. RUN from here
# ----------------------------

# Group values by ranges
ranges <- c(0, 0.1, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 2048, Inf)
labels <- c("0-0.1", "0.1-0.5", "0.5-1", "1-2", "2-4","4-8", "8-16", "16-32", "32-64", "64-128", "128-256", "256-2048", "2048+")

final_table$mic_group <- cut(final_table$value, breaks = ranges, labels = labels, include.lowest = TRUE)

# Make values as numeric
final_table$value <- as.numeric(sapply(strsplit(final_table$value, "/"), "[[", 1))

# Run shiny App
shinyApp(
  fluidPage(
  titlePanel("Antibiogram MIC Distribution"),
  sidebarLayout(
    sidebarPanel(
      selectInput("antibiotics", "Select Antibiotics:",
                  choices = c("Select All", names(sort(table(total_antibiograms$antibiotic), decreasing = T))),
                  multiple = TRUE, selected = "Select All"),
      actionButton("plotBtn", "Generate Plot"),
      width = 3
    ),
    mainPanel(
      plotlyOutput("barPlot"),
      dataTableOutput("summaryTable"),
      width = 9
    )
  )
),
function(input, output) {
  sum_of_rows <- reactive({
    df_sum <- aggregate(value ~ antibiotic, data = total_antibiograms, FUN = length)
    colnames(df_sum)[2] <- "Sum of Rows"
    df_sum
  })
  
  
  output$barPlot <- renderPlotly({
    req(input$plotBtn)
    
    if ("Select All" %in% input$antibiotics) {
      selected_data <- total_antibiograms
    } else {
      selected_data <- total_antibiograms[total_antibiograms$antibiotic %in% input$antibiotics, ]
    }
    
    ggplotly(
      ggplot(selected_data, aes(x = interaction(sign, " (", mic_group, "]", sep = ""))) +
      geom_bar(fill = "skyblue") +
      geom_text(stat = "count", aes(label = after_stat(count)), vjust = 1, size = 3) +
      labs(x = "MIC Signs and Ranges", y = "MIC Counts", title = "Antibiogram MIC Barplot") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
            plot.title = element_text(hjust = 0.5))
    )
  })
  
  output$summaryTable <- renderDataTable({
    req(input$plotBtn)

    if ("Select All" %in% input$antibiotics) {
      selected_data <- total_antibiograms
    } else {
      selected_data <-
        total_antibiograms[total_antibiograms$antibiotic %in% input$antibiotics,]
    }

    summary_table <-
      aggregate(value ~ antibiotic, data = selected_data, FUN = length)
    colnames(summary_table)[2] <- "Sum of Rows"
    summary_table
  })
})


# ----------------------------
# 8. Additional Notes
# ----------------------------

# Split antibiotic combination rows to a second dataframe
final_table_double <- final_table |>
  filter(antibiotic %in% c("trimethoprim-sulfamethoxazole", 
                           "piperacillin-tazobactam", 
                           "amoxicillin-clavulanic acid",
                           "ceftazidime-avibactam",
                           "ampicillin-sulbactam",
                           "ceftolozane-tazobactam",
                           "cefotaxime-clavulanic acid",
                           "ceftazidime-clavulanic acid"))

  
nrow(final_table[which(final_table$sign == "<" & final_table$value == 16),])

# Normalize values, based on sign (for further analysis)
final_table$value[which(final_table$sign == "<")] <- 
  final_table$value[which(final_table$sign == "<")] / 2

