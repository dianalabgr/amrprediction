library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggmap)

# Load the filtered data after the quality control
strain_metadata_filtered <- read.csv("/mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/strains_metadata.csv")

# ---------------------
# Plot 1: Location Distribution
# ---------------------
install.packages(c("maps", "ggmap"))

data <- data.frame("location" = strain_metadata_filtered$location, "genus" = strain_metadata_filtered$genus)

city_dict <- list(
  "" = "",
  "Region 10" = "Region 10",
  "Region 7" = "Region 7",
  "Region 5" = "Region 5",
  "Region 9" = "Region 9",
  "Region 1" = "Region 1",
  "Region 4" = "Region 4",
  "Region 8" = "Region 8",
  "Region 6" = "Region 6",
  "Region 2" = "Region 2",
  "Region 3" = "Region 3",
  "USA: Boston, MA" = "Boston",
  "USA: Boston" = "Boston",
  "USA: California" = "California",
  "USA" = "USA",
  "USA:Boston" = "Boston",
  "USA:TN" = "Tennessee",
  "USA:FL" = "Florida",
  "USA:OH" = "Ohio",
  "USA:CO" = "Colorado",
  "USA:ID" = "Idaho",
  "USA:WI" = "Wisconsin",
  "USA:SC" = "South Carolina",
  "USA:CA" = "California",
  "USA:PA" = "Pennsylvania",
  "USA:TX" = "Texas",
  "USA:KS" = "Kansas",
  "USA:NE" = "Nebraska",
  "USA:NY" = "New York",
  "USA: New York" = "New York",
  "USA:IA" = "Iowa",
  "USA:NJ" = "New Jersey",
  "Brazil: Santos" = "Santos",
  "Taiwan: Tainan" = "Tainan",
  "United Arab Emirates: Tawam Hospital, Al Ain" = "Al Ain",
  "Thailand" = "Thailand",
  "Thailand: Bangkok" = "Bangkok",
  "USA:MS" = "Mississippi",
  "USA:GA" = "Georgia",
  "USA:WA" = "Washington",
  "USA:IN" = "Indiana",
  "USA:ND" = "North Dakota",
  "USA:MI" = "Michigan",
  "USA:OK" = "Oklahoma",
  "USA:MN" = "Minnesota",
  "USA:NC" = "North Carolina",
  "Sweden: Gothenburg" = "Gothenburg",
  "USA:KY" = "Kentucky",
  "USA:OR" = "Oregon",
  "USA:SD" = "South Dakota",
  "USA:MD" = "Maryland",
  "USA:MT" = "Montana",
  "USA:VA" = "Virginia",
  "USA:AR" = "Arkansas",
  "USA:MO" = "Missouri",
  "USA: Rochester, NY" = "Rochester",
  "USA:IL" = "Illinois",
  "Canada" = "Canada",
  "USA:VT" = "Vermont",
  "USA:AL" = "Alabama",
  "USA:MA" = "Massachusetts",
  "USA:Cincinnati" = "Cincinnati",
  "India" = "India",
  "USA: OH" = "Ohio",
  "USA: WA" = "Washington",
  "Canada: Quebec" = "Quebec",
  "Netherlands" = "Netherlands",
  "USA: MO" = "Missouri",
  "USA: VA" = "Virginia",
  "USA:CT" = "Connecticut",
  "USA: IN" = "Indiana",
  "USA:WV" = "West Virginia",
  "USA: KY" = "Kentucky",
  "USA: PA" = "Pennsylvania",
  "USA: NY" = "New York",
  "USA: MA" = "Massachusetts",
  "USA: IL" = "Illinois",
  "Canada:New Brunswick" = "New Brunswick",
  "Canada:Prince Edward Island" = "Prince Edward Island",
  "Canada:Nova Scotia" = "Nova Scotia",
  "Canada:Newfoundland" = "Newfoundland",
  "USA: GA" = "Georgia",
  "USA: ID" = "Idaho",
  "Canada: 0" = "Canada",
  "USA: SC" = "South Carolina",
  "USA: Scotsdale, AZ" = "Scotsdale",
  "USA: Temple, TX" = "Temple",
  "USA: Minnesota" = "Minnesota",
  "USA: Unknown" = "",
  "OUTPATIENT" = "",
  "USA: MI" = "Michigan",
  "USA: Ohio" = "Ohio",
  "Brazil: Sao Paulo" = "Sao Paulo",
  "USA: District of Columbia" = "District of Columbia",
  "USA:Kansas" = "Kansas",
  "USA:Kentucky" = "Kentucky",
  "Botswana: Kasane" = "Kasane",
  "Portugal: Aveiro" = "Aveiro",
  "Greece" = "Greece",
  "Portugal: Lisbon" = "Lisbon",
  "India: Anna Nagar, Chennai, TN" = "Chennai",
  "USA: Alaska" = "Alaska",
  "USA:Washington, DC" = "Washington, DC",
  "USA: Washington, DC" = "Washington, DC",
  "USA:Alaska" = "Alaska",
  "USA: Maryland, Bethesda" = "Bethesda",
  "Turkey: Bursa" = "Bursa",
  "Germany: Cologne" = "Cologne",
  "Germany: Bochum" = "Bochum",
  "Germany: Leverkusen" = "Leverkusen",
  "Belgium: Verviers" = "Verviers",
  "Germany: Berlin" = "Berlin",
  "Germany: Braunschweig" = "Braunschweig",
  "Germany: somewhere in BW" = "Baden-Württemberg",
  "Germany: somewhere in SN" = "Saxony",
  "Germany: somewhere in TH" = "Thuringia",
  "Germany: Dresden" = "Dresden",
  "Germany: Cottbus" = "Cottbus",
  "Germany: Schwerin" = "Schwerin",
  "Germany: Munich" = "Munich",
  "Germany: Mainz" = "Mainz",
  "Germany: NRW" = "North Rhine-Westphalia",
  "Germany: somewhere in NRW" = "North Rhine-Westphalia",
  "Australia: Melbourne" = "Melbourne",
  "Pakistan" = "Pakistan",
  "Germany" = "Germany",
  "Ukraine" = "Ukraine",
  "Afghanistan" = "Afghanistan",
  "Peru" = "Peru",
  "Honduras" = "Honduras",
  "Israel:Tel-Aviv" = "Tel-Aviv",
  "Thailand: northern region" = "Thailand",
  "Thailand: northeastern region" = "Thailand",
  "USA: Cleveland, OH" = "Cleveland",
  "USA: Irvine, CA" = "Irvine",
  "USA:Maryland,Bethesa" = "Bethesa",
  "Germany:Landstuhl" = "Landstuhl",
  "USA: Connecticut" = "Connecticut",
  "USA: Indiana" = "Indiana",
  "USA: Oklahoma" = "Oklahoma",
  "USA: Illinois" = "Illinois",
  "USA: Alabama" = "Alabama",
  "USA: Washington" = "Washington",
  "USA: Virginia" = "Virginia",
  "USA: Pennsylvania" = "Pennsylvania",
  "USA: Texas" = "Texas",
  "Taiwan: Taipei" = "Taipei",
  "USA: Houston" = "Houston",
  "USA: MN" = "Minnesota",
  "Russia: Kazan" = "Kazan",
  "South Africa: Gauteng" = "Gauteng",
  "Israel" = "Israel",
  "Brazil: Parana, Curitiba" = "Curitiba",
  "China" = "China",
  "Australia" = "Australia",
  "Australia: Victoria" = "Victoria",
  "USA: Nevada" = "Nevada",
  "Ukraine: Vinnitsa" = "Vinnitsa",
  "South Korea" = "South Korea",
  "Kenya: Nairobi" = "Nairobi",
  "Philippines" = "Philippines",
  "USA: North Carolina" = "North Carolina",
  "Italy: Sienna" = "Sienna",
  "USA: Arizona" = "Arizona",
  "Germany: Landstuhl" = "Landstuhl",
  "Jordan" = "Jordan",
  "Israel: Jerusalem" = "Jerusalem",
  "USA: Maryland" = "Maryland",
  "USA: New Jersey" = "New Jersey",
  "USA: Conneticut" = "Connecticut",
  "Guam: Hagatna" = "Hagatna",
  "Peru: Iquitos" = "Iquitos",
  "Turkey: Ankara" = "Ankara",
  "USA: Hawaii" = "Hawaii",
  "South Africa: Pretoria" = "Pretoria",
  "USA: Ann Arbor, Michigan" = "Ann Arbor",
  "Portugal" = "Portugal",
  "Brazil: Curitiba" = "Curitiba",
  "Italy: Milan" = "Milan",
  "Canada: Hamilton, Ontario" = "Hamilton",
  "Ecuador: Quito" = "Quito",
  "Mexico" = "Mexico",
  "China: Shanghai" = "Shanghai",
  "Saudi Arabia: Riyadh" = "Riyadh",
  "Iran: Sanandaj" = "Sanandaj",
  "USA: Rochester, MN" = "Rochester",
  "USA: Cleveland, Ohio" = "Cleveland",
  "USA: CT" = "Connecticut",
  "USA:PR" = "Puerto Rico",
  "USA: MT" = "Montana",
  "USA: OR" = "Oregon",
  "USA: VT" = "Vermont",
  "USA: Cambridge, MA" = "Cambridge",
  "USA:Nebraska" = "Nebraska",
  "USA:California" = "California",
  "Russia: Moscow" = "Moscow",
  "USA: Detroit" = "Detroit",
  "Australia: Queensland" = "Queensland",
  "Thailand: Pathum Thani" = "Pathum Thani",
  "USA: Detroit, MI" = "Detroit",
  "USA: Madison, WI" = "Madison"
)


# based on city_dict, take each strain_metadata_filtered$location and create a new list with the respective city
cleaned_locations <- lapply(data$location, function(x) {
  if (is.null(city_dict[[x]])) {
    return(x)
  } else {
    return(city_dict[[x]])
  }
})

data$location <- unlist(cleaned_locations)

# Group the data by location and genus and count the number of occurrences
data <- data |>
  group_by(location, genus) |>
  summarise(count = n()) |>
  ungroup()

colnames(data) <- c("Location", "Genus", "Frequency")

# Remove empty rows and rows starting with "Region"
data <- data |>
  filter(Location != "" & !grepl("Region", Location))

# Make some arrangements to the location names, to find them on the map
data$Location[data$Location == "Bethesda"] <- "Bethesda, Maryland"
data$Location[data$Location == "Sienna"] <- "Siena, Italy"
data$Location[data$Location == "Santos"] <- "General Santos City, Philippines"
data$Location[data$Location == "Scotsdale"] <- "Scottsdale, Arizona"
data$Location[data$Location == "Victoria"] <- "Victoria, Australia"
data$Location[data$Location == "Hamilton"] <- "Hamilton, Ontario"
data$Location[data$Location == "Milan"] <- "milan, metropolitan city of milan, italy"
data$Location[data$Location == "Sao Paulo"] <- "são paulo, state of são paulo, brazil"
data$Location[data$Location == "Tainan"] <- "tainan, east district, tainan city, taiwan"
data$Location[data$Location == "Washington"] <- "washington, dc, usa"


register_google(key = "API_KEY")

# Get coordinates for each unique location
location_coords <- geocode(data$Location, output = "latlon", source = "google")
location_data <- cbind(data, location_coords)

# Load world map data
library(maps)
world_map <- map_data("world")
# World map
world_map <- world_map[world_map$region != "Antarctica", ]

# Plot the data
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray93", color = "white") +
  geom_point(data = location_data, aes(x = lon, y = lat, size = Frequency, color = Genus), alpha = 0.7) +
  theme_minimal() +
  labs(title = "Isolation Locations of Strains",
    x = NULL, y = NULL,
    subtitle = "The size of the points represents the frequency of the genus in the respective location") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Genus", nrow = 1, override.aes = list(size = 5))) +
  theme(legend.title = element_text(family = "Times New Roman", face = "bold", size = 12),
        legend.text = element_text(family = "Times New Roman", size = 11),
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.title = element_text(family = "Times New Roman", size = 20, face = "bold"),
        plot.subtitle = element_text(family = "Times New Roman", size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())


library(echarts4r)
library(dplyr)
strain_metadata_filtered <- read.csv("/mnt/raid1b/kdan_data/Paper/Data_Acquisition/Assemblies/strains_metadata.csv")

strain_metadata_filtered |>
  count(database) |>
  e_charts(database) |>
  e_pie(n, radius = c("40%", "80%")) |>
  e_title("Distribution of Used Databases",
          "The pie chart shows the distribution of the databases used.\n",
          textStyle = list(color = "black", fontFamily = "Times New Roman", fontSize = 24),
          subtextStyle = list(color = "black", fontFamily = "Times New Roman", fontSize = 18, lineHeight = 0)) |>
  e_theme("walden") |>
  e_legend(show = FALSE) |>
  e_color(color = c("#6be6c1", "#c4ebad", "#3fb1e3")) |>
  e_labels(normal = list(
    show = TRUE,
    position = "outside",
    formatter = htmlwidgets::JS("function(params) {
      return params.name + '\\n' + params.percent + '%';
    }"),
    textStyle = list(fontFamily = "Times New Roman", fontSize = 17, color = "black")
  ))


strain_metadata_filtered |>
  count(genus) |>
  e_charts(genus) |>
  e_pie(n, radius = c("40%", "80%")) |>
  e_title(
    "Genus Distribution of Strains",
    "The pie chart shows the distribution of the genera in the strains dataset.\n",
    textStyle = list(color = "black", fontFamily = "Times New Roman", fontSize = 24),
    subtextStyle = list(color = "black", fontFamily = "Times New Roman", fontSize = 18, lineHeight = 0)) |>
  e_theme("walden") |>
  e_legend(show = FALSE) |>
  e_color(color = c("#3fb1e3", "#6be6c1", "#626c91", "#a0a7e6", "#c4ebad", "#ff9f7f", "#96dee8")) |>
  e_labels(normal = list(
    show = TRUE,
    position = "outside",
    formatter = htmlwidgets::JS("function(params) {
      return params.name + '\\n' + params.percent + '%' + '\\n' + params.value;
    }"),
    textStyle = list(fontFamily = "Times New Roman", fontSize = 17, color = "black")
  ))



# ---------------------
# Plot 2: Antibiotic Distribution
# ---------------------
# I combined the QC'ed metadata of all antibiotics into one file and I was trying to gain insights from them
combined_metadata <- NULL

for (folder in list.files("/mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/")) {
  if (folder != "combined_metadata.csv") {
    file_name <- tolower(folder)
    
    if (folder == "amoxicillin-clavulanicAcid") {
      file_name <- "amoxicillin-clavulanic"
    }

    metadata <- read.csv(paste0("/mnt/raid1b/kdan_data/Paper/Data_Acquisition/Antibiograms/QC_antibiotics/used_antibiotics/", folder, "/", gsub("-", "_", file_name), "_QCed_metadata.csv"))

    if (exists("combined_metadata")) {
      combined_metadata <- rbind(combined_metadata, metadata)
    } else {
      combined_metadata <- metadata
    }
  }
}

# Remove rows where antibiotic == vancomycin from combined_metadata, due to poor quality of results
combined_metadata <- combined_metadata[!(combined_metadata$antibiotic == "vancomycin"), ]

ggplot(combined_metadata) +
  aes(x = antibiotic, fill = phenotype) +
  geom_bar() +
  labs(title = "Antibiotic Distribution of Strains",
       subtitle = "The bar chart shows the distribution of the antibiotic resistance phenotypes in the strains dataset.\n",
       x = "Antibiotic",
       y = "Count") +
  scale_fill_viridis_d(option = "viridis", direction = -1, 
                    breaks = c("S", "I", "R", "N"),
                    labels = c("Susceptible", "Intermediate", "Resistant", "Not Stated"),
                    name = "Phenotype") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45L, hjust = 1L),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 12),
        legend.text = element_text(family = "Times New Roman", size = 10),
        plot.title = element_text(family = "Times New Roman", size = 20, face = "bold"),
        plot.subtitle = element_text(family = "Times New Roman", size = 14),
        axis.text = element_text(family = "Times New Roman", size = 10, color = "gray24"),
        axis.title = element_text(family = "Times New Roman", size = 12))


# ---------------------
# Plot 3: Venn Diagram
# ---------------------
library(tidyverse)

circle <- function(center_x, center_y, radius, label = "", resolution = 1000) {
  x <- seq(-radius, radius, length.out = resolution)
  y_upper <- sqrt(radius^2 - x^2)
  y_lower <- -sqrt(radius^2 - x^2)

  tibble(
    x = x + center_x,
    ymin = y_lower + center_y,
    ymax = y_upper + center_y,
    label = label
  )
}

bind_rows(
  circle(0, 0, 1, "A"),
  circle(1, 0, 1, "B"),
  circle(0.5, -0.9, 1, "C")
) %>%
  ggplot(aes(x = x, ymin = ymin, ymax = ymax, group = label, fill = label, alpha = label)) +
  geom_ribbon(color = "black", linewidth = 1, show.legend = FALSE) +
  annotate("text",
    x = c(-0.5, 0.5, 1.5), y = c(0.2, 0.4, 0.2),
    label = c("1,081", "32", "35"), fontface = "bold"
  ) +
  annotate("text",
    x = c(-0.1, 0.5, 1.1), y = c(-0.6, -0.3, -0.6),
    label = c("782", "211", "14"), fontface = "bold"
  ) +
  annotate("text",
    x = 0.5, y = -1.3,
    label = "399", fontface = "bold"
  ) +
  annotate("text",
    x = c(-0.5, 0.5, 1.5), y = c(1.2, -2, 1.2),
    label = c(
      "AMRfinderPlus\n8,157",
      "CARD\n5,078",
      "ResFinder\n3,150"
    ),
    fontface = "bold", lineheight = 0.8, vjust = 1
  ) +
  coord_equal() +
  scale_fill_manual(
    breaks = c("A", "B", "C"),
    values = c("#96dee8", "#ff9f7f", "#3fb1e3")
  ) +
  scale_alpha_discrete(
    breaks = c("A", "B", "C"),
    range = c(1, 0.25)
  ) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

c("#3fb1e3", "#6be6c1", "#626c91", "#a0a7e6", "#c4ebad", "#ff9f7f", "#96dee8")
