# Load libraries
library(UpSetR)
library(dplyr)
library(tidyr)
# Install (if not installed yet)
#install.packages("patchwork")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
# Load patchwork
library(patchwork)

# Your data list
data <- list(
  escherichia = c("amoxicillin-clavulanicAcid", "ampicillin", "ampicillin-sulbactam", "azithromycin", "aztreonam", 
                  "cefazolin", "cefepime", "cefotaxime", "cefoxitin", "cefpodoxime", "ceftazidime", "ceftriaxone", 
                  "chloramphenicol", "ciprofloxacin", "doripenem", "ertapenem", "gentamicin", "imipenem", 
                  "levofloxacin", "meropenem", "nalidixic_acid", "piperacillin_tazobactam", "sulfisoxazole", 
                  "tetracycline", "tobramycin", "trimethoprim", "trimethoprim_sulfamethoxazole"),
  klebsiella = c("amikacin", "amoxicillin-clavulanicAcid", "cefotetan", "cefoxitin", "ceftazidime", "ceftriaxone", 
                 "chloramphenicol", "ciprofloxacin", "doripenem", "ertapenem", "gentamicin", "imipenem", 
                 "levofloxacin", "meropenem", "nalidixic_acid", "piperacillin_tazobactam", "tetracycline", 
                 "tobramycin", "trimethoprim", "trimethoprim_sulfamethoxazole"),
  acinetobacter = c("amikacin", "ampicillin-sulbactam", "cefepime", "ceftazidime", "ceftriaxone", "ciprofloxacin", 
                    "doripenem", "gentamicin", "imipenem", "levofloxacin", "meropenem", 
                    "piperacillin_tazobactam", "tetracycline", "tobramycin", "trimethoprim_sulfamethoxazole"),
  pseudomonas = c("amikacin", "aztreonam", "cefepime", "ceftazidime", "ciprofloxacin", "doripenem", "gentamicin", 
                  "levofloxacin", "meropenem", "piperacillin_tazobactam", "tobramycin"),
  enterobacter = c("cefepime", "ciprofloxacin", "doripenem", "ertapenem", "gentamicin", "imipenem", 
                   "levofloxacin", "meropenem", "tetracycline", "tobramycin", "trimethoprim_sulfamethoxazole"),
  staphylococcus = c("cefoxitin", "ciprofloxacin", "clindamycin", "erythromycin", "gentamicin", "levofloxacin", 
                     "moxifloxacin", "oxacillin", "penicillin", "tetracycline")
)

# Convert list to long dataframe
df_long <- bind_rows(lapply(names(data), function(genus) {
  data.frame(Genus = genus, Antibiotic = data[[genus]])
}))

# Calculate number of antibiotics per genus
genus_counts <- df_long %>%
  count(Genus)

# Ensure Genus order is fixed
df_long$Genus <- factor(df_long$Genus, levels = genus_counts$Genus)
genus_counts$Genus <- factor(genus_counts$Genus, levels = genus_counts$Genus)


# Main dot-line plot (top)
main_plot <- ggplot(df_long, aes(x = Genus, y = Antibiotic, group = Antibiotic, color = Antibiotic)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_point(size = 2) +
  scale_y_discrete(limits = rev(unique(df_long$Antibiotic))) +
  labs(x = NULL, y = "Antibiotic") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# Middle empty plot (genus labels, nudged up)
middle_plot <- ggplot(genus_counts, aes(x = Genus, y = 0)) +
  geom_blank() +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 15) +
  theme(axis.text.x = element_text(angle = 60, hjust = 0.6, margin = margin(t = -50)),
        axis.ticks.x = element_blank())

# Bar plot (below, reversed bars, smaller height, slightly lower)
bar_plot <- ggplot(genus_counts, aes(x = Genus, y = -n)) +
  geom_bar(stat = "identity", fill = "grey70", width = 0.4) +
  geom_text(aes(label = n), vjust = 1, size = 5) +  # adjust label vertical position
  labs(y = "Number of antibiotics", x = NULL) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_y_continuous(labels = function(x) abs(x))  # show positive labels

# Combine plots: top, middle (labels), bottom (bar) with reduced height
combined_plot <- main_plot / middle_plot / bar_plot + plot_layout(heights = c(4, 0.2, 1.5))

# Print combined plot
print(combined_plot)

# Save as SVG
ggsave("antibiotics_per_genus.svg", combined_plot, width = 10, height = 15, dpi = 300)

# Save as PNG
ggsave("antibiotics_per_genus.png", combined_plot, width = 10, height = 15, dpi = 300)

