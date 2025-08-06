# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Convert list to long dataframe
df_long <- bind_rows(lapply(names(data), function(genus) {
  data.frame(Genus = genus, Antibiotic = data[[genus]])
}))

# Calculate number of antibiotics per genus
genus_counts <- df_long %>%
  count(Genus)

# Ensure order is fixed (reverse for top-to-bottom)
df_long$Genus <- factor(df_long$Genus, levels = rev(genus_counts$Genus))
genus_counts$Genus <- factor(genus_counts$Genus, levels = rev(genus_counts$Genus))

# Bar plot (left, flipped left, very narrow)
bar_plot <- ggplot(genus_counts, aes(x = -n, y = Genus)) +
  geom_bar(stat = "identity", fill = "grey70", width = 0.3) +
  geom_text(aes(label = n), hjust = 0.5, vjust=0.5, size = 3.5) +
  labs(x = "Number of antibiotics", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_x_continuous(labels = function(x) abs(x))

# Main dot-line plot (right, wide)
main_plot <- ggplot(df_long, aes(x = Antibiotic, y = Genus, group = Antibiotic, color = Antibiotic)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_point(size = 2) +
  scale_x_discrete(limits = unique(df_long$Antibiotic)) +
  labs(x = "Antibiotic", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.y = element_blank(),
        legend.position = "none")

main_plot <- main_plot +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  # increase from 10 → 14
    axis.text.y = element_text(size = 14)  # increase from 10 → 14
  )

# Combine plots with narrow bar plot and wide main plot
combined_plot <- bar_plot + main_plot + plot_layout(widths = c(0.3, 1))
# Print combined plot
print(combined_plot)

# Save as SVG and PNG
ggsave("antibiotics_per_genus_leftbars.svg", combined_plot, width = 15, height = 8, dpi = 300)
ggsave("antibiotics_per_genus_leftbars.png", combined_plot, width = 15, height = 8, dpi = 300)

