library(tidyverse)

circle <- function(center_x, center_y, radius, label = "", resolution = 1000) {
  x <- seq(-radius, radius, length.out = resolution)
  y_upper <- sqrt(radius^2 - x^2)
  y_lower <- -sqrt(radius^2 - x^2)
  
  tibble(x = x + center_x,
         ymin = y_lower + center_y,
         ymax = y_upper + center_y,
         label = label)
}

bind_rows(circle(0, 0, 1, "A"),
          circle(1, 0, 1, "B"),
          circle(0.5, -0.9, 1, "C")) %>%
  ggplot(aes(x = x, ymin = ymin, ymax = ymax, group = label, fill = label, alpha = label)) +
  geom_ribbon(color = "black", linewidth = 1, show.legend = FALSE) +
  annotate("text", x = c(-0.5, 0.5, 1.5), y = c(0.2, 0.4, 0.2),
           label = c("1,105", "32", "35"), fontface = "bold") +
  annotate("text", x = c(-0.1, 0.5, 1.1), y = c(-0.6, -0.3, -0.6),
           label = c("782", "211", "14"), fontface = "bold") +
  annotate("text", x = 0.5, y = -1.3,
           label = "399", fontface = "bold") +
  annotate("text", x = c(-0.5, 0.5, 1.5), y = c(1.2, -2, 1.2),
           label = c("ReferemceGemeCatalog\n8,157",
                     "CARD\n5,078",
                     "ResFinder\n3,150"),
           fontface = "bold", lineheight = 0.8, vjust = 1) +
  coord_equal() +
  scale_fill_manual(breaks = c("A", "B", "C"),
                    values = c("#96dee8", "#ff9f7f", "#3fb1e3")) +
  scale_alpha_discrete(breaks = c("A", "B", "C"),
                       range = c(1, 0.25)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

c("#3fb1e3", "#6be6c1", "#626c91", "#a0a7e6", "#c4ebad", "#ff9f7f", "#96dee8")
