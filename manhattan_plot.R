library(ggrepel)
library(tidyverse)
library(readr)
library(RColorBrewer)


# Define a custom color palette with 17 distinct colors
custom_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")

phewas_results <- fread("input_filename.csv")

phewas_results$p = as.numeric(phewas_results$p)

threshold = 0.05 / nrow(phewas_results)

phewas_results <- phewas_results %>%
  filter(OR != 1) %>%
  mutate(
    or_direction = case_when(
      OR > 1 ~ "up",
      OR < 1 ~ "down"
    )
  )

PheWAS <- ggplot(phewas_results,
                 aes(x = str_to_title(group), y = -log10(p))) +
  geom_point(aes(shape = or_direction, fill = group), size = 3, position = "jitter") +
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  scale_fill_manual(values = custom_palette) +
  theme_classic() +
  theme(
    plot.title = element_text(vjust = 0.3, hjust = .4, size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0, size = 10),
    axis.text.y = element_text(angle = 0, vjust = 0.3, hjust = 0, size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "-log10(P)", title = "FOG2 PheWAS") +
  coord_cartesian(ylim = c(0, 5), clip = "off") +  # prevents clipping
  ylim(c(0, 10)) +
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed", linewidth = .6, alpha = 0.8)



ggsave("output_filename.pdf", plot = PheWAS, width = 12, height = 10, units = "in")
