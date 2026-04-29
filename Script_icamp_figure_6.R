setwd("/...")

library(dplyr)
library(ggplot2)
library(tidyr)

##Overall processes
df <- read.csv("CH4-continuous.ProcessImportance_EachTurnover.csv")

process_columns <- c("HeS", "HoS", "DL", "HD", "DR")
overall_contributions <- colMeans(df[, process_columns])

print(overall_contributions)


df_long <- df %>%
  select(all_of(process_columns)) %>%
  pivot_longer(cols = everything(), names_to = "Process", values_to = "Contribution") %>%
  mutate(Contribution = Contribution * 100)  # Convert to %

# Rename levels for readability
df_long$Process <- factor(df_long$Process,
                          levels = c("HeS", "HoS", "DL", "HD", "DR"),
                          labels = c("Heterogeneous selection", "Homogeneous selection", "Dispersal limitation", "Homogenizing dispersal", "Drift"))

# Plot
p1 <- ggplot(df_long, aes(x = Process, y = Contribution)) +
  geom_boxplot(fill = "#4292c6", color = "black", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white") +  # diamond mean
  labs(
    x = "Assembly processes",
    y = "Contribution (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
p1


###With different colors for each process

fill_colors_all <- c(
  "Homogeneous selection" = "#1f78b4",
  "Heterogeneous selection" = "#a6cee3",
  "Dispersal limitation" = "#fdbf6f",
  "Drift" = "#e31a1c",
  "Homogenizing dispersal" = "#fb9a99"
)

# Plot
p2 <- all_procs <- ggplot(df_long, aes(x = Process, y = Contribution, fill = Process)) +
  geom_boxplot(color = "black", outlier.shape = 16) +
  scale_fill_manual(values = fill_colors_all) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +  # Add line breaks
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white") +
  labs(
    x = "Assembly processes",
    y = "Contribution (%)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 17),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", size = 19, hjust = 0.5, margin = margin(t=20)),
    axis.title.y = element_text(face = "bold", size = 19),
    legend.position = "none"
  )

all_procs
p2


###Processes contributing to community assembly in each methane range

df2 <- read.delim("CH4-continuous.iCAMP.BootSummary.CH4cat3.txt", header = TRUE, sep = "\t") %>%
  filter(!grepl("vs", Group))

# Reshape to long format
df_long2 <- df2 %>%
  select(Group, Process, Mean) %>%
  pivot_wider(names_from = Process, values_from = Mean) %>%
  mutate(
    `Homogeneous selection` = Homogeneous.Selection,
    `Heterogeneous selection` = Heterogeneous.Selection,
    `Dispersal limitation` = Dispersal.Limitation,
    Drift = Drift.and.Others,
    `Homogenizing dispersal` = Homogenizing.Dispersal
  ) %>%
  select(Group, `Homogeneous selection`, `Heterogeneous selection`,
         `Dispersal limitation`, Drift, `Homogenizing dispersal`) %>%
  pivot_longer(cols = -Group, names_to = "Process", values_to = "Contribution")%>%
  mutate(Contribution = Contribution * 100)  # Convert to %

# Define all levels including fake headers
process_levels <- c(
  "Deterministic processes",    #fake
  "Heterogeneous selection",
  "Homogeneous selection",
  "Stochastic processes",             # fake
  "Dispersal limitation",
  "Homogenizing dispersal",
  "Drift"
)

# Assign full factor with unused levels
df_long2$Process <- factor(df_long2$Process, levels = process_levels)

# Add fake rows for legend grouping
fake_legend_rows <- expand.grid(
  Group = unique(df_long2$Group),
  Process = c("Deterministic processes", "Stochastic processes")
)
fake_legend_rows$Contribution <- 0
fake_legend_rows$Process <- factor(fake_legend_rows$Process, levels = process_levels)

# Combine real and fake entries
df_plot <- bind_rows(df_long2, fake_legend_rows)

df_plot$Group <- factor(
  df_plot$Group,
  levels = c(
    "Atmospheric Equilibrium",
    "Range_4.9_20",
    "Range_20_40",
    "Range_40_80",
    "Range_80_490",
    "Hotspot"
  ),
  labels = c(
    "Atm. Equilibrium",
    "4.9–20",
    "20–40",
    "40–80",
    "80–490",
    "Hotspot"
  )
)


# Color scale (NA for fake legend headers)
fill_colors <- c(
  "Deterministic processes" = NA,
  "Heterogeneous selection" = "#a6cee3",
  "Homogeneous selection" = "#1f78b4",
  "Stochastic processes" = NA,
  "Dispersal limitation" = "#fdbf6f",
  "Homogenizing dispersal" = "#fb9a99",
  "Drift" = "#e31a1c"
)

library(ggtext)

p3 <- ggplot(df_plot, aes(x = Group, y = Contribution, fill = Process)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = fill_colors,
    drop = FALSE,
    guide = guide_legend(
      title = NULL,
      override.aes = list(fill = c(
        NA, # Deterministic processes (fake)
        "#a6cee3",             # Heterogeneous selection
        "#1f78b4",             # Homogeneous selection             
        NA,                    # Stochastic processes (fake)
        "#fdbf6f",            # Dispersal limitation
        "#fb9a99",              # Homogenizing dispersal
        "#e31a1c"             # Drift
      ))
    )
  ) +
  labs(
    x = "<b>CH<sub>4</sub> range</b>",
    y = "Contribution (%)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 17),
    legend.text = element_text(face = "plain", size = 16),
    legend.key.height = unit(0.8, "cm"),
    legend.spacing.y = unit(0.6, "cm"),
    legend.margin = margin(5,5,5,5),
    axis.title.x = ggtext::element_markdown(size = 19,hjust = 0.5, margin = margin(t=20)),
    axis.title.y = element_text(face = "bold", size = 19),
    plot.margin = margin(10, 10, 30, 25)
  )
p3


library(patchwork)
library(ggpubr)
print(p2 + p3)

p_combined <- ggarrange(p2, p3, labels = c("A", "B"), ncol = 1, nrow = 2)
p_combined

####Shared legends for both plots

library(cowplot)

p3 <- p3 + theme(plot.margin = margin(10, 10, 10, 25))

p3_noleg <- p3 + theme(legend.position = "none")

legend_p3 <- get_legend(
  p3 + theme(legend.position = "right")
)

plots_col <- plot_grid(
  p2, p3_noleg,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

final_graph <- plot_grid(
  plots_col, legend_p3,
  ncol = 2,
  rel_widths = c(1, 0.28)
)

final_graph


final_graph2 <- plot_grid(   #graph with 2 columns and 1 row 
  p2, p3,
  labels = c("A", "B"),
  rel_widths = c(1, 1),
  align = "h"
)

final_graph2
