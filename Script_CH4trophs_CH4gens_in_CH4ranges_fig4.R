setwd("/...")

#Load phyloseq object

#load all libraries necessary for now and for further analyses

library(phyloseq)
library(tidyverse)
library(ggplot2) 
theme_set(theme_bw())
library(dplyr)

ps_rare

sample_data(ps_rare)

###Arrange CH4 concentrations as ranges

# Make sure CH4 is numeric
sample_data(ps_rare)$CH4_range <- as.numeric(as.character(sample_data(ps_rare)$CH4))

sample_data(ps_rare)$CH4_range <- cut(
  sample_data(ps_rare)$CH4,
  breaks = c(-Inf, 4.9, 20, 40, 80, 490, Inf),
  labels = c(
    "Atmospheric Equilibrium",
    "4.9–20",
    "20–40",
    "40–80",
    "80–490",
    "Hotspot"
  ),
  include.lowest = TRUE,
  right = TRUE
)

sample_data(ps_rare)
ps_rare@sam_data$CH4_range

ps_methanotrophs <- subset_taxa(ps_rare, Order %in% c("Hyphomicrobiales", "Methanosarcinales", "Methylacidiphilales", "Methylococcales", "Methylomirabilales"))

ps_methanotrophs <- subset_taxa(ps_methanotrophs, !(Family %in% c("A0839", "Amb-16S-1323", "Beijerinckiaceae", "Devosiaceae", "Hyphomicrobiaceae", "Hyphomicrobiales Order", "Methanosaetaceae", "Pleomorphomonadaceae", "Rhizobiaceae", "Xanthobacteraceae")))


ps_methanotrophs
tax_table(ps_methanotrophs)


df <- psmelt(ps_methanotrophs) %>%
  group_by(Sample, CH4_range) %>%
  summarise(
    Methanotroph_reads = sum(Abundance),
    .groups = "drop"
  )

ggplot(df, aes(x = CH4_range, y = Methanotroph_reads)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) + 
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "Methanotroph read counts"
  )

ch4_cols <- c(
  "Atmospheric Equilibrium" = "#1b9e77",
  "4.9–20"    = "#d95f02",
  "20–40"     = "#7570b3",
  "40–80" = "#E7298AFF",
  "80–490" = "#66A61EFF",
  "Hotspot" = "#E6AB02FF"
  
)

ggplot(df, aes(x = CH4_range, y = Methanotroph_reads, fill = CH4_range)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) + 
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "Methanotroph read counts"
  )+
  scale_fill_manual(values = ch4_cols)


#phyloseq object with methanogens

ps_methanogens <- subset_taxa(ps_rare, Genus %in% c("Methanoregula", "Methanobacterium", "Rice Cluster II Family", "Methanosphaerula", "Methanomicrobiales Order", "Methanothrix", "Methanomassiliicoccaceae Family"))

ps_methanogens

df_gens <- psmelt(ps_methanogens) %>%
  group_by(Sample, CH4_range) %>%
  summarise(
    Methanogen_reads = sum(Abundance),
    .groups = "drop"
  )


ggplot(df_gens, aes(x = CH4_range, y = Methanogen_reads, fill = CH4_range)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) + #how many samples in each category
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "Methanogen read counts"
  )+
  scale_fill_manual(values = ch4_cols)


ggplot(df_gens, aes(x = CH4_range, y = Methanogen_reads, fill = CH4_range)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "Methanogen read counts"
  )+
  scale_fill_manual(values = ch4_cols)

######### Plot proportions

# Identify methanotrophs by Genus
methanotroph_genera <- c("Methylobacter", "Crenothrix", "Methylomonadaceae Family", "Sh765B-TzT-35", "Methyloglobulus", "Candidatus Methylomirabilis", "Methylacidiphilaceae Family", "Methyloligellaceae Family", "Candidatus Methanoperedens")

# Extract OTU table and tax table
otu <- as(otu_table(ps_rare), "matrix")
tax <- as.data.frame(tax_table(ps_rare))

# Ensure samples are rows
if (taxa_are_rows(ps_rare)) {
  otu <- t(otu)
}

# Compute total reads per sample
total_reads <- rowSums(otu)

# Subset to methanotrophs
methano_otus <- rownames(tax)[tax$Genus %in% methanotroph_genera]
otu_methano <- otu[, colnames(otu) %in% methano_otus, drop = FALSE]

# Compute methanotroph reads per sample
methano_reads <- rowSums(otu_methano)

# Compute proportion
methano_prop <- methano_reads / total_reads

# Create data frame
methano_df <- data.frame(
  OriginalID = rownames(otu),
  MethanoProp = methano_prop,
  log_trophs = log10(methano_prop*100)
)
methano_df

##Import table with CH4 data
CH4 <- readr::read_tsv("/.../metadata.tsv")  
CH4
str(CH4)


# Merge
merged_df <- left_join(methano_df, CH4, by = "OriginalID")
merged_df


#Order the CH4 ranges

merged_df$Methane_range <- factor(merged_df$Methane_range, levels = c("Atmospheric Equilibrium", "4.9-20", "20-40", "40-80", "80-490", "Hotspot"))


ch4_cols2 <- c(
  "Atmospheric equilibrium" = "#1b9e77",
  "4.9-20"    = "#d95f02",
  "20-40"     = "#7570b3",
  "40-80" = "#E7298AFF",
  "80-490" = "#66A61EFF",
  "Hotspot" = "#E6AB02FF"
  
)


ggplot(merged_df, aes(x = Methane_range, y = MethanoProp, fill = Methane_range)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) + #how many samples in each category
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "Relative abundance of methanotrophs (%)"
  )+
  scale_fill_manual(values = ch4_cols2, name = expression("CH"[4]~"range (nM)"))+
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

library(stringr)

trophs_abd <- ggplot(merged_df, aes(x = Methane_range, y = log_trophs, fill = Methane_range)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) + #how many samples in each category
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "% Methanotrophs\n(log10 scale)"
  )+
  scale_fill_manual(values = ch4_cols2)+
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 15))
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))#plot the atm equil in 2 lines
  
trophs_abd
 

####Methanogens

# Identify methanogens by Genus
methanogenic_genera <- c("Methanoregula", "Methanobacterium", "Rice Cluster II Family", "Methanosphaerula", "Methanomicrobiales Order", "Methanothrix", "Methanomassiliicoccaceae Family")

# Subset to methanogens
methanogenic_otus <- rownames(tax)[tax$Genus %in% methanogenic_genera]
otu_methanogens <- otu[, colnames(otu) %in% methanogenic_otus, drop = FALSE]

# Compute methanogen reads per sample
methanogen_reads <- rowSums(otu_methanogens)

# Compute proportion
methanogen_prop <- methanogen_reads / total_reads

# Create data frame
methanogen_df <- data.frame(
  OriginalID = rownames(otu),
  MethanogenProp = methanogen_prop,
  log_gens = log10(methanogen_prop*100)
)
methanogen_df

# Merge
merged_df_gens <- left_join(methanogen_df, CH4, by = "OriginalID")
merged_df_gens

#Order the CH4 ranges

merged_df_gens$Methane_range <- factor(merged_df_gens$Methane_range, levels = c("Atmospheric Equilibrium", "4.9-20", "20-40", "40-80", "80-490", "Hotspot"))


ch4_cols2 <- c(
  "Atmospheric Equilibrium" = "#1b9e77",
  "4.9-20"    = "#d95f02",
  "20-40"     = "#7570b3",
  "40-80" = "#E7298AFF",
  "80-490" = "#66A61EFF",
  "Hotspot" = "#E6AB02FF"
  
)

gens_abd <- ggplot(merged_df_gens, aes(x = Methane_range, y = log_gens, fill = Methane_range)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) + #how many samples in each category
  theme_bw() +
  labs(
    x = expression(CH[4]~range~"(nM)"),
    y = "% Methanogens\n(log10 scale)"
  )+
  scale_fill_manual(values = ch4_cols2)+
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 15))
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) #plot the atm equil in 2 lines

gens_abd


#combine plots

library(patchwork)
print(trophs_abd + gens_abd)

library(cowplot)
library(ggpubr)


p_combined <- ggarrange(gens_abd, trophs_abd, labels = c("A", "B"), ncol = 1, nrow = 2)
p_combined

