setwd("/Users/...")

#Load dada2 data

load("final_march4g2.RData")

#Plotting and exploring the dada2 object
#load all libraries necessary for now and for further analyses

library(phyloseq)
library(readr)
library(dada2)
library(tidyverse)
library(vegan)
library(dendextend)
library(viridis)  
library(ggplot2) 
theme_set(theme_bw())
library(microViz)
library(shiny)
library(Biostrings)
library(ape)
library(phytools)
library(microbiome)


ASV=otu_table(asv_tab_sampleNames, taxa_are_rows = TRUE)

TAX=tax_table(taxa2)
ASV
TAX
TAX[222,]

metadata <- readr::read_tsv("metadata.tsv")
class(metadata)
row.names(metadata)
colnames(metadata)
metadata
#row names did not match sample names, they were in numeric order, so had to set column as row names:
metadata_corrected = metadata %>% remove_rownames %>% column_to_rownames(var = "OriginalID")

row.names(metadata_corrected)

sampledata = sample_data(metadata_corrected)
sampledata
dim(sampledata)
class(sampledata)

#load tree
load("16S-tree.RData")

tree <- fitGTR$tree
tree

all(rownames(ASV) %in% taxa_names(tree)) # should return TRUE
all(taxa_names(tree) %in% rownames(ASV)) # should return TRUE
#it returns false, showing that the names in the ASV table don't match with the names in the tree

# Check which taxa names are in ASV but not in the tree
setdiff(rownames(ASV), taxa_names(tree))

# Check which taxa names are in the tree but not in ASV
setdiff(taxa_names(tree), rownames(ASV))


#Create the phyloseq object
physeq = phyloseq(ASV, TAX, sampledata, tree)

physeq

sample_names(physeq)

samdat_tbl(physeq)
sample_variables(physeq)
otu_get(physeq)
nsamples(physeq)
ntaxa(physeq)
phy_tree(physeq)

sample_names(physeq)
taxa_names(physeq)

tax_table(physeq) [1:5, ]
#look at the taxonomy table interactively
tax_fix_interactive(physeq)

tax_table(physeq) [1:56, ]

#remove contaminants, including chloroplast and mitochondria. Also remove NA phylum
asvs_to_remove <- readLines("conts-NAs.txt")
taxa_decontam <- taxa2[!(rownames(taxa2) %in% asvs_to_remove),]
dim(taxa_decontam)
dim(taxa)

#Do the same for asv_tab_sampleNames
asv_tab_decontam <- asv_tab_sampleNames[!(rownames(asv_tab_sampleNames) %in% asvs_to_remove), ]
dim(asv_tab_decontam)
dim(asv_tab_sampleNames)

#Now recreate phyloseq object without the contaminants

ASV=otu_table(asv_tab_decontam, taxa_are_rows=TRUE)
TAX=tax_table(taxa_decontam)
asv_tab_decontam
taxa_decontam
ps = phyloseq(ASV, TAX, sampledata, tree)


ps
sample_names(ps)
sample_data(ps)
taxa_names(ps)
ntaxa(ps)
otu_table(ps)
otu_get(ps)
tax_fix_interactive(ps)

#before using subset_taxa it is important to remove all the NAs from the ASV table, otherwise the function will delete/omit ASVs that have NA, rendering the analysis wrong 

ps <- tax_fix(ps, min_length = 3,
              unknowns = c("NA", "Incertae Sedis", "Hyphomicrobiales Incertae Sedis", "Gammaproteobacteria Incertae Sedis"),
              sep = " ", anon_unique = TRUE,
              suffix_rank = "classified")

tax_table(ps) [1:56, ]

ps
sample_data(ps)
tax_fix_interactive(ps) #names are fixed

summarize_phyloseq(ps)
reads_sample <- readcount(ps)
reads_sample

#remove blanks, positive controls,  samples that are not to be included in the analyses

ps_sterivexes <- subset_samples(ps, OrSiteName != "UP-ICE" & ConsSiteName != "KA-05-b" & ConsSiteName != "KA-04-c" & Sample_or_Control != "PositiveControl" & ConsSiteName != "KA-08" & ConsSiteName != "KA-07" & Sample_or_Control != "Control" & OrSiteName != "KA-03-A" & OrSiteName != "KA-05-B" & OrSiteName != "NU02-C")

ps_sterivexes

phyloseq_validate(ps_sterivexes, remove_undetected = TRUE, verbose = TRUE)
ps_sterivexes

#remove taxa that have 0 counts 
ps_sterivex2 = ps_sterivexes
ps_sterivex2 = prune_taxa(taxa_sums(ps_sterivexes) >= 1, ps_sterivexes)
ps_sterivex2

summarize_phyloseq(ps_sterivex2)

reads_sample2 <- readcount(ps_sterivex2)
reads_sample2

#Samples 16, 17, 25, 34, 61, 62 have fewer than 200 reads - before discarding them from the phyloseq object, check what ASVs were assigned to them

tax_fix_interactive(ps_sterivex2)

ps_sterivex2 %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(tax_level = "Genus", sample_order = "asis", n_taxa = 41)+
  facet_wrap(facets = vars(ConsSiteName), scales = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Microbial composition by location",
    x = "Samples", y = "Relative Abundance"
  ) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) 


ps_sterivex3 <- subset_samples(ps_sterivex2, ConsSiteName != "KA-03-a" & OrSiteName != "KA-05-C" & ConsSiteName != "NU-02" & OrSiteName != "QA-5-1")

ps_sterivex3
phyloseq_validate(ps_sterivex3, remove_undetected = TRUE, verbose = TRUE)
ps_sterivex3
sample_data(ps_sterivex3)

#remove taxa that have 0 counts 
ps_sterivex4 = ps_sterivex3
ps_sterivex4 = prune_taxa(taxa_sums(ps_sterivex3) >= 1, ps_sterivex3)
ps_sterivex4

summarize_phyloseq(ps_sterivex4)

reads_sample <- readcount(ps_sterivex4)
reads_sample

write.table(tax_table(ps_sterivex4), "/Users/.../tax_table_fulldataset-clean.txt", sep = "\t", quote = F, col.names = NA)
write.table(otu_table(ps_sterivex4), "/Users/.../asv_table_fulldataset-clean.txt", sep = "\t", quote = F, col.names = NA)


#rarefaction curve before rarefying dataset
taxa_are_rows(ps_sterivex4)
mat <- t(otu_table(ps_sterivex4))
class(mat) <- "matrix"
class(mat)
mat <- as(t(otu_table(ps_sterivex4)), "matrix")
class(mat)
raremax <- min(rowSums(mat))
rarecurve(mat, step = 500, sample = raremax,xlab = "Sequencing depth", ylab = "Number of ASVs", col = "blue")
rarecurve(mat, step = 500, sample = raremax,xlab = "Sequencing depth", ylab = "Number of ASVs", col = "blue", xlim=c(0,15000))


#rarefy reads to 7678 (the smallest number of reads in the dataset) - majority os samples didn't increase the number of ASVs after this threshold
ps_rare <- rarefy_even_depth(ps_sterivex4, sample.size = 7678, rngseed = 500, replace = FALSE)

ps_rare

sample_sums(ps_rare)

tax_table(ps_rare)
otu_table(ps_rare)

tax_table(ps_rare) [1400:1450, ]
tax_fix_interactive(ps_rare)


###Root tree in ps_rare
ps_rare

##check if tree is rooted
is.rooted(phy_tree(ps_rare)) 

#Try to root tree choosing an outgroup based on the longest tree branch terminating in a tip

pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }

tree_rare <- phy_tree(ps_rare)
out.group_rare <- pick_new_outgroup(tree_rare)

out.group_rare ##ASV_758

new.tree_rare <- ape::root(tree_rare, outgroup=out.group_rare, resolve.root=TRUE)
phy_tree(ps_rare) <- new.tree_rare
phy_tree(ps_rare) ##

##check if tree is rooted
is.rooted(phy_tree(ps_rare)) 

#save tree
tree_rooted_rare = phy_tree(ps_rare)
ape::write.tree(tree_rooted_rare, "/Users/.../tree_rooted_rare.nwk")

###Root tree in ps_sterivex4
is.rooted(phy_tree(ps_sterivex4))

tree_ps4 <- phy_tree(ps_sterivex4)
out.group_ps4 <- pick_new_outgroup(tree_ps4)

out.group_ps4 ##ASV_959

new.tree_ps_4 <- ape::root(tree_ps4, outgroup=out.group_ps4, resolve.root=TRUE)
phy_tree(ps_sterivex4) <- new.tree_ps_4
phy_tree(ps_sterivex4) ##

##check if tree is rooted
is.rooted(phy_tree(ps_sterivex4)) 

tree_rooted_fulldataset = phy_tree(ps_sterivex4)
ape::write.tree(tree_rooted_fulldataset, "/Users/.../tree_rooted_fulldataset.nwk")

write.table(tax_table(ps_rare), "/Users/.../tax_table_rarefieddataset.txt", sep = "\t", quote = F, col.names = NA)
write.table(otu_table(ps_rare), "/Users/.../asv_table_rarefieddataset.txt", sep = "\t", quote = F, col.names = NA)


ps_sterivex4 %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(tax_level = "Genus", sample_order = "asis", n_taxa = 41)+
  facet_wrap(facets = vars(Location), scales = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Microbial composition by location",
    x = "Samples", y = "Relative Abundance"
  ) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

ps_sterivex4 %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(tax_level = "Genus", sample_order = "asis", n_taxa = 41)+
  facet_wrap(facets = vars(Site), scales = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Microbial composition by location",
    x = "Samples", y = "Relative Abundance"
  ) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) 


ps_rare %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(tax_level = "Genus", sample_order = "asis", n_taxa = 41)+
  facet_wrap(facets = vars(Site), scales = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Microbial composition by location",
    x = "Samples", y = "Relative Abundance"
  ) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))


########### Prepare new phyloseq object for iCAMP analyses - average triplicates from the same sampling location
##Average from object ps_sterivex4 and only then rarefy

sample_data(ps_sterivex4)

ps_merged <- merge_samples(ps_sterivex4, group = "ConsSiteName", fun = mean)
ps_merged
otu_table(ps_merged)
otu_get(ps_merged)


ps_merged %>%
  tax_agg("Genus") %>%
  ps_seriate(dist = "bray", method = "OLO_ward") %>%
  comp_barplot(tax_level = "Genus", sample_order = "asis", n_taxa = 41) +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Microbial composition by location",
    x = "Samples", y = "Relative Abundance"
  ) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

nsamples(ps_sterivex4)          # original number of samples (with replicates)
nsamples(ps_merged)   # number of sites (unique methane values)

sample_data(ps_merged)$CH4
sample_data(ps_merged)
sample_names(ps_merged)

#Correct the metadata table 
md <- as.data.frame(sample_data(ps_sterivex4))

library(dplyr)

md_avg <- md %>%
  group_by(ConsSiteName) %>%
  summarise(across(where(is.numeric), mean),
            across(where(is.character), ~ first(.x)))
md_avg

md_avg_df <- as.data.frame(md_avg)
rownames(md_avg_df) <- md_avg_df$ConsSiteName

rownames(md_avg_df)

sample_names(ps_merged)
rownames(md_avg_df)

sample_data(ps_merged) <- sample_data(md_avg_df)
ps_merged
sample_data(ps_merged)

#Correct the DOC value in IL-02-a
ps_merged@sam_data[["IL-02-a", "DOC"]] <- 18

sample_data(ps_merged)
ps_merged

summarize_phyloseq(ps_merged)

reads_merged <- readcount(ps_merged)
reads_merged

#rarefaction curve before rarefying dataset
taxa_are_rows(ps_merged)
mat2 <- t(otu_table(ps_merged))
class(mat2) <- "matrix"
class(mat2)
mat2 <- as(t(otu_table(ps_merged)), "matrix")
class(mat2)
raremax2 <- min(rowSums(mat2))
rarecurve(mat2, step = 500, sample = raremax2,xlab = "Sequencing depth", ylab = "Number of ASVs", col = "blue")
rarecurve(mat2, step = 500, sample = raremax2,xlab = "Sequencing depth", ylab = "Number of ASVs", col = "blue", xlim=c(0,15000))


#rarefy reads to 50225 (the smallest number of reads in the dataset) 
ps_merged_rare <- rarefy_even_depth(ps_merged, sample.size = 50225, rngseed = 500, replace = FALSE)

ps_merged_rare

sample_sums(ps_merged_rare)

tax_table(ps_merged_rare)
otu_table(ps_merged_rare)

tax_table(ps_merged_rare) [1400:1450, ]

write.table(tax_table(ps_merged_rare), "/Users/.../tax_table_merged_sites_Rarefied.txt", sep = "\t", quote = F, col.names = NA)
write.table(otu_table(ps_merged_rare), "/Users/.../asv_table_merged_sites_Rarefied.txt", sep = "\t", quote = F, col.names = NA)

otu_table_mat <- as.matrix(otu_table(ps_merged_rare))
otu_table_mat

##check if tree is rooted
is.rooted(phy_tree(ps_merged_rare)) 

#save tree
tree_merged_rare = phy_tree(ps_merged_rare)
ape::write.tree(tree_merged_rare, "/Users/.../tree_merged_rare.nwk")


save.image("physeq-merged-replicates-rareAndFull-dataset.RData")
