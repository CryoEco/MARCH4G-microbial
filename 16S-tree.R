setwd("/...")

load("final_march4g2.RData")

#Plotting and exploring the dada2 object

library(phyloseq)
library(readr)
library(dada2)

ASV=otu_table(asv_tab_sampleNames, taxa_are_rows = TRUE)

TAX=tax_table(taxa2)
ASV
TAX
TAX[222,]

#Create a phylogenetic tree before creating the phyloseq object
#code taken from https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#Paper ref: https://f1000research.com/articles/5-1492/v1#ref-3
#install package DECIPHER to align seqs and phangorn to construct tree

library(DECIPHER)
library(phangorn)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE )

phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

save.image("/.../16S-tree.RData")
