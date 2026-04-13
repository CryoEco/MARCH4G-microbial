
#### DADA2 pipeline

path1 <-"/storage/praha1/home/liawtz/final-sterivex-march4g/run1"
path2 <-"/storage/praha1/home/liawtz/final-sterivex-march4g/run1_repeat"
path3 <-"/storage/praha1/home/liawtz/final-sterivex-march4g/run2"
list.files(path1)
list.files(path2)
list.files(path3)

#set up working environment, install packages and load libraries needed
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#Separate the steps per run
#####Run 1 

fnFs_run1 <- sort(list.files(path1, pattern = "_R1_trimmed.fq.gz", full.names = TRUE))
fnRs_run1 <- sort(list.files(path1, pattern = "_R2_trimmed.fq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name_1 <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names_run1 <- unname(sapply(fnFs_run1, get.sample.name_1))
head(sample.names_run1)

##Filter and trim
filtFs_run1 <- file.path(path1, "filtered", basename(fnFs_run1))
filtRs_run1 <- file.path(path1, "filtered", basename(fnRs_run1))

out_run1 <- filterAndTrim(fnFs_run1, filtFs_run1, fnRs_run1, filtRs_run1, maxN = 0, maxEE = c(2,2), minLen = 200, truncQ = 2, truncLen=c(223, 222), rm.phix = TRUE, compress = TRUE, multithread = FALSE)  
head(out_run1)

#Learn the error rates
errF_run1 <- learnErrors(filtFs_run1, multithread = TRUE)
errR_run1 <- learnErrors(filtRs_run1, multithread = TRUE)

#Dereplication
derepFs_run1 <- derepFastq(filtFs_run1, verbose = TRUE)
derepRs_run1 <- derepFastq(filtRs_run1, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs_run1) <- sample.names_run1
names(derepRs_run1) <- sample.names_run1

#Sample inference
dadaFs_run1 <- dada(derepFs_run1, err = errF_run1, multithread = TRUE)
dadaRs_run1 <- dada(derepRs_run1, err = errR_run1, multithread = TRUE)

dadaFs_run1[1]

#Merge paired reads
mergers_run1 <- mergePairs(dadaFs_run1, derepFs_run1, dadaRs_run1, derepRs_run1, verbose=TRUE)

head(mergers_run1[[1]])

#Construct sequence table
seqtab_run1 <- makeSequenceTable(mergers_run1)
dim(seqtab_run1)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_run1))) #do all contigs have a length of 250bp or around this?

#Remove seqs that are much longer or shorter than expected - amplicon is 253 bp
seqtab_run1_2 <- seqtab_run1[,nchar(colnames(seqtab_run1)) %in% seq(250,256)]

#inspect distribution of sequence lengths from new table
table(nchar(getSequences(seqtab_run1_2)))


#####Run 2 - run1_repeat

fnFs_run2 <- sort(list.files(path2, pattern = "_R1_trimmed.fq.gz", full.names = TRUE))
fnRs_run2 <- sort(list.files(path2, pattern = "_R2_trimmed.fq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name_2 <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names_run2 <- unname(sapply(fnFs_run2, get.sample.name_2))
head(sample.names_run2)

##Filter and trim
filtFs_run2 <- file.path(path2, "filtered", basename(fnFs_run2))
filtRs_run2 <- file.path(path2, "filtered", basename(fnRs_run2))

out_run2 <- filterAndTrim(fnFs_run2, filtFs_run2, fnRs_run2, filtRs_run2, maxN = 0, maxEE = c(2,2), minLen = 200, truncQ = 2, truncLen=c(223,222), rm.phix = TRUE, compress = TRUE)  
head(out_run2)

#Learn the error rates
errF_run2 <- learnErrors(filtFs_run2, multithread = TRUE)
errR_run2 <- learnErrors(filtRs_run2, multithread = TRUE)

#Dereplication
derepFs_run2 <- derepFastq(filtFs_run2, verbose = TRUE)
derepRs_run2 <- derepFastq(filtRs_run2, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs_run2) <- sample.names_run2
names(derepRs_run2) <- sample.names_run2

#Sample inference
dadaFs_run2 <- dada(derepFs_run2, err = errF_run2, multithread = TRUE)
dadaRs_run2 <- dada(derepRs_run2, err = errR_run2, multithread = TRUE)

dadaFs_run2[1]

#Merge paired reads
mergers_run2 <- mergePairs(dadaFs_run2, derepFs_run2, dadaRs_run2, derepRs_run2, verbose=TRUE)

head(mergers_run2[[1]])

#Construct sequence table
seqtab_run2 <- makeSequenceTable(mergers_run2)
dim(seqtab_run2)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_run2))) #do all contigs have a length of 250bp or around this?

#Remove seqs that are much longer or shorter than expected - amplicon is 253 bp
seqtab_run2_2 <- seqtab_run2[,nchar(colnames(seqtab_run2)) %in% seq(250,256)]

#inspect distribution of sequence lengths from new table
table(nchar(getSequences(seqtab_run2_2)))

#####Run 3 - namely run2

fnFs_run3 <- sort(list.files(path3, pattern = "_R1_trimmed.fq.gz", full.names = TRUE))
fnRs_run3 <- sort(list.files(path3, pattern = "_R2_trimmed.fq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name_3 <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names_run3 <- unname(sapply(fnFs_run3, get.sample.name_3))
head(sample.names_run3)

##Filter and trim
filtFs_run3 <- file.path(path3, "filtered", basename(fnFs_run3))
filtRs_run3 <- file.path(path3, "filtered", basename(fnRs_run3))

out_run3 <- filterAndTrim(fnFs_run3, filtFs_run3, fnRs_run3, filtRs_run3, maxN = 0, maxEE = c(2,2), minLen = 200, truncQ = 2, truncLen=c(223,222), rm.phix = TRUE, compress = TRUE)  
head(out_run3)

#Learn the error rates
errF_run3 <- learnErrors(filtFs_run3, multithread = TRUE)
errR_run3 <- learnErrors(filtRs_run3, multithread = TRUE)

#Dereplication
derepFs_run3 <- derepFastq(filtFs_run3, verbose = TRUE)
derepRs_run3 <- derepFastq(filtRs_run3, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs_run3) <- sample.names_run3
names(derepRs_run3) <- sample.names_run3

#Sample inference
dadaFs_run3 <- dada(derepFs_run3, err = errF_run3, multithread = TRUE)
dadaRs_run3 <- dada(derepRs_run3, err = errR_run3, multithread = TRUE)

dadaFs_run3[1]

#Merge paired reads
mergers_run3 <- mergePairs(dadaFs_run3, derepFs_run3, dadaRs_run3, derepRs_run3, verbose=TRUE)

head(mergers_run3[[1]])

#Construct sequence table
seqtab_run3 <- makeSequenceTable(mergers_run3)
dim(seqtab_run3)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_run3))) #do all contigs have a length of 250bp or around this?

#Remove seqs that are much longer or shorter than expected - amplicon is 253 bp
seqtab_run3_2 <- seqtab_run3[,nchar(colnames(seqtab_run3)) %in% seq(250,256)]

#inspect distribution of sequence lengths from new table
table(nchar(getSequences(seqtab_run3_2)))

#Merge sequence tables
seqtab_all <- mergeSequenceTables(seqtab_run1_2, seqtab_run2_2, seqtab_run3_2)
seqtab_all[1]

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab_all, method = "consensus", multithread=TRUE, verbose = TRUE)

sum(seqtab.nochim)/sum(seqtab_all)

# bimeras out of  seqs
dim(seqtab.nochim)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))



#Assigning taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread = TRUE)


#Inspect the assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print[1:56]
taxa


#Extracting the info from DADA2
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_seqs
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
asv_headers #do not include in batch job
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
asv_headers


# making and writing out a fasta of my final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_final_march4g2.fa")

#count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASV_counts_final_march4g2.tsv", sep = "\t", quote = F, col.names = NA)

#create another count table but with sample names as column names (instead of 01_filtered.fq)
asv_tab_sampleNames <- t(seqtab.nochim)
row.names(asv_tab_sampleNames) <- sub(">", "", asv_headers)
colnames(asv_tab_sampleNames) <- gsub("_filtered.fq", "", colnames(asv_tab_sampleNames))
asv_tab_sampleNames

write.table(asv_tab_sampleNames, "seq_counts_final_march4g2.tsv", sep = "\t", quote = F, col.names = NA)

#create a count table with sample names as column names and sequences instead of 'ASV_01'
seq_tab_sampleNames <- t(seqtab.nochim)
colnames(seq_tab_sampleNames) <- gsub("_filtered.fq", "", colnames(seq_tab_sampleNames))
seq_tab_sampleNames

write.table(seq_tab_sampleNames, "seq_counts_sampleNames_final_march4g2.tsv", sep = "\t", quote = F, col.names = NA)

#export the tax table
#with the ASVs sequences
write.table(taxa, "ASVs_tax_seqs_final_march4g2.tsv", sep = "\t", quote = F, col.names = NA)

#with ASVs as the names
taxa2 <- taxa
taxa2
row.names(taxa2) <- sub(">", "", asv_headers)
taxa2[222,]
write.table(taxa2, "ASVs_tax-names-final_march4g2.tsv", sep = "\t", quote = F, col.names = NA)


##export tax table with taxanomy
write.csv(cbind(t(seqtab.nochim), taxa), "new-count-tax-table.csv", quote = FALSE)

save.image("/storage/praha1/home/liawtz/final-sterivex-march4g/final_march4g2.RData")
