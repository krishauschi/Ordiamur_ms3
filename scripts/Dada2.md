#load required packages
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)

#upload raw amplicon sequence files for all samples in the format "sample_name_R1_001.fq.gz" (forward) and "sample_name_R2_001.fq.gz" (reverse) to the folder you are working in. 

path <- "~/my/path" # change to the right location
list.files(path)

#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fq.gz and SAMPLENAME_R2_001.fq.gz
fnFs <- sort(list.files(path, pattern = "_R1_001.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fq.gz", full.names = TRUE))
length(fnFs) #shows you the number of forward files --> if labelling has errors for one or more samples, those won't appear here
length(fnRs) #shows you the number of reverse files

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

#plotQualityProfile(system.file, n = 5e+05, aggregate = FALSE)
pdf("~/my/path/plot_Quality_Forward.pdf")
plotQualityProfile(fnFs[1:6])
dev.off()

pdf("~/my/path/plot_Quality_Reverse.pdf")
plotQualityProfile(fnRs[1:6])
dev.off()

#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))

#Filtering on Windows,set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress = TRUE, truncQ = 2, minLen=100, trimLeft=c(17,20), maxN=0, maxEE=c(2,2), rm.phix=TRUE, matchIDs = TRUE, multithread=TRUE)
head(out)

#plot quality after filter
pdf("~/my/path/plot_Quality_Forward2.pdf")
plotQualityProfile(filtFs[1:6])
dev.off()

pdf("~/my/path/plot_Quality_Reverse2.pdf")
plotQualityProfile(filtRs[1:6])
dev.off()

#set run
set.seed(2022)  #set any number for to initilize the random seed. Note this number.

#predict error and plot
errF <- learnErrors(filtFs, multithread=TRUE)
pdf("~/my/path/plot_error_Forward.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

errR <- learnErrors(filtRs, multithread=TRUE)
pdf("~/my/path/plot_error_Reverse.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#dereplicate sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#run dada
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

save.image(file = "~/KH10/KH10_dada_16S/dada2_output/KH10_ITSseq.RData")
load("~/KH10/KH10_dada_16S/dada2_output/KH10_ITSseq.RData")

#merge paired sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)  # first number is samples, second number is ASVs 
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#get sequences, remove chimera, and write tables
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)                # first number is samples, second number is ASVs
sum(seqtab.nochim)/sum(seqtab)     
seqtab.nochim.t<-t(seqtab.nochim)        

#Export results
write.csv(seqtab.nochim.t, "~/my/path/seq_Phy_Table.csv")

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

##If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "~/my/path/seq_Quality_Control.csv")

saveRDS(seqtab.nochim, file = "~/my/path/seqtab.nochim_seq.rds")
save.image(file = "~/my/path/dada2seq.RData")

#for 16S
#Assign Taxonomy at species level based on the minBoot of bootstrap confidence - minBoot=50 is default
taxa <- assignTaxonomy(seqtab.nochim, "~/my/path/database/silva_nr99_v138.1_train_set.fa.gz", minBoot = 0, outputBootstraps = TRUE, multithread=TRUE)
write.csv(taxa, "~/my/path/taxa_silva_v138.1_table.csv")

#for ITS
#Assign Taxonomy at species level based on the minBoot of bootstrap confidence - minBoot=50 is default
taxa <- assignTaxonomy(seqtab.nochim, "~/my/path/sh_general_release_dynamic_04.04.2024.fasta", minBoot = 0, outputBootstraps = TRUE, multithread=TRUE)
write.csv(taxa, "~/my/path/taxa_UNITE_10.0_table.csv")
