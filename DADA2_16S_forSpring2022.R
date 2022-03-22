#First we load the necessary libraries. 
#The ShortRead and ggplot2 packages are available from Bioconductor
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

#Tell me where is the working directory (this command its like pwd in UNIX)
getwd()

#Change the working directory (like cd)
setwd("~/data/MiSeq_SOP/")

#double check where are you working
getwd()

#Define the following path variable so that it points to the extracted directory on your machine:
#create a variable called path --> so this will be your working directory for the further script, when the word path occurs in a command. so you have to define where is your input and where has to be made the output)
path <- getwd() #CHANGE ME to the directory containing the fastq files after unzipping. (path = getwd)
path

#List the files in your current/working directory
fns <- list.files(path)
fns

# Forward and reverse fastq filenames
#If you files are a different format you should change the commands below
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Visualize the quality profile of the forward reads:
plot1 <- plotQualityProfile(fnFs[[1]])
plot2 <- plotQualityProfile(fnFs[[2]])
#Visualize the quality profile of the reverse reads:
plot3 <- plotQualityProfile(fnRs[[1]])
plot4 <- plotQualityProfile(fnRs[[2]])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(6,6), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#The standard filtering parameters are starting points, not set in stone. 
#If you want to speed up downstream computation, consider tightening maxEE. 
#If too few reads are passing the filter, consider relaxing maxEE, 
#perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing 
#the truncLen to remove low quality tails. #The maxEE parameter sets the maximum number 
#of expected errors allowed in a read, 
#Remember though, when choosing #truncLen for paired-end reads you must maintain overlap 
#after truncation in order to merge them later.

# Save this output as RDS file:
saveRDS(out, "filter_and_trim_out.rds")

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plot_errF <- plotErrors(errF, nominalQ=TRUE)

# save error calculation as RDS files:
saveRDS(errF, "errF.rds")
saveRDS(errR, "errR.rds")

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# save dereplication as RDS files:
saveRDS(derepFs, "derepF.rds")
saveRDS(derepRs, "derepR.rds")

#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)
#Inspecting the returned dada-class object:
dadaFs[[1]]

# Save sequence-variant inference as RDS files which may be uploaded in case R crashes (or save session, see above): 
saveRDS(dadaFs, "dadaFs.rds")
saveRDS(dadaRs, "dadaRs.rds")

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#save mergers
saveRDS(mergers,"mergers.rds")


#Constructing the sequence table
#We can now construct a sequence table of the samples that is analogous to the OTU table
#produced by classical methods.
#Construct sequence table:
seqtab <- makeSequenceTable(mergers)
#Check the dimensions of the sequence table
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#The sequence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, 
#and may be worth removing, e.g.
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)])
#This is analogous to cutting a band in-silico to get amplicons of the targeted length.

# Save sequence table
saveRDS(seqtab, "seqtab.rds")

#Remove chimeras
#The core dada method removes substitution and indel errors, 
#but chimeras remain. Fortunately, the accuracy of the sequences after denoising makes 
#identifying chimeras easier than it is when dealing with fuzzy OTUs: 
#all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequences.
#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. Here chimeras make up about XX% of the inferred sequence variants, 
#but those variants account for only about XX% of the total sequence reads.
#Most of your reads should remain after chimera removal 
#(it is not uncommon for a majority of sequence variants to be removed though). 
#If most of your reads were removed as chimeric, upstream processing may need to be revisited. 
#In almost all cases this is caused by primer sequences with ambiguous nucleotides that 
#were not removed prior to beginning the DADA2 pipeline.

# Save table with the non-chimeric sequences as rds-file:
saveRDS(seqtab.nochim, "seqtab_nochim.rds")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#This is a great place to do a last sanity check. 
#Outside of filtering (depending on how stringent you want to be) 
#there should no step in which a majority of reads are lost. 
#If a majority of reads failed to merge, you may need to revisit the  
#truncLen parameter used in the filtering step and make sure that the 
#truncated reads span your amplicon. If a majority of reads were removed 
#as chimeric, you may need to revisit the removal of primers, as the ambiguous 
#nucleotides in unremoved primers interfere with chimera identification.
write.table(track,"track.txt",sep="\t",col.names = NA)

#Assign taxonomy
#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, 
#to classify sequence variants taxonomically. The DADA2 package provides a native implementation 
#of the RDPs naive Bayesian classifier for this purpose. 
#The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, 
#and outputs the taxonomic assignments with at least minBoot bootstrap confidence.
#Assign taxonomy:

taxa_silva <- assignTaxonomy(seqtab.nochim, "~/data/taxonomy/silva_nr99_v138.1_train_set.fa.gz")
taxa_silva_revcomp <- assignTaxonomy(seqtab.nochim, "~/data/taxonomy/silva_nr99_v138.1_train_set.fa.gz", tryRC=TRUE)

taxa_GTDB <- assignTaxonomy(seqtab.nochim, "~/data/taxonomy/GTDB_bac-arc_ssu_r86.fa.gz")
taxa_GTDB_revcomp <- assignTaxonomy(seqtab.nochim, "~/data/taxonomy/GTDB_bac-arc_ssu_r86.fa.gz", tryRC=TRUE)

#taxa_silva
#colnames(taxa_silva)
#taxa_silva_revcomp
#colnames(taxa_silva_revcomp)

taxa_silva <- addSpecies(taxa_silva, "~/data/taxonomy/silva_species_assignment_v138.1.fa.gz")
taxa_silva_revcomp <- addSpecies(taxa_silva_revcomp, "~/data/taxonomy/silva_species_assignment_v138.1.fa.gz")

taxa_GTDB <- addSpecies(taxa_GTDB, "~/data/taxonomy/GTDB_dada2_assignment_species.fa.gz")
taxa_GTDB_revcomp <- addSpecies(taxa_GTDB_revcomp, "~/data/taxonomy/GTDB_dada2_assignment_species.fa.gz")

colnames(taxa_silva) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
unname(head(taxa_silva))
#head(taxa_silva)
colnames(taxa_silva_revcomp) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
unname(head(taxa_silva_revcomp))
#head(taxa_silva_revcomp)
colnames(taxa_GTDB) <- c("domain", "phylum", "class", "order", "family", "genus", "species", "subspecies")
unname(head(taxa_GTDB))
colnames(taxa_GTDB_revcomp) <- c("domain", "phylum", "class", "order", "family", "genus", "species", "subspecies")
unname(head(taxa_GTDB_revcomp))


#Save sequence file:
sink("seqs.fa");cat(paste(">","ASV_",seq(1:dim(seqtab.nochim)[2]),"\n",paste(colnames(seqtab.nochim),"\n",sep=""),sep=""),sep="");sink()
seqtab.final<-seqtab.nochim
colnames(seqtab.final)<-paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep="")
seqtab.final
#Save the sequence table:
write.csv(t(seqtab.final),"seq_table.csv",quote=F)

#Save the sequence taxonomy:
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_silva)),"seq_Taxonomy_silva.csv",sep=",",col.names=F,quote=F,na="")
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_silva_revcomp)),"seq_Taxonomy_silva_revcomp.csv",sep=",",col.names=F,quote=F,na="")
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_GTDB)),"seq_Taxonomy_GTDB.csv",sep=",",col.names=F,quote=F,na="")
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_GTDB_revcomp)),"seq_Taxonomy_GTDB_revcomp.csv",sep=",",col.names=F,quote=F,na="")

#Come out of R and check the files. DADA2 is done. 