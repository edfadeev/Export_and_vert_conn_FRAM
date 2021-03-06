#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files("Clipped", pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files("Clipped", pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_clip_R1.fastq"), `[`, 1)

filt_path <- file.path("Report")
if(!file_test("-d", filt_path)) dir.create(filt_path)

# quality check
QualityProfileFs <- list()
for(i in 1:length(fnFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs[i])
}
pdf(file.path("Report","RawProfileForward.pdf"))
for(i in 1:length(fnFs)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(fnRs)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(fnRs[i])
}
pdf(file.path("Report","RawProfileReverse.pdf"))
for(i in 1:length(fnRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Make directory and filenames for the filtered fastqs
filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#separate the different runs
#MiSeq 1 - 493:291
fnFs_mi1 <- sort(file.path("Clipped",paste(c(9:13,15:21,23:29,30:35,37:47,49,53:59,61:65), "_clip_R1.fastq", sep = "")))
fnRs_mi1 <- sort(file.path("Clipped",paste(c(9:13,15:21,23:29,30:35,37:47,49,53:59,61:65), "_clip_R2.fastq", sep = "")))
filtFs_mi1 <- sort(file.path("Filtered",paste(c(9:13,15:21,23:29,30:35,37:47,49,53:59,61:65), "_F_filt.fastq.gz", sep = "")))
filtRs_mi1 <- sort(file.path("Filtered",paste(c(9:13,15:21,23:29,30:35,37:47,49,53:59,61:65), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi1 <- filterAndTrim(fnFs_mi1, filtFs_mi1, fnRs_mi1, filtRs_mi1, truncLen=c(230,195),
                         maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi1 <- learnErrors(filtFs_mi1, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi1 <- learnErrors(filtRs_mi1, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi1 <- dada(filtFs_mi1, err=errF_mi1, multithread=TRUE, verbose = TRUE)
dadaRs_mi1 <- dada(filtRs_mi1, err=errR_mi1, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi1 <- mergePairs(dadaFs_mi1, filtFs_mi1, dadaRs_mi1, filtRs_mi1, verbose=TRUE, minOverlap = 10)

#MiSeq 2 - 493:295
fnFs_mi2 <- sort(file.path("Clipped",paste(c(14,48,50:52,60), "_clip_R1.fastq", sep = "")))
fnRs_mi2 <- sort(file.path("Clipped",paste(c(14,48,50:52,60), "_clip_R2.fastq", sep = "")))
filtFs_mi2 <- sort(file.path("Filtered",paste(c(14,48,50:52,60), "_F_filt.fastq.gz", sep = "")))
filtRs_mi2 <- sort(file.path("Filtered",paste(c(14,48,50:52,60), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi2 <- filterAndTrim(fnFs_mi2, filtFs_mi2, fnRs_mi2, filtRs_mi2, truncLen=c(230,195),
                         maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi2 <- learnErrors(filtFs_mi2, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi2 <- learnErrors(filtRs_mi2, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi2 <- dada(filtFs_mi2, err=errF_mi2, multithread=TRUE, verbose = TRUE)
dadaRs_mi2 <- dada(filtRs_mi2, err=errR_mi2, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi2 <- mergePairs(dadaFs_mi2, filtFs_mi2, dadaRs_mi2, filtRs_mi2, verbose=TRUE, minOverlap = 10)

#MiSeq 3 - MISEQ:357
fnFs_mi3 <- sort(file.path("Clipped",paste(c(1:7), "_clip_R1.fastq", sep = "")))
fnRs_mi3 <- sort(file.path("Clipped",paste(c(1:7), "_clip_R2.fastq", sep = "")))
filtFs_mi3 <- sort(file.path("Filtered",paste(c(1:7), "_F_filt.fastq.gz", sep = "")))
filtRs_mi3 <- sort(file.path("Filtered",paste(c(1:7), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi3 <- filterAndTrim(fnFs_mi3, filtFs_mi3, fnRs_mi3, filtRs_mi3, truncLen=c(230,195),
                         maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi3 <- learnErrors(filtFs_mi3, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi3 <- learnErrors(filtRs_mi3, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi3 <- dada(filtFs_mi3, err=errF_mi3, multithread=TRUE, verbose = TRUE)
dadaRs_mi3 <- dada(filtRs_mi3, err=errR_mi3, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi3 <- mergePairs(dadaFs_mi3, filtFs_mi3, dadaRs_mi3, filtRs_mi3, verbose=TRUE, minOverlap = 10)

#MiSeq 4 - 493:292
fnFs_mi4 <- sort(file.path("Clipped",paste(c(8,22,36), "_clip_R1.fastq", sep = "")))
fnRs_mi4 <- sort(file.path("Clipped",paste(c(8,22,36), "_clip_R2.fastq", sep = "")))
filtFs_mi4 <- sort(file.path("Filtered",paste(c(8,22,36), "_F_filt.fastq.gz", sep = "")))
filtRs_mi4 <- sort(file.path("Filtered",paste(c(8,22,36), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi4 <- filterAndTrim(fnFs_mi4, filtFs_mi4, fnRs_mi4, filtRs_mi4, truncLen=c(230,195),
                         maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi4 <- learnErrors(filtFs_mi4, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi4 <- learnErrors(filtRs_mi4, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi4 <- dada(filtFs_mi4, err=errF_mi4, multithread=TRUE, verbose = TRUE)
dadaRs_mi4 <- dada(filtRs_mi4, err=errR_mi4, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi4 <- mergePairs(dadaFs_mi4, filtFs_mi4, dadaRs_mi4, filtRs_mi4, verbose=TRUE, minOverlap = 10)

#MiSeq 4 - MISEQ:364
fnFs_mi5 <- sort(file.path("Clipped",paste(c(66:72), "_clip_R1.fastq", sep = "")))
fnRs_mi5 <- sort(file.path("Clipped",paste(c(66:72), "_clip_R2.fastq", sep = "")))
filtFs_mi5 <- sort(file.path("Filtered",paste(c(66:72), "_F_filt.fastq.gz", sep = "")))
filtRs_mi5 <- sort(file.path("Filtered",paste(c(66:72), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi5 <- filterAndTrim(fnFs_mi5, filtFs_mi5, fnRs_mi5, filtRs_mi5, truncLen=c(230,195),
                         maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi5 <- learnErrors(filtFs_mi5, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi5 <- learnErrors(filtRs_mi5, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi5 <- dada(filtFs_mi5, err=errF_mi5, multithread=TRUE, verbose = TRUE)
dadaRs_mi5 <- dada(filtRs_mi5, err=errR_mi5, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi5 <- mergePairs(dadaFs_mi5, filtFs_mi5, dadaRs_mi5, filtRs_mi5, verbose=TRUE, minOverlap = 10)

save.image("vert_conn_dada2_sep_runs.Rdata")

##summary of filtering
# quality check
QualityProfileFs <- list()
for(i in 1:length(filtFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(filtFs[i])
}
pdf(file.path("Report","FiltProfileForward.pdf"))
for(i in 1:length(filtFs)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(filtRs)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(filtRs[i])
}
pdf(file.path("Report","FiltProfileReverse.pdf"))
for(i in 1:length(filtRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Plot error profiles
pdf(file.path("Report","ErrorProfiles.pdf"))
plotErrors(errF_mi1, nominalQ = TRUE)+ggtitle("F-493:291")
plotErrors(errF_mi2, nominalQ = TRUE)+ggtitle("F-493:295")
plotErrors(errF_mi3, nominalQ = TRUE)+ggtitle("F-MISEQ:357")
plotErrors(errF_mi4, nominalQ = TRUE)+ggtitle("F-493:292")
plotErrors(errF_mi5, nominalQ = TRUE)+ggtitle("F-MISEQ:364")
plotErrors(errR_mi1, nominalQ = TRUE)+ggtitle("R-493:291")
plotErrors(errR_mi2, nominalQ = TRUE)+ggtitle("R-493:295")
plotErrors(errR_mi3, nominalQ = TRUE)+ggtitle("R-MISEQ:357")
plotErrors(errR_mi4, nominalQ = TRUE)+ggtitle("R-493:292")
plotErrors(errR_mi5, nominalQ = TRUE)+ggtitle("R-MISEQ:364")
dev.off()


#write out filtered read counts
write.csv(rbind(out_mi1,out_mi2,out_mi3,out_mi4,out_mi5), 
          file= file.path("Report","dada2_filterAndTrim_output.csv"))



#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(table1= makeSequenceTable(mergers_mi5), 
                             table2 = makeSequenceTable(mergers_mi1),
                             table3 = makeSequenceTable(mergers_mi2),
                             table4 = makeSequenceTable(mergers_mi3),
                             table5 = makeSequenceTable(mergers_mi4))

save.image("vert_conn_dada2_sep_runs.Rdata")

#Combine together sequences that are identical 
seqtab1 <- collapseNoMismatch(seqtab, verbose = TRUE)

dim(seqtab1)

save.image("vert_conn_dada2_sep_runs.Rdata")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))

seqtab.nochim <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#proportion of chimeras
sum(seqtab.nochim)/sum(seqtab1)

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V45 amplicons / based upon output table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(365:385) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "../tax/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE, verbose = TRUE)
taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz", tryRC = TRUE, verbose = TRUE)

save.image("vert_conn_dada2_sep_runs.Rdata")

# get summary tables 
sample.order <- names(c(dadaFs_mi1,dadaFs_mi2,dadaFs_mi3,dadaFs_mi4,dadaFs_mi5))

getN <- function(x) sum(getUniques(x))
track <- cbind(rbind(out_mi1,out_mi2,out_mi3,out_mi4,out_mi5), 
               sapply(c(dadaFs_mi1,dadaFs_mi2,dadaFs_mi3,dadaFs_mi4,dadaFs_mi5), getN),
               sapply(c(mergers_mi1,mergers_mi2,mergers_mi3,mergers_mi4,mergers_mi5), getN),
               rowSums(seqtab.nochim[sample.order,]),
               rowSums(seqtab.nochim2[sample.order,]))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nochim", "tabled")
rownames(track) <- gsub("_F_filt.fastq.gz","",sample.order)
track <- data.frame(track)

#add unclassified levels of taxonomy 
TAX <- taxa
k <- ncol(TAX) - 1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) > 1) {
    test <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(test)) {
      if (sum(is.na(test[j, i:(k + 1)])) == length(test[j, i:(k + 1)])) {
        test[j, i] <- paste(test[j, (i - 1)], "_uncl", sep = "")
        test[j, (i + 1):(k + 1)] <- test[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- test
  }
  if (sum(is.na(TAX[, i])) == 1) {
    test <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(test[i:(k + 1)])) == length(test[i:(k + 1)])) {
      test[i] <- paste(test[(i - 1)], "_uncl", sep = "")
      test[(i + 1):(k + 1)] <- test[i]
    }
    TAX[is.na(TAX[, i]),] <- test
  }
}
TAX[is.na(TAX[, (k + 1)]), (k + 1)] <- paste(TAX[is.na(TAX[, (k + 1)]), k], "_uncl", sep = "")

# write output
write.table(t(seqtab.nochim2), "dada2/dada2_seqtab_nochim2.txt", quote = F, sep = "\t")
write.table(TAX, "dada2/dada2_taxonomy_table.txt", sep = "\t", quote = F)
write.table(track, "dada2/libs_summary_table.txt", sep = "\t", quote = F)

save.image("vert_conn_dada2_sep_runs.Rdata")
