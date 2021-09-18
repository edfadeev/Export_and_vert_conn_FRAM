#######################################
## DADA2 ANALYSIS OF AMPLICON DATA ##
#######################################
require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)

# list files
path <- "/scratch2/efadeev/PS99/dada2/Clipped/"
fns <- list.files(path)
fns

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

# Extract sample names - provide table with one sample-ID per line
# assuming filenames have format: SAMPLENAME_XXX.fastq (e.g. 1_clip as from swarm-clipping output)
sample.names <- sort(read.table("./sample_names.txt", h = F, stringsAsFactors = F)$V1)
sample.names <- gsub("_clip_R.*","",sample.names)
sample.names <- unique(paste("X",sample.names, sep =""))

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# quality check
QualityProfileFs <- list()
for(i in 1:length(fnFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs[i])
}
pdf("QualityProfileForward.pdf")
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
pdf("QualityProfileReverse.pdf")
for(i in 1:length(fnRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)
# expected max length: 400 // min overlap: 30

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "../Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Filter: depending on the expected overlap, filtering parameters can be stricter
# calculate expected error rates: # sum(10^(-x/10))  #  x: quality scores of each read
# truncLen was adjusted according to QualityProfile output (low quality esp of rev-reads)
out <- filterAndTrim(
  fnFs, 
  filtFs, 
  fnRs, 
  filtRs,
  truncLen = c(230, 195),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 10
)
head(out)
summary(out[, 2]/out[, 1])
# should be retaining >70% (0.7) OK here!

# Quality check 
QualityProfileFs.filt <- list()
for(i in 1:length(filtFs)) {
  QualityProfileFs.filt[[i]] <- list()
  QualityProfileFs.filt[[i]][[1]] <- plotQualityProfile(filtFs[i])
}
pdf("QualityProfileForwardFiltered.pdf")
for(i in 1:length(filtFs)) {
  do.call("grid.arrange", QualityProfileFs.filt[[i]])  
}
dev.off()
rm(QualityProfileFs.filt)

QualityProfileRs.filt <- list()
for(i in 1:length(filtRs)) {
  QualityProfileRs.filt[[i]] <- list()
  QualityProfileRs.filt[[i]][[1]] <- plotQualityProfile(filtRs[i])
}
pdf("QualityProfileReverseFiltered.pdf")
for(i in 1:length(filtRs)) {
  do.call("grid.arrange", QualityProfileRs.filt[[i]])  
}
dev.off()
rm(QualityProfileRs.filt)

# Learn errors 
errF <- learnErrors(filtFs, multithread = 10, randomize = TRUE, verbose = 1, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread = 10, randomize = TRUE, verbose = 1, MAX_CONSIST = 20)

# Plot error profiles
pdf("ErrorProfiles.pdf")
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()
# should not be too many outliers outside the black line - but still ok here

# Dereplication 
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# denoising
dadaFs <- dada(derepFs, err = errF, multithread = 10, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = 10, pool = TRUE)

# Read merging
mergers <- mergePairs(
  dadaFs, 
  derepFs, 
  dadaRs,
  derepRs,
  minOverlap = 10,
  verbose = TRUE,
  propagateCol = c("birth_fold", "birth_ham")
)

# create sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #identified xxx sequences

# removing chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 10, verbose = TRUE)
dim(seqtab.nochim) 
summary(rowSums(seqtab.nochim)/rowSums(seqtab))
table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim)))

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V4-5 amplicons / based upon oputput table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(365:385) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

# get summary overview 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim2))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nochim", "tabled")
rownames(track) <- sample.names
track <- data.frame(track)
head(track)
summary(track$tabled/track$input) 
summary(track$filtered/track$input)
summary(track$denoised/track$filtered)
summary(track$merged/track$denoised)
summary(track$nochim/track$merged)
summary(track$tabled/track$nochim)
summary(rowSums(seqtab.nochim2))

# Assign taxonomy
tax <- assignTaxonomy(
  seqtab.nochim2, 
  "silva_nr_v132_train_set.fa.gz", 
  tryRC = TRUE,
  multithread = 10)
summary(tax)

# select BACTERIA AND ARCHAEA and remove OTUs unclassified on phylum level 
# also removes Mitochondria and Chloroplast seqs, if desired
table(tax[, 1])   
sum(is.na(tax[, 2]))  
tmp.arch <- tax[!is.na(tax[, 2]) & tax[, 1] %in% c("Bacteria", "Archaea"),]
tax.good.arch <- tmp.arch[-c(grep("Chloroplast", tmp.arch[, 4]), grep("Mitochondria", tmp.arch[, 5])), ]
seqtab.nochim2.good.arch <- seqtab.nochim2[, rownames(tax.good.arch)]
summary(rowSums(seqtab.nochim2.good.arch))

# format output
seqtab.nochim2.print.arch <- t(seqtab.nochim2.good.arch)
tax.print.arch <- tax.good.arch
all.equal(rownames(seqtab.nochim2.print.arch), rownames(tax.print.arch))  #TRUE
rownames(seqtab.nochim2.print.arch) <- paste("sq", 1:ncol(seqtab.nochim2.good.arch), sep = "")
rownames(tax.print.arch) <- rownames(seqtab.nochim2.print.arch)

# write output
write.table(seqtab.nochim2.print.arch, "./Data/dada2_output/bac-arch_seqtab_nochim2.txt", quote = F, sep = "\t")
write.table(tax.print.arch, "./Data/dada2_output/bac-arch_taxonomy_table.txt", sep = "\t", quote = F)
