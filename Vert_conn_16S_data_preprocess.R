#load libraries
library(phyloseq); packageVersion("phyloseq")

#####################################
# Parse for Phyloseq
#####################################
OTU<- read.csv("./dada2/dada2_seqtab_nochim2.txt", h=T, sep="\t")
#correct sample names in ASV table
names(OTU)<- gsub("_F_filt.fastq.gz","",names(OTU))

TAX<- as.matrix(read.csv("./dada2/dada2_taxonomy_table.txt", h=T,sep = "\t"))

#metadata
ENV <- read.csv("./dada2/PS99_samples_meta.txt", sep = "\t" , h = T, row.names = 1, fill = T, na.strings=c("","NA"))
#reorder according to the ASV table
ENV<-ENV[names(OTU),]

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))
all.equal(names(OTU), rownames(ENV))

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
PS99_phy <- phyloseq(OTU, TAX, meta)

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(PS99_phy))
names(dna) <- taxa_names(PS99_phy)
PS99_phy <- merge_phyloseq(PS99_phy, dna)
taxa_names(PS99_phy) <- paste0("ASV", seq(ntaxa(PS99_phy)))

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
PS99_phy0 <- subset_taxa(PS99_phy, !Kingdom %in% c("Eukaryota") & !Phylum %in% c("NA_uncl" )& !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )

#####################################
# Fix categories 
#####################################
sample_data(PS99_phy0)$Type <- factor(
  sample_data(PS99_phy0)$Type, levels = c("SRF", "EPI", "MESO", "BATHY","Sediment"))

sample_data(PS99_phy0)$Fraction<- factor(
  sample_data(PS99_phy0)$Fraction, 
  levels = c("FL", "PA"))

sample_data(PS99_phy0)$StationName<- factor(
  sample_data(PS99_phy0)$StationName, 
  levels = c("EG1","EG4","N5","N4","HG9", "HG4","HG2","HG1","S3"))

#####################################
# Prevalence filter for water communities
#####################################
PS99_phy0_water<-subset_samples(PS99_phy0, Type!="Sediment")
PS99_phy0_water<- prune_taxa(taxa_sums(PS99_phy0_water)>0,PS99_phy0_water)

prevdf <-  apply(X = otu_table(PS99_phy0_water),
                 MARGIN = ifelse(taxa_are_rows(PS99_phy0_water), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
#prevdf.tax  <-  data.frame(Prevalence = prevdf,
#                      TotalAbundance = taxa_sums(PS99_phy0_water),
#                       tax_table(PS99_phy0_water))
#summarize
#prevdf.tax.summary <- plyr::ddply(prevdf.tax, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#plot
#prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(PS99_phy0),color=Phylum)) +
# Include a guess for parameter
# geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
#scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 4% of total samples
prevalenceThreshold <- round(0.04 * nsamples(PS99_phy0_water))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
PS99_water.prev <-  prune_taxa((prevdf > prevalenceThreshold), PS99_phy0_water)

#####################################
# export R-native serialized RDS file
#####################################
saveRDS(PS99_phy0, "./Data/PS99_phy0.rds")
saveRDS(PS99_water.prev, "./Data/PS99_water_prev.rds")

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))

#unload libraries
detach("package:phyloseq")
