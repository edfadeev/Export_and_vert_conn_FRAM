#load libraries
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")


PS99_phy0<- readRDS("Data/PS99_phy0.rds")
enriched_deep<- read.csv("Tables/enriched_deep.txt")

#####################################
#Barplot of enriched ASV sequence proportion
#####################################
#extract enriched ASVs along the water column 
PS99_phy0.prev.ra <- transform_sample_counts(PS99_phy0, function(x) x / sum(x))

PS99_phy0_sed.ra_EGC<- subset_samples(PS99_phy0.prev.ra, Region == "EGC")
PS99_phy0_sed.ra_EGC <- prune_taxa(enriched_deep$ASV[enriched_deep$region =="EGC"],PS99_phy0_sed.ra_EGC)
PS99_phy0_sed.ra_EGC.long <- psmelt(PS99_phy0_sed.ra_EGC)

PS99_phy0_sed.ra_WSC<- subset_samples(PS99_phy0.prev.ra, Region == "WSC")
PS99_phy0_sed.ra_WSC <- prune_taxa(enriched_deep$ASV[enriched_deep$region =="WSC"],PS99_phy0_sed.ra_WSC)
PS99_phy0_sed.ra_WSC.long <- psmelt(PS99_phy0_sed.ra_WSC)

BAC_pruned.ra.long <-rbind(PS99_phy0_sed.ra_EGC.long,PS99_phy0_sed.ra_WSC.long)

PS99_phy0_sed_sub<- merge_phyloseq(PS99_phy0_sed.ra_EGC,PS99_phy0_sed.ra_WSC)

#transform data
BAC_pruned.ra.long <- psmelt(PS99_phy0_sed_sub)
BAC_pruned.ra.long$Abundance <- BAC_pruned.ra.long$Abundance*100
BAC_pruned.ra.long<- BAC_pruned.ra.long[BAC_pruned.ra.long$Abundance>0,]

#calculate abundance for each taxa
BAC_pruned.ra.long.agg <- aggregate(Abundance~StationName+Fraction+Type+Class, BAC_pruned.ra.long, FUN = "sum")



#order of stations
levels(BAC_pruned.ra.long.agg$StationName) <- c("EG1","EG4","N5","N4","HG9", "HG4","HG2","HG1","S3")
BAC_pruned.ra.long.agg$Class <- as.character(BAC_pruned.ra.long.agg$Class)

#remove below 2% ra
taxa_classes <- unique(BAC_pruned.ra.long.agg$Class[!BAC_pruned.ra.long.agg$Abundance<2])

BAC_pruned.ra.long.agg$Class[BAC_pruned.ra.long.agg$Abundance<2] <- "Other taxa"

BAC_pruned.ra.long.agg$Class <- factor(BAC_pruned.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))

#Plot 
barplots_sediment <- ggplot(BAC_pruned.ra.long.agg, 
                            aes(x = StationName, y = Abundance,
                                fill = Class)) + 
  facet_grid(Type~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = phyla.col )+ 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  ylim(0,100)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle =90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


ggsave("./figures/sediment.png", 
       plot = barplots_sediment,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)




#####################################
#Explore results
#####################################
BAC_pruned.ra.long.tax <- aggregate(Abundance~Region+StationName+Fraction+Type+Class+Order+Family+Genus, BAC_pruned.ra.long, FUN = "sum")


BAC_pruned.ra.long.tax<-as.data.frame(as.list(aggregate(Abundance~Region+StationName+Fraction+Type+Class+Order+Family,
                                                        BAC_pruned.ra.long, 
                                                        FUN = function(x) c(sum = sum(x), count=length(x)))))


BAC_pruned.ra.long.samp <- aggregate(Abundance~Region+StationName+Fraction+Type, BAC_pruned.ra.long, FUN = "sum")


#summary of sequence proportion of the shared ASVs
BAC_pruned.ra.long.summary <- as.data.frame(as.list(aggregate(Abundance~Region+Fraction+Type, BAC_pruned.ra.long.samp, FUN = function(x) c(mean = mean(x), sd = sd(x)))))


sediment_ASVs<- taxa_names(PS99_phy0_sed)

length(intersect(sediment_ASVs, enriched_deep$ASV[enriched_deep$region =="EGC"]))
length(intersect(sediment_ASVs, enriched_deep$ASV[enriched_deep$region =="WSC"]))

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))

#unload libraries
detach("package:phyloseq")
detach("package:dplyr")
detach("package:DESeq2")
detach("package:ggplot2")
detach("package:reshape2")