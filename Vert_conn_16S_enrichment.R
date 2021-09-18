#load libraries
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")


#load dataset
PS99_phy0_water<- readRDS("Data/PS99_water_prev.rds")

source("scripts/col_palette.R")
source("scripts/functions.R")


#####################################
# Enrichment test using DESeq2
####################################
#subset communities
PS99_PA <- subset_samples(PS99_phy0_water, Fraction == "PA")
PS99_PA <- prune_taxa(taxa_sums(PS99_PA)>0,PS99_PA)

#DCM EPI
PS99_DCM_EPI <- subset_samples(PS99_PA, Type %in% c("SRF","EPI"))
PS99_DCM_EPI <- prune_taxa(taxa_sums(PS99_DCM_EPI)>0,PS99_DCM_EPI)

#EPI MESO
PS99_EPI_MESO <- subset_samples(PS99_PA, Type %in% c("EPI","MESO"))
PS99_EPI_MESO <- prune_taxa(taxa_sums(PS99_EPI_MESO)>0,PS99_EPI_MESO)

#MESO BATHY
PS99_MESO_BATHY <- subset_samples(PS99_PA, Type %in% c("MESO","BATHY"))
PS99_MESO_BATHY <- prune_taxa(taxa_sums(PS99_MESO_BATHY)>0,PS99_MESO_BATHY)

#run DEseq2
type <- c("PS99_DCM_EPI","PS99_EPI_MESO","PS99_MESO_BATHY")
deseq_res_all <- data.frame()
enriched_agg_all <- data.frame()

for (n in c("EGC","WSC")){
  for (i in 1:3){
    #run DEseq
    sub <- subset_samples(get(type[i]), Region == n)
    BAC_sed.ddsMat <- phyloseq_to_deseq2(sub, ~Type)
    geoMeans = apply(counts(BAC_sed.ddsMat), 1, gm_mean)
    BAC_sed.ddsMat = estimateSizeFactors(BAC_sed.ddsMat, geoMeans = geoMeans)
    BAC_sed.ddsMat <- estimateDispersions(BAC_sed.ddsMat)
    BAC_sed.DEseq <- DESeq(BAC_sed.ddsMat, fitType="local")
    BAC_sed.DEseq.res <- results(BAC_sed.DEseq)
    
    #extract only significant OTU
    BAC_sed.DEseq.res.sig <- BAC_sed.DEseq.res[which(BAC_sed.DEseq.res$padj < 0.1), ]
    BAC_sed.DEseq.res.sig <- cbind(as(BAC_sed.DEseq.res.sig, "data.frame"),
                                   as(tax_table(sub)[rownames(BAC_sed.DEseq.res.sig), ], "matrix"))
    BAC_sed.DEseq.res.sig.sub <- BAC_sed.DEseq.res.sig[abs(BAC_sed.DEseq.res.sig$log2FoldChange)>1,]
    BAC_sed.DEseq.res.sig.sub$merging <- type[i]
    BAC_sed.DEseq.res.sig.sub$region <- n
    BAC_sed.DEseq.res.sig.sub$ASV <- rownames(BAC_sed.DEseq.res.sig.sub)
    rownames(BAC_sed.DEseq.res.sig.sub) <- c()
    deseq_res_all <- rbind(deseq_res_all,BAC_sed.DEseq.res.sig.sub, make.row.names = TRUE)
  }
}

#####################################
#Explore the results
#####################################
#separate the results by enrichment towards surface and to depth
enriched_shallow <- deseq_res_all[deseq_res_all[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","ASV","merging","region") ]
enriched_deep  <- deseq_res_all[deseq_res_all[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","ASV","merging","region")]

#export the table of enriched with depth for Fig 6
write.csv(enriched_deep, "Tables/enriched_deep.txt")

#Aggregate on Order level 
enriched_shallow.agg <-as.data.frame(as.list(aggregate(log2FoldChange~region+merging+Phylum+Class+Order+Family,
                                                       enriched_shallow, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

enriched_deep.agg <-as.data.frame(as.list(aggregate(log2FoldChange~region+merging+Phylum+Class+Order+Family,
                                                    enriched_deep, 
                                                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#merge data for ploting
enriched_agg <- rbind(enriched_shallow.agg,enriched_deep.agg)

#exclude taxa with less than 3 ASV
enriched_agg_top <- enriched_agg[enriched_agg$log2FoldChange.count>3,]

#order the x axis by classes
enriched_agg_top$Class <- factor(enriched_agg_top$Class, ordered = TRUE,
                                 levels= sort(unique(as.character(enriched_agg_top$Class)),decreasing=FALSE))
enriched_agg_top$Family <- factor(enriched_agg_top$Family, ordered = TRUE,
                                  levels= unique(enriched_agg_top$Family[order(enriched_agg_top$Class)]))
enriched_agg_top$merging <- factor(enriched_agg_top$merging,
                                   levels= c("PS99_DCM_EPI","PS99_EPI_MESO","PS99_MESO_BATHY"))
enriched_agg_top$region <- factor(enriched_agg_top$region,
                                  levels= c("EGC","WSC"))


#plot
PS99_daOTU.p <- ggplot(data=enriched_agg_top, aes(y=log2FoldChange.mean , x=Family, fill = Class, label = log2FoldChange.count))+ 
  #geom_point(data=enriched_agg_top, aes(y=log2FoldChange.mean , x=Order, fill = Class, label = log2FoldChange.count), size = 0, shape = 21)+
  geom_text(data=enriched_agg_top,aes(y=log2FoldChange.mean , x=Family), nudge_y= -1.5, nudge_x= 0)+
  geom_errorbar(data=enriched_agg_top,aes(ymin = log2FoldChange.mean-log2FoldChange.se, ymax = log2FoldChange.mean +log2FoldChange.se), width = 0.2) +   
  ylab("log2foldchange")+
  geom_point(data=enriched_agg_top[enriched_agg_top$log2FoldChange.mean<0,], aes(y=log2FoldChange.mean , x=Family, fill = Class, label = log2FoldChange.count), size = 5, shape = 24)+
  geom_point(data=enriched_agg_top[enriched_agg_top$log2FoldChange.mean>0,], aes(y=log2FoldChange.mean , x=Family, fill = Class, label = log2FoldChange.count), size = 5, shape = 25)+
  scale_y_reverse()+
  #scale_x_discrete("Order")+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_manual(values = phyla.col)+
  guides(shape = 22)+
  facet_grid(merging~region)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle =90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


ggsave("./figures/enriched_ASVs.png", 
       plot = PS99_daOTU.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Explore results
#####################################
#overall enriched ASVs in Families by region
dASVs.summary_by_region <-as.data.frame(as.list(aggregate(ASV~region+Class+Order+Family, 
                                                          enriched_deep, 
                                                          FUN = function(x) c(count=length(x)))))

#overall ASVs in Families
dASVs.summary_total.phylum <-as.data.frame(as.list(aggregate(ASV~region+Phylum+Class, 
                                                             unique(enriched_deep[,c("region", "Phylum","Class","Order","Family","Genus","ASV")]), 
                                                             FUN = function(x) c(count=length(x)))))



dASVs.summary_total <-as.data.frame(as.list(aggregate(ASV~Phylum+Class+Order+Family, 
                                                      unique(enriched_deep[,c("Phylum","Class","Order","Family","Genus","ASV")]), 
                                                      FUN = function(x) c(count=length(x)))))


#Generate overview table

#summarize enriched with depth
enriched_EGC.agg.depth <-as.data.frame(as.list(aggregate(ASV~Phylum+Class+Order+Family+Genus, 
                                                         enriched_deep[enriched_deep$region=="EGC",], 
                                                         FUN = function(x) c(count=length(x)))))

enriched_WSC.agg.depth <-as.data.frame(as.list(aggregate(ASV~Phylum+Class+Order+Family+Genus,
                                                         enriched_deep[enriched_deep$region=="WSC",], 
                                                         FUN = function(x) c(count=length(x)))))


#merge together
enriched_ASVs_tab <-rbind(enriched_EGC.agg.depth[,c("Phylum","Class","Order","Family","Genus")], 
                          enriched_WSC.agg.depth[,c("Phylum","Class","Order","Family","Genus")]) %>%
  distinct()%>%
  left_join(enriched_EGC.agg.depth[,c("Genus","ASV")], by = "Genus")%>%
  left_join(enriched_WSC.agg.depth[,c("Genus","ASV")], by = "Genus")%>%
  select("Phylum","Class","Order","Family","Genus","ASV.x","ASV.y")

names(enriched_ASVs_tab) <- c("Phylum","Class","Order","Family","Genus","ice-covered","ice-free")

write.table(enriched_ASVs_tab, "Tables/Enriched_ASVs.txt", sep = "\t", quote = F)



#explore the results
#EGC
BAC_pruned.ra.EGC_dASVs<- BAC_pruned.ra.long %>%
  filter(Region =="EGC",
         Abundance>0,
         OTU %in% enriched_deep$ASV[enriched_deep$region =="EGC"])

EGC_dASVs_agg <-as.data.frame(as.list(aggregate(Abundance~Region+StationName+Fraction+Type+Class+Order+Family,
                                                BAC_pruned.ra.EGC_dASVs, 
                                                FUN = function(x) c(sum = sum(x)*100, count=length(x)))))

EGC_dASVs_agg.class <-as.data.frame(as.list(aggregate(Abundance~Region+StationName+Fraction+Type+Class,
                                                      BAC_pruned.ra.EGC_dASVs, 
                                                      FUN = function(x) c(sum = sum(x)*100, count=length(x)))))

EGC_dASVs_agg.class.mean <-as.data.frame(as.list(aggregate(Abundance.sum~Region+Fraction+Type+Class,
                                                           EGC_dASVs_agg.class, 
                                                           FUN = function(x) c(mean = mean(x), se = se(x),count=length(x)))))


#WSC
BAC_pruned.ra.WSC_dASVs<- BAC_pruned.ra.long %>%
  filter(Region =="WSC",
         Abundance>0,
         OTU %in% enriched_deep$ASV[enriched_deep$region =="WSC"])

WSC_dASVs_agg <-as.data.frame(as.list(aggregate(Abundance~Region+StationName+Fraction+Type+Class+Order+Family,
                                                BAC_pruned.ra.WSC_dASVs, 
                                                FUN = function(x) c(sum = sum(x)*100, count=length(x)))))

WSC_dASVs_agg.class <-as.data.frame(as.list(aggregate(Abundance~Region+StationName+Fraction+Type+Class,
                                                      BAC_pruned.ra.WSC_dASVs, 
                                                      FUN = function(x) c(sum = sum(x)*100, count=length(x)))))

WSC_dASVs_agg.class.mean <-as.data.frame(as.list(aggregate(Abundance.sum~Region+Fraction+Type+Class,
                                                           WSC_dASVs_agg.class, 
                                                           FUN = function(x) c(mean = mean(x), se = se(x),count=length(x)))))

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

