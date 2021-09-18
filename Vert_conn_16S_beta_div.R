#load libraries
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")
library(reshape2); packageVersion("reshape2")
library(purrr); packageVersion("purrr")
library(ggpubr); packageVersion("ggpubr")


#load dataset
PS99_phy0<- readRDS("Data/PS99_phy0.rds")

source("scripts/col_palette.R")
source("scripts/functions.R")

#####################################
# Community composition
#####################################
#transform data
BAC_pruned.ra <- transform_sample_counts(PS99_phy0, function(x) x / sum(x))
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)

#calculate abundance for each taxa
BAC_pruned.ra.long.agg <- aggregate(Abundance~StationName+Fraction+Type+Class, BAC_pruned.ra.long, FUN = "sum")
BAC_pruned.ra.long.agg$Abundance <- BAC_pruned.ra.long.agg$Abundance*100
BAC_pruned.ra.long.agg<- BAC_pruned.ra.long.agg[BAC_pruned.ra.long.agg$Abundance>0,]


#order of stations
levels(BAC_pruned.ra.long.agg$StationName) <- c("EG1","EG4","N5","N4","HG9", "HG4","HG2","HG1","S3")
BAC_pruned.ra.long.agg$Class <- as.character(BAC_pruned.ra.long.agg$Class)

#remove below 3% ra
taxa_classes <- unique(BAC_pruned.ra.long.agg$Class[!BAC_pruned.ra.long.agg$Abundance<3])

BAC_pruned.ra.long.agg$Class[BAC_pruned.ra.long.agg$Abundance<3] <- "Other taxa"

BAC_pruned.ra.long.agg$Class <- factor(BAC_pruned.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))

BAC_pruned.ra.long.agg$Type <- factor(BAC_pruned.ra.long.agg$Type,
                                       levels=c("SRF","EPI","MESO","BATHY","Sediment"))

#Plot 
barplots_water <- ggplot(BAC_pruned.ra.long.agg, 
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
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggsave("./figures/barplot_water.png", 
       plot = barplots_water,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
# Ordination plot
#####################################
PS99_water.prev<- readRDS("Data/PS99_water_prev.rds")

#stabilize the dataset using gemetric mean 
# calculate geometric means prior to estimate size factors
PS99_phy0_water.dds <- phyloseq_to_deseq2(PS99_water.prev, ~1)
geoMeans = apply(counts(PS99_phy0_water.dds), 1, gm_mean)
PS99_phy0_water.dds = estimateSizeFactors(PS99_phy0_water.dds, geoMeans = geoMeans)
PS99_phy0_water.dds <- estimateDispersions(PS99_phy0_water.dds)
otu.vst <- getVarianceStabilizedData(PS99_phy0_water.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(PS99_water.prev))

PS99_phy0_water.vst<-PS99_water.prev
otu_table(PS99_phy0_water.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)


#calculate ordination
PS99_phy0_water.ord <- ordinate(PS99_phy0_water.vst, method = "RDA", distance = "eucledian")
PS99_phy0_water.ord.df <- plot_ordination(PS99_phy0_water.vst, PS99_phy0_water.ord, axes = c(1,2,3),justDF = TRUE)

#adjust grouping for clustering
PS99_phy0_water.ord.df$new_ordination <- paste(PS99_phy0_water.ord.df$Type, PS99_phy0_water.ord.df$Fraction, sep= ".")
PS99_phy0_water.ord.df$Region<- factor(PS99_phy0_water.ord.df$Region, levels = c("EGC","WSC"))

#extract explained variance
PS99.ord.evals <- 100 * summary(PS99_phy0_water.ord)$cont$importance[2, c("PC1","PC2")]

PS99.ord.p <- ggplot(data = PS99_phy0_water.ord.df, aes(x =PC1, y =PC2, shape = Fraction, colour = Region))+
  geom_point(colour="black",size = 6)+
  geom_point(size = 5)+
  #geom_text(aes(label = Type), colour = "black", nudge_y= -1,  size=3)+
  labs(x = sprintf("PC1 [%s%%]", round(PS99.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(PS99.ord.evals[2], 2)), shape = "Fraction", color = "Origin")+
  stat_ellipse(aes(group = interaction(Type,Fraction)))+
  scale_color_manual(values = c("EGC" = "blue", "WSC"="red")) +
  #coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggsave("./figures/PCA_water.png", 
       plot = PS99.ord.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#PERMANOVA test
df <- as(sample_data(PS99_phy0_water.vst), "data.frame")
df$group<- paste(df$Region,df$Type, sep ="_")
d <- phyloseq::distance(PS99_phy0_water.vst, "euclidean")
adonis_all <- adonis2(d ~ Type*Fraction*Region , data= df, perm = 999)
adonis_all

#posthoc to check all groups of samples are different
groups <- df[["group"]]
mod <- betadisper(d, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

#PERMANOVA test
df <- as(sample_data(subset_samples(PS99_phy0_water.vst, Type  %in% c("SRF","EPI"))), "data.frame")
d <- phyloseq::distance(subset_samples(PS99_phy0_water.vst, Type  %in% c("SRF","EPI")),"euclidean")
adonis_SRF <- adonis2(d ~ Fraction*Region , data= df, perm = 999)
adonis_SRF

#PERMANOVA test
df <- as(sample_data(subset_samples(PS99_phy0_water.vst, Type  %in% c("MESO","BATHY"))), "data.frame")
d <- phyloseq::distance(subset_samples(PS99_phy0_water.vst, Type  %in% c("MESO","BATHY")),"euclidean")
adonis_deep <- adonis2(d ~ Fraction*Region , data= df, perm = 999)
adonis_deep

#####################################
#plot distances distribution between fractions
#####################################
sample_data(PS99_phy0_water.vst)$SampleID <- sample_names(PS99_phy0_water.vst)
frac_distances_all <- data.frame()

frac_distances = phyloseq::distance(PS99_phy0_water.vst, "euclidean")
frac_distances.m = melt(as.matrix(frac_distances))

# remove self-comparisons
frac_distances.m = frac_distances.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
samples_data = as(sample_data(PS99_phy0_water.vst),"data.frame") %>%
  select(SampleID, Fraction, Type,Region)

# combined distances with sample data
colnames(samples_data) = c("Var1", "Fraction1", "Type1","Region1")
frac_distances.m = left_join(frac_distances.m, samples_data, by = "Var1")

colnames(samples_data) = c("Var2", "Fraction2","Type2","Region2")
frac_distances.m = left_join(frac_distances.m, samples_data, by = "Var2")


#remove duplicates
frac_distances.m <- frac_distances.m%>%
  mutate(Var = map2_chr(Var1, Var2, ~toString(sort(c(.x, .y))))) %>%
  distinct(Var, .keep_all = TRUE) %>%
  select(-Var)

#compare distances in eacg fraction along the water column
frac_distances_all_depth<- frac_distances.m%>%
  filter(Type1 == Type2,
         Region1==Region2,
         Fraction1 == "FL",
         Fraction2 == "PA")

#define region and depth category 
frac_distances_all_depth$group = factor(interaction(frac_distances_all_depth$Type1,
                                                    frac_distances_all_depth$Region1),
                                           levels= rev(c("SRF.EGC","SRF.WSC",
                                           "EPI.EGC","EPI.WSC",
                                           "MESO.EGC","MESO.WSC",
                                           "BATHY.EGC","BATHY.WSC")))
                                    
#test for normal distribution
shapiro.test(frac_distances_all_depth$value)

kruskal.test(value ~ group, data = frac_distances_all_depth)

Region_dist_Wilcox <- frac_distances_all_depth   %>%
  rstatix::wilcox_test(value ~ group, p.adjust.method = "BH") %>%
  rstatix::add_significance() 

#define the compared groups
type_comparisons <- list(c("SRF.EGC","SRF.WSC"),
                         c("EPI.EGC","EPI.WSC"),
                         c("MESO.EGC","MESO.WSC"),
                         c("BATHY.EGC","BATHY.WSC"))

#plot
dist_reg.p<- ggplot(frac_distances_all_depth, aes(x = group, y = value, fill = Region1)) +
  labs(x = "Water layer")+
  geom_boxplot(outlier.color = NULL, notch = FALSE, color = "black", size = 2)+
  geom_signif(comparisons = type_comparisons, map_signif_level=TRUE, test = "wilcox.test")+
  theme_classic() +
  scale_fill_manual(values = c("EGC" = "blue", "WSC"="red")) +
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggsave("./figures/dist_frac.pdf", 
       plot = dist_reg.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#compare distances in eacg fraction along the water column
frac_distances_all_type<- frac_distances.m%>%
  filter(Type1 != Type2,
         Region1==Region2,
         Fraction1 == "FL",
         Fraction2 == "PA")

#define region and depth category 
frac_distances_all_type$group = factor(interaction(frac_distances_all_type$Type1,
                                                   frac_distances_all_type$Region1),
                                        levels= rev(c("SRF.EGC","SRF.WSC",
                                                      "EPI.EGC","EPI.WSC",
                                                      "MESO.EGC","MESO.WSC",
                                                      "BATHY.EGC","BATHY.WSC")))

depth_comparisons <- list(c("SRF.EGC","EPI.EGC"),
                          c("EPI.EGC","MESO.EGC"), c("MESO.EGC","BATHY.EGC"),
                          c("SRF.WSC","EPI.WSC"),
                          c("EPI.WSC","MESO.WSC"), c("MESO.WSC","BATHY.WSC"))

frac_distances_all_type_sub<- frac_distances_all_type %>% filter(Region1=="EGC")

depth_comparisons <- list(c("SRF.EGC","EPI.EGC"),
                          c("EPI.EGC","MESO.EGC"), 
                          c("MESO.EGC","BATHY.EGC"))

frac_distances_all_type_sub<- frac_distances_all_type %>% 
  filter(Region1=="WSC")

depth_comparisons <- list(c("SRF.WSC","EPI.WSC"),
                          c("EPI.WSC","MESO.WSC"), 
                          c("MESO.WSC","BATHY.WSC"))

#plot
dist_frac.p<- ggplot(frac_distances_all_type_sub, aes(x = group, y = value,  fill= Region1)) +
  labs(x = "Fraction")+
  geom_boxplot(outlier.color = NULL, notch = FALSE, color = "black", size = 2)+
  geom_signif(comparisons = depth_comparisons, map_signif_level=TRUE, test = "wilcox.test")+
  scale_fill_manual(values = c("EGC" = "blue", "WSC"="red")) +
  #facet_grid(.~Region1)+
  coord_flip()+
  theme_classic() +
  theme(legend.position = "none")

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
detach("package:vegan")
detach("package:reshape2")
detach("package:purrr")
detach("package:ggpubr")
