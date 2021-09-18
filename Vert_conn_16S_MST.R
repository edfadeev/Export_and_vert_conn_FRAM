#load libraries
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")

source("scripts/col_palette.R")
source("scripts/functions.R")

#script that will run the SourceTracker on the dataset (no need to run again, just load the results)
#source("scripts/PS99_sourcetracker.R)

#load the results
ST_results<- readRDS("Data/ST_results_5K.rds")
ST_results.train <- readRDS("Data/ST_results_train_5K.rds")

####################################
# Source tracker results
#####################################
#make labels
if(is.element('StationName',colnames(metadata))) desc <- metadata$StationName
if(is.element('Region',colnames(metadata))) regi <- metadata$Region
if(is.element('StationName',colnames(metadata))) depth <- metadata$Type
if(is.element('StationName',colnames(metadata))) frac <- metadata$Fraction
if(is.element('StationName',colnames(metadata))) no <- rownames(metadata)
labels <- sprintf('%s %s %s %s %s', regi, desc, depth, frac, no)

#Mark ice-covered stations
metadata$ice <- "no"
metadata$ice[metadata$StationName == "EG1"|
               metadata$StationName == "EG4"|
               metadata$StationName == "N5"|
               metadata$StationName == "N4"] <- "yes"

#rename the the results tables
rownames(ST_results$proportions) <- labels[test.ix]
rownames(ST_results.train$proportions) <- labels[train.ix]

#melt and combined the train and test datasets
results_melted <- melt(ST_results$proportions)
results.train_melted <- melt(ST_results.train$proportions)

#merge tables and define variables order
resluts_merged <- rbind(results_melted,results.train_melted)
resluts_merged <- resluts_merged %>%
  separate(Var1, c("Region", "StationName", "Type", "Fraction", "SampleID"), " ")


#define levels
resluts_merged$Type <- factor(resluts_merged$Type, levels = c("SRF", "EPI", "MESO", "BATHY", "Sediment"))
resluts_merged$Var2<- factor(resluts_merged$Var2,levels = c("SRF", "EPI", "MESO", "BATHY","Unknown"))
resluts_merged$Fraction<- factor(resluts_merged$Fraction,levels = c("FL", "PA"))
resluts_merged$StationName<- factor(resluts_merged$StationName, levels = c("EG1","EG4","N5","N4","HG9", "HG4", "HG2","HG1","S3"))


resluts_merged_sub<-resluts_merged[resluts_merged$Type != "Sediment",]

#plot
resluts_merged_plot <- ggplot() + 
  geom_bar(aes(y = value, x = StationName, fill = Var2), colour = "black", 
           data = resluts_merged, stat="identity")+
  facet_grid(Type~Fraction)+
  labs(y="Average source contribution")+
  scale_fill_manual(values = c("SRF" = "lightskyblue1",
                               "MESO"= "blue", 
                               "EPI"="deepskyblue2",
                               "sediment"= "green",
                               "BATHY"= "darkblue", 
                               "Unknown"="gray50"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggsave("./figures/MST.pdf", 
       plot = resluts_merged_plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



#####################################
# Correlation between predicitions
#####################################
#predictions in train dataset
source_prop.train <- as.data.frame(ST_results.train$proportions)
source_prop.train$sample <- rownames(source_prop.train)
source_prop.train.trans<- source_prop.train %>%
  separate(sample, c("StationName", "Type", "Fraction", "SampleID"), " ")
rownames(source_prop.train.trans) <- source_prop.train.trans$SampleID


df.source_dist <- vegdist(source_prop.train.trans[,c("BATHY","SRF","MESO","EPI","Unknown")], "euclidean")
d.meta <- as(metadata[train.ix,],"data.frame")
adonis_source <- adonis2(df.source_dist ~ Source.Alt, d.meta)
adonis_source


#predictions in test dataset
source_prop <- as.data.frame(ST_results$proportions)
source_prop$sample <- rownames(source_prop)
source_prop.trans<- source_prop %>%
  separate(sample, c("StationName", "Type", "Fraction", "SampleID"), " ")
rownames(source_prop.trans) <- source_prop.trans$SampleID

df.sink_dist <- vegdist(source_prop.trans[,c("BATHY","EPI","MESO","SRF","Unknown")], "euclidean")
d.sink.meta <- as(metadata[test.ix,],"data.frame")
d.sink.meta<-d.sink.meta[d.sink.meta$Type!="Sediment",]

adonis_sink <- adonis2(df.sink_dist ~ Source.Alt, d.sink.meta)
adonis_sink


#calculate mean contribution
cont_test <- resluts_merged

#merge to surface and deep ocean
cont_test$Type[cont_test$Type == "EPI"]<- "SRF"
cont_test$Type[cont_test$Type == "BATHY"]<- "MESO"

cont_test$Var2[cont_test$Var2 == "EPI"]<- "SRF"
cont_test$Var2[cont_test$Var2 == "BATHY"]<- "MESO"

#sum surface and epipelgaic contribution together
cont_test_sum<- as.data.frame(as.list(aggregate(value~StationName+Type+Fraction+Var2,cont_test, FUN = sum)))
cont_test_sum_agg<- as.data.frame(as.list(aggregate(value~Type+Fraction+Var2,cont_test_sum, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#test contribution according to the different regions
cont_test_Region<- as.data.frame(as.list(aggregate(value~Region+StationName+Type+Fraction+Var2,cont_test, FUN = sum)))
test_agg_Region<- as.data.frame(as.list(aggregate(value~Region+Type+Fraction+Var2,cont_test_Region, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#####################################
# MST Summary table
#####################################
source_prop_mean <- ST_results$proportions
colnames(source_prop_mean) <- paste(colnames(source_prop_mean), "mean", sep = "_")


source_prop_RSD <- ST_results$proportions_sd
colnames(source_prop_RSD) <- paste(colnames(source_prop_RSD), "RSD", sep = "_")

#merge all into one table
as.data.frame(cbind(source_prop_mean,source_prop_RSD)) %>%
  rownames_to_column(var = "ID") %>%
  separate(ID, c("Region", "StationName","Type", "Fraction", "ID"), " ") %>%
  select(ID, StationName, Type, Fraction, SRF_mean,EPI_mean, MESO_mean, BATHY_mean,
         Unknown_mean, SRF_RSD,EPI_RSD, MESO_RSD, BATHY_RSD, Unknown_RSD) %>%
  #round numbers        
  mutate_if(is.numeric, ~round(., 2)) %>% 
  #merged mean and SD
  unite("SRF",SRF_mean,SRF_RSD, sep = "\u00B1") %>%
  unite("EPI",EPI_mean,EPI_RSD, sep = "\u00B1") %>%
  unite("MESO", MESO_mean, MESO_RSD, sep = "\u00B1") %>%
  unite("BATHY", BATHY_mean,BATHY_RSD, sep = "\u00B1")%>%
  unite("Unknown", Unknown_mean, Unknown_RSD, sep = "\u00B1") -> source_prop_table

write.table(source_prop_table,"Tables/MST_results.txt", sep = "\t", quote = F)

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))

#unload libraries
detach("package:phyloseq")
detach("package:dplyr")
detach("package:ggplot2")
detach("package:vegan")
detach("package:reshape2")