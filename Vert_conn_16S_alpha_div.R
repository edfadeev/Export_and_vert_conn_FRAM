#load libraries
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(iNEXT); packageVersion("iNEXT")

PS99_phy0<- readRDS("Data/PS99_phy0.rds")

#####################################
# Overview table
#####################################
#metadata
meta <- as(sample_data(PS99_phy0),"data.frame")
meta$Row.names<-sample_names(sample_data(PS99_phy0))

#dada2 workflow
reads.tab <- read.csv2("./dada2/libs_summary_table.txt",header = TRUE, sep = "\t")
reads.tab$Row.names <-paste("X",row.names(reads.tab), sep ="")

#alpha diversity indeces
PS99_phy0_alpha <- estimate_richness(PS99_phy0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
PS99_phy0_alpha$Row.names<-rownames(PS99_phy0_alpha)

#merge all togather
PS99_phy0_summary_table <- merge(meta,reads.tab,by ="Row.names") %>%
  merge(PS99_phy0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(tabled/input,2)) %>%
  select("StationID","StationName","Region","Type", "Longitude..degrees_east.","Latitude..degrees_north.","Depth","Fraction", "input","merged", "tabled","Seq.prop", "Observed","Chao1","Shannon","InvSimpson")%>%
  rename("Event ID" = "StationID",
         "Longitude [degE]"="Latitude..degrees_north.",
         "Latitude [degN]"="Longitude..degrees_east.",
         "Sampling depth [m]" = "Depth",
         "Raw seq." = "input",
         "Final seq." = "tabled",
         "Seq. proportions" = "Seq.prop",
         "Observed ASVs" = "Observed",
         "Chao1 richness est." = "Chao1",
         "Shannon Index" = "Shannon")


write.table(PS99_phy0_summary_table, "./Tables/Micro_overview_table.txt" , sep = "\t", quote = F)


#####################################
# Rarefaction curves
#####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(PS99_phy0)), q=0, datatype="abundance")
rare <-fortify(iNEXT.out, type=1)

meta <- as(sample_data(PS99_phy0), "data.frame")
meta$site <- rownames(meta)
rare$Fraction <- as.character(meta$Fraction[match(rare$site, meta$site)])
rare$StationName <- meta$StationName[match(rare$site, meta$site)] 
rare$Type <- meta$Type[match(rare$site, meta$site)]
rare$label <- paste(rare$StationName,rare$Type, sep = "-")

rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))


rare.line$Fraction[rare.line$Type=="Sediment"] <- "Sediment"
rare.point$Fraction[rare.point$Type=="Sediment"] <- "Sediment"

rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Type, fill=Fraction), size =3, data= rare.point)+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  facet_wrap(~Fraction, scales = "free", ncol  = 2)+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  coord_cartesian(clip="off")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggsave("./figures/rarefactions.png", 
       plot = rare.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#summarize communities in the water column
#####################################
PS99_phy0_water<-subset_samples(PS99_phy0, Type!="Sediment")
PS99_phy0_water<- prune_taxa(taxa_sums(PS99_phy0_water)>0,PS99_phy0_water)

prevdf.water <-  apply(X = otu_table(PS99_phy0_water),
                       MARGIN = ifelse(taxa_are_rows(PS99_phy0_water), yes = 1, no = 2),
                       FUN = function(x){sum(x > 0)})

# # Add taxonomy and total read counts to this data.frame
prevdf.tax.water  <-  data.frame(Prevalence = prevdf.water,
                                 TotalAbundance = taxa_sums(PS99_phy0_water),
                                 tax_table(PS99_phy0_water))
tax.summary <- plyr::ddply(prevdf.tax.water, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$TotalAbundance), length(df1$Prevalence))})
names(tax.summary)<- c("Class","No.samples","Total.reads","No.ASVs")

order.summary <- plyr::ddply(prevdf.tax.water, "Order", function(df1){cbind(mean(df1$Prevalence),sum(df1$TotalAbundance), length(df1$Prevalence))})
names(order.summary)<- c("Order","No.samples","Total.reads","No.ASVs")


#####################################
# Alpha diversity statistics
#####################################
# Create new ps object with diversity estimates added to sample_data
PS99_ps0_div <- merge_phyloseq(PS99_phy0, sample_data(PS99_phy0_alpha))

#exclude sediment samples
PS99_ps0_div <- subset_samples(PS99_ps0_div, Type != "Sediment")
PS99_ps0_div <- prune_taxa(taxa_sums(PS99_ps0_div)>0, PS99_ps0_div)


shapiro.test(sample_data(PS99_ps0_div)$Chao1)
#Chao1 richness did not show normal distribution (p < 0.01), thus will be analyzed using Kruskal Wallis test
kruskal.test(Chao1 ~ Type, data = data.frame(sample_data(PS99_ps0_div)))


shapiro.test(sample_data(PS99_ps0_div)$Shannon)
#Shanonn richness did not show normal distribution (p < 0.01), thus will be analyzed using Kruskal Wallis test
kruskal.test(Shannon ~ Type, data = data.frame(sample_data(PS99_ps0_div)))


# separate by fractions fractions
PS99_ps0_div_FL <- subset_samples(PS99_ps0_div, Fraction == "FL")
PS99_ps0_div_PA <- subset_samples(PS99_ps0_div, Fraction == "PA")


#FL
#Chao1 richness
kruskal.test(Chao1 ~ Type, data = data.frame(sample_data(PS99_ps0_div_FL)))

Chao1_Wilcox_FL <- as(sample_data(PS99_ps0_div_FL),"data.frame")   %>%
  rstatix::wilcox_test(Chao1 ~ Type, p.adjust.method = "BH") %>%
  rstatix::add_significance()


#Shannon diversity index
kruskal.test(Shannon ~ Type, data = data.frame(sample_data(PS99_ps0_div_FL)))

Shannon_Wilcox_FL <- as(sample_data(PS99_ps0_div_FL),"data.frame")   %>%
  rstatix::wilcox_test(Shannon ~ Type, p.adjust.method = "BH") %>%
  rstatix::add_significance()

#PA
#Chao1 richness
kruskal.test(Chao1 ~ Type, data = data.frame(sample_data(PS99_ps0_div_PA)))

Chao1_Wilcox_PA <- as(sample_data(PS99_ps0_div_PA),"data.frame")   %>%
  rstatix::wilcox_test(Chao1 ~ Type,p.adjust.method = "BH") %>%
  rstatix::add_significance()


#Shannon diversity index
kruskal.test(Shannon ~ Type, data = data.frame(sample_data(PS99_ps0_div_PA)))

Shannon_Wilcox_PA <- as(sample_data(PS99_ps0_div_PA),"data.frame")   %>%
  rstatix::wilcox_test(Shannon ~ Type, p.adjust.method = "BH") %>%
  rstatix::add_significance()


#plot alpha diversity
PS99_phy0_alpha.m <- merge(meta,PS99_phy0_alpha,by ="Row.names") %>%
  select(Region,Type, Fraction, Observed, Chao1, Shannon,InvSimpson)%>%
  melt(id.vars = c("Type","Fraction","Region"))%>%
  filter(Type!="Sediment")



alpha.p<- ggplot(PS99_phy0_alpha.m, aes(x = Type, y = value, fill= Fraction)) +
  labs(x = "Water layer", y = "Alpha diversity")+
  scale_fill_manual(values =c("yellow","darkgreen"))+
  geom_boxplot(outlier.color = NULL, notch = FALSE)+
  facet_wrap(~variable, ncol =2, scales = "free")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  coord_cartesian(clip="off")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggsave("./figures/alpha_div.png", 
       plot = alpha.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))

#unload libraries
detach("package:phyloseq")
detach("package:dplyr")


