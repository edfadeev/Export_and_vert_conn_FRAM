#load libraries
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")

source("scripts/functions.R")
source("scripts/col_palette.R")
#####################################
#import the long term sediment trap data
#####################################
long_term_flux <- read.csv("./Tables/Long_term_sed_traps_EF.csv", header =TRUE)

long_term_flux_month<- long_term_flux %>%  group_by(Year, Month) %>% 
  summarise(POC_month=mean(POC))

ggplot(long_term_flux_month,aes(x=as.factor(Month), y= POC_month, colour = as.factor(Year), group = as.factor(Year)))+ 
  #geom_boxplot() +
  #geom_jitter(size = 5)+
  geom_line()+
  scale_color_manual(values = unique(sample(tol21rainbow,14)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        legend.position = "right", 
        axis.title.x = element_blank())

#explore mean values
flux_mean <- long_term_flux %>% group_by(Month) %>% 
  dplyr::summarise(ice_distance.mean= mean(ice_distance),
                   ice_distance.se= se(ice_distance), 
                   POC.mean=mean(POC),
                   POC.se=se(POC),
                   n=n()) 

#####################################
# Extract flux peaks
#####################################
#extract only the spring (March-April)
long_term_flux_spring<- long_term_flux%>%
  filter(Month %in% 3:5) %>%
  mutate(Year = factor(Year),
         Month = factor(case_when(Month==3~"March",Month==4~"April",Month==5~"May"),
                        levels = c("March", "April","May")),
         Season = "spring")

#extract only the summer (July-August)
long_term_flux_summer<- long_term_flux%>%
  filter(Month %in% 7:9) %>% 
  mutate(Year = factor(Year),
         Month = factor(case_when(Month==6~"June",Month==7~"July",Month==8~"August",Month==9~"September"),
                        levels = c("June","July", "August","September")),
         Season = "summer")

#extract only the winter (October-February)
long_term_flux_winter<- long_term_flux%>%
  filter(Month %in% c(1,2,10:12)) %>% 
  mutate(Year = factor(Year),
         Month = factor(case_when(Month==10~"October",Month==11~"November",Month==12~"December",Month==1~"January",Month==2~"February"),
                        levels = c("October","November", "December","January","February")),
         Season = "Winter")

#####################################
# Correspondance between Years, ice-distance and POC
#####################################
#both peaks
kruskal.test(ice_distance ~ Season, data = rbind(long_term_flux_spring,long_term_flux_summer))
kruskal.test(POC ~ Season, data = rbind(long_term_flux_spring,long_term_flux_summer))

#spring
kruskal.test(ice_distance ~ Year, data = long_term_flux_spring)
kruskal.test(POC ~ Year, data = long_term_flux_spring)

#summer
kruskal.test(ice_distance ~ Year, data = long_term_flux_summer)
kruskal.test(POC ~ Year, data = long_term_flux_summer)

#####################################
# Fit an exponential model to spring flux
#####################################
# Estimate the rest parameters using a linear model
model.0 <- lm(log(POC) ~ ice_distance, data=long_term_flux_spring)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
summary(model.0)

######## Calculate the linear model predictions for the figure
icevalues <- long_term_flux_spring$ice_distance
POC.exponential2 <- exp(predict(model.0,list(ice_distance= icevalues)))
long_term_flux_spring$POC_reg<-POC.exponential2

#####################################
# Generate the combined figure
#####################################

# Final plot
long_term_flux.p<-ggplot(data=long_term_flux_spring)+
  geom_point(aes(y=POC, x=ice_distance, color=Year, shape=Month), size = 10)+
  geom_line(aes(x=ice_distance,y=POC_reg), color="red")+
  annotate("text", x=70, y=50, label = paste("y ==", signif(alpha.0,3),"*e^{",signif(beta.0,3),"*x}"), parse =TRUE)+
  annotate("text", x=70, y=47, label = paste("R^{2}==",signif(summary(model.0)$adj.r.squared, 3)), parse =TRUE)+
  annotate("text", x=70, y=44, label = paste("p ==",signif(summary(model.0)$coefficients[2,4],3)), parse =TRUE)+
  xlab("ice distance (km)")+
  ylab(expression(paste("POC flux (mg.m"^"-2"*".d"^"-1"*")")))+
  scale_color_manual(values = unique(sample(tol21rainbow,14)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        legend.position = "right", 
        axis.title.x = element_blank())

ggsave("./figures/sed_trap_flux.pdf", 
       plot = long_term_flux.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



