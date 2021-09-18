#####################################
#import trajectories of in situ measured velocities
#####################################
traj_EGC <- as.data.frame(read.csv("./Tables/trajectories/drifter_Arc40_start2016_EG_bottom_speed52md_ice_dist_len.txt",
                                   sep = " ", header = F))
names(traj_EGC) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")

traj_WSC <- read.csv("./Tables/trajectories/drifter_Arc40_start2016_HG_bottom_speed29md_ice_dist_len.txt",
                     sep = " ", header = F)
names(traj_WSC) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")
traj_N <- read.csv("./Tables/trajectories/drifter_Arc40_start2016_N_bottom_speed52md_ice_dist_len.txt",
                   sep = " ", header = F)
names(traj_N) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")

#subset to the productive season (March-July)
traj_EGC_sub <- traj_EGC[traj_EGC$day %in% c(59:213),]
traj_EGC_sub$Region <- "EGC"
traj_WSC_sub <- traj_WSC[traj_WSC$day %in% c(59:213),]
traj_WSC_sub$Region <- "WSC"
traj_N_sub <- traj_N[traj_N$day %in% c(59:213),]
traj_N_sub$Region <- "N"

#merge togather
traj_sub <- rbind(traj_EGC_sub,traj_WSC_sub,traj_N_sub)

#adjust presence of ice to 15% conc.
traj_sub[traj_sub$ice.conc>15,]$ice.pres <- 1
traj_sub[traj_sub$ice.conc<15,]$ice.pres <- 0

#Generate a dataframe that keeps a summary of the data
df_summary <- data.frame(Region=levels(factor(traj_sub$Region)), N=tapply(traj_sub$length, traj_sub$Region, length), length=tapply(traj_sub$length, traj_sub$Region, mean), length.median=tapply(traj_sub$length, traj_sub$Region, median),dist=tapply(traj_sub$dist, traj_sub$Region, mean),dist.median=tapply(traj_sub$dist, traj_sub$Region, median))

overview_by_region <- merge(df_summary, aggregate(ice.pres~Region, traj_sub,sum), by = "Region")
overview_by_region$velo <- "measured"

# Or in short:
prop.test(c(overview_by_region[overview_by_region$Region =="EGC",]$ice.pres,
            overview_by_region[overview_by_region$Region =="WSC",]$ice.pres,
            overview_by_region[overview_by_region$Region =="N",]$ice.pres), 
          c(overview_by_region[overview_by_region$Region =="EGC",]$N,
            overview_by_region[overview_by_region$Region =="WSC",]$N,
            overview_by_region[overview_by_region$Region =="N",]$N),
          correct = TRUE)

#####################################
#import trajectories of 20 m/d velocity
#####################################
#import files
traj_EGC_20 <- as.data.frame(read.csv("./Data/trajectories/drifter_Arc40_start2016_EG_bottom_speed20md_ice_dist_len.txt",
                                      sep = " ", header = F))
names(traj_EGC_20) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")
traj_WSC_20 <- read.csv("./Data/trajectories/drifter_Arc40_start2016_HG_bottom_speed20md_ice_dist_len.txt",
                        sep = " ", header = F)
names(traj_WSC_20) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")
traj_N_20 <- read.csv("./Data/trajectories/drifter_Arc40_start2016_N_bottom_speed20md_ice_dist_len.txt",
                      sep = " ", header = F)
names(traj_N_20) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")

#subset to the spring and summer 
traj_EGC_20_sub <- traj_EGC_20[traj_EGC_20$day %in% c(59:213),]
traj_EGC_20_sub$Region <- "EGC"
traj_WSC_20_sub <- traj_WSC_20[traj_WSC_20$day %in% c(59:213),]
traj_WSC_20_sub$Region <- "WSC"
traj_N_20_sub <- traj_N_20[traj_N_20$day %in% c(59:213),]
traj_N_20_sub$Region <- "N"

traj_sub_20 <- rbind(traj_EGC_20_sub,traj_WSC_20_sub,traj_N_20_sub)

#adjust presence of ice to 15% conc.
traj_sub_20[traj_sub_20$ice.conc>15,]$ice.pres <- 1
traj_sub_20[traj_sub_20$ice.conc<15,]$ice.pres <- 0

#Generate a dataframe that keeps a summary of the data
df_summary_20 <- data.frame(Region=levels(factor(traj_sub_20$Region)), N=tapply(traj_sub_20$length, traj_sub_20$Region, length), length=tapply(traj_sub_20$length, traj_sub_20$Region, mean), length.median=tapply(traj_sub_20$length, traj_sub_20$Region, median),dist=tapply(traj_sub_20$dist, traj_sub_20$Region, mean),dist.median=tapply(traj_sub_20$dist, traj_sub_20$Region, median))

overview_by_region_20 <- merge(df_summary_20, aggregate(ice.pres~Region, traj_sub_20,sum), by = "Region")
overview_by_region_20$velo <- 20

# Or in short:
prop.test(c(overview_by_region_20[overview_by_region_20$Region =="EGC",]$ice.pres,overview_by_region_20[overview_by_region_20$Region =="WSC",]$ice.pres,overview_by_region_20[overview_by_region_20$Region =="N",]$ice.pres), 
          c(overview_by_region_20[overview_by_region_20$Region =="EGC",]$N,overview_by_region_20[overview_by_region_20$Region =="WSC",]$N,overview_by_region_20[overview_by_region_20$Region =="N",]$N),
          correct = TRUE)

#####################################
#import trajectories of 60 m/d velocity
#####################################
#import files
traj_EGC_60 <- as.data.frame(read.csv("./Data/trajectories/drifter_Arc40_start2016_EG_bottom_speed60md_ice_dist_len.txt",
                                      sep = " ", header = F))
names(traj_EGC_60) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")
traj_WSC_60 <- read.csv("./Data/trajectories/drifter_Arc40_start2016_HG_bottom_speed60md_ice_dist_len.txt",
                        sep = " ", header = F)
names(traj_WSC_60) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")
traj_N_60 <- read.csv("./Data/trajectories/drifter_Arc40_start2016_N_bottom_speed60md_ice_dist_len.txt",
                      sep = " ", header = F)
names(traj_N_60) <- c("long","lat","N","year","day","year2","ice.conc","ice.pres","dist","length")

#subset to the spring and summer 
traj_EGC_60_sub <- traj_EGC_60[traj_EGC_60$day %in% c(59:213),]
traj_EGC_60_sub$Region <- "EGC"
traj_WSC_60_sub <- traj_WSC_60[traj_WSC_60$day %in% c(59:213),]
traj_WSC_60_sub$Region <- "WSC"
traj_N_60_sub <- traj_N_60[traj_N_60$day %in% c(59:213),]
traj_N_60_sub$Region <- "N"

traj_sub_60 <- rbind(traj_EGC_60_sub,traj_WSC_60_sub,traj_N_60_sub)

#adjust presence of ice to 15% conc.
traj_sub_60[traj_sub_60$ice.conc>15,]$ice.pres <- 1
traj_sub_60[traj_sub_60$ice.conc<15,]$ice.pres <- 0

#Generate a dataframe that keeps a summary of the data
df_summary_60 <- data.frame(Region=levels(factor(traj_sub_60$Region)), N=tapply(traj_sub_60$length, traj_sub_60$Region, length), length=tapply(traj_sub_60$length, traj_sub_60$Region, mean), length.median=tapply(traj_sub_60$length, traj_sub_60$Region, median),dist=tapply(traj_sub_60$dist, traj_sub_60$Region, mean),dist.median=tapply(traj_sub_60$dist, traj_sub_60$Region, median))

overview_by_region_60 <- merge(df_summary_60, aggregate(ice.pres~Region, traj_sub_60,sum), by = "Region")
overview_by_region_60$velo <- 60

# Or in short:
prop.test(c(overview_by_region_60[overview_by_region_60$Region =="EGC",]$ice.pres,overview_by_region_60[overview_by_region_60$Region =="WSC",]$ice.pres,overview_by_region_60[overview_by_region_60$Region =="N",]$ice.pres), 
          c(overview_by_region_60[overview_by_region_60$Region =="EGC",]$N,overview_by_region_60[overview_by_region_60$Region =="WSC",]$N,overview_by_region_60[overview_by_region_60$Region =="N",]$N),
          correct = TRUE)


#####################################
#merge all results together
#####################################
particles_trajectories <- rbind(overview_by_region_60,overview_by_region_20,overview_by_region)
particles_trajectories$from.ice.pro <- particles_trajectories$ice.pres/particles_trajectories$N

particles_trajectories

write.csv(particles_trajectories, "./Tables/trajectories_summary.csv")
