###Utilization distributions
##Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
##Purpose: The code calculates kernel densities curves and utilization distributions of seals foraging patches
##Output: For each seal UD maps (95%, 75%, 50% and 25%) and overlap in ARS UD between consecutive trips (95%, 75%, 50% and 25%)

#Created on: 17/12/2020
#Updated on: 22/12/2020

library(pacman)
p_load(ggplot2, sf, rgdal, sp, adehabitatHR, tidyverse, magrittr, lubridate, ggspatial)

## Data preparation ----------------------------------------------------------------------------
#Read in trip summaries
trip.summaries <- read.csv(here::here("Dataset", "2017","pv64-2017_trip_summaries.csv"), header=TRUE)
trip.summaries$Trip_Code <- format(trip.summaries$Trip_Code, nsmall=3)

#Read dataset with the interpolated 30 minutes location and state classification from the HMM
#dat <- read.csv(here::here("Dataset", "HMM_classified_interpolated_trips.csv"), header=TRUE)
#dat <- dat[,c(2,5:8,10,11)]

dat <- read.csv(here::here("Output", "Model 1 - HMM dive batches classified.csv"), header=TRUE)
dat <- dat[,c(2,5:15,17,18)]

dat$ID <- format(dat$ID, nsmall=3)

#Extract the foraging trips from 2017
dat$year <- year(dat$start.time)
dat <- dat[which(dat$year==2017),]

#read in the coastaline shapefile as you will need it for plotting
coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")

#Create grid for which you are going to calculate the kernel UD. 
#This grid is roughly 250mx250m and fits within the original 1x1km MF grid
#Try with three different grid sizes: 100x100 m, 500x500 m and 1x1 km
grid <- readOGR("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/MF_grid/MF_grid_1km_UTM30.shp")
x <- seq(grid@bbox[1]+250, grid@bbox[3]-250, by=500)
y <- seq(grid@bbox[2]+250, 6600733, by=500)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

## UD and overlap -----------------------------------------------------------------------------
#Create utilization distribution and calculate foraging trips overlap wihtin individuals
seal.Ids <- unique(dat$seal_ID)
#Error at seal number 7 because she is the one that goes to Orkney
#Error at seal 8 because outlier GPS locations and because she goes to orkney again
for(y in 1:length(seal.Ids)){
  tmp <- dat[which(dat$seal_ID==seal.Ids[y]),]
  
  #Extract all locations classified as ARS
  tmpfr <- tmp[which(tmp$state=="ARS"),]
  
  #Create list of Trip ID that don't have enough observations
  splst <- tmpfr %>% group_by(ID) %>% summarize(n = n()) %>% arrange(n) %>%
  filter(n > 9) %>% dplyr::select(ID) %>% extract2(1) %>% as.character()
  #Remove trips that don't have eough observations
  tmpfr <- tmpfr %>% filter(ID %in% splst) %>% droplevels()
  
  #Make ARS location a spatial point dataframe to maintain the 
  spfor <- SpatialPointsDataFrame(coordinates(cbind(tmpfr$x, tmpfr$y)),  data=tmpfr, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
  
  #Calculate kernel density curve
  #h <- c(50,250,500,1000,1500,2000)
  kd <- kernelUD(spfor[,1], h=2000, grid=xy, kern=c("bivnorm"))
  #kd1 <- kernelUD(spfor[,1], h=5000, grid=xy3, kern=c("bivnorm"))
  
  #For each trip extract the 95% utilization distribution
  utilizations.list <- c(95,75,50,25)
  for(z in 1:length(utilizations.list)){
    
    utilizations <- list()
    for(x in 1:length(kd)){
      ud1 <- getverticeshr.estUD(kd[[x]], utilizations.list[z])
      utilizations[[x]] <- ud1
      }
    
    #Create kd name list
    kd_names <- names(kd)
    
    #Extract x% UD polygons list
     #Extract the first polygon
    if(length(utilizations[[1]]@polygons[[1]]@Polygons)==1){
      udsp <- utilizations[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
        as.data.frame(.) %>%
         mutate(Id = rep(kd_names[1], length(V1)), poly = rep(1, length(V1)))
      } else {
        n <- length(utilizations[[1]]@polygons[[1]]@Polygons)
        udsp <- utilizations[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
          as.data.frame(.) %>%
           mutate(Id = rep(kd_names[1], length(V1)), poly = rep(1, length(V1)))
        for(i in 2:n){
          udsp1 <- utilizations[[1]]@polygons[[1]]@Polygons[[i]]@coords %>% 
            as.data.frame(.) %>%
            mutate(Id = rep(kd_names[1], length(V1)), poly = rep(i, length(V1)))
          udsp <- rbind(udsp, udsp1)
          }
        }
  
    #in a loop extract and combine all the other polygons
    for(x in 2:length(kd_names)){
      if(length(utilizations[[x]]@polygons[[1]]@Polygons)==1){
        udtmp <- utilizations[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
          as.data.frame(.) %>%
          mutate(Id = rep(kd_names[x], length(V1)), poly = rep(1, length(V1)))
        }  else {
          n <- length(utilizations[[x]]@polygons[[1]]@Polygons)
          udtmp <- utilizations[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
                    as.data.frame(.) %>%
             mutate(Id = rep(kd_names[x], length(V1)), poly = rep(1, length(V1)))
          for(i in 2:n){
          udtmp1 <- utilizations[[x]]@polygons[[1]]@Polygons[[i]]@coords %>% 
            as.data.frame(.) %>%
            mutate(Id = rep(kd_names[x], length(V1)), poly = rep(i, length(V1)))
          udtmp <- rbind(udtmp, udtmp1)
          }
        }
      udsp <- rbind(udsp, udtmp)
      }

    colnames(udsp) <- c("x","y","Id","poly")
    udsp$group <- paste(udsp$Id, udsp$poly)

  #Create and save x% UD maps 
      ARS.UD.plot <- ggplot()+
         annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
         geom_polygon(data= udsp, aes(x=x,y=y, fill=as.factor(Id), group=as.factor(group)), alpha=.4)+
         geom_point(data= tmpfr, aes(x=x, y=y, col=as.factor(HMMstate)), col="black")+
         xlim(min(udsp$x)-1000,max(udsp$x)+1000)+ylim(min(udsp$y)-1000,max(udsp$y)+1000)+
         xlab("Longitude")+ylab("Latitude")+
         theme_classic()+
         theme(legend.position="none")
       ggsave(plot=ARS.UD.plot , filename=here::here("Figures","ARS utilization distributions by indviduals", paste0("ARS_", utilizations.list[z],"UD"), paste0("ARS_",utilizations.list[z],"UD_seal",seal.Ids[y],"_2017.png")),
         device="png", height=6, width = 8, unit="in")
  
  #Calculate overlap between consecutive trips
       overlap <- kerneloverlaphr(kd, method="BA", percent=utilizations.list[z], conditional = TRUE)
       overlap <- as.data.frame(overlap)
  
  #Create dataframe with overlap between consecutive trips
       trips.ID <- unique(tmpfr$ID)
       df <- data.frame(trip1 = trips.ID[1:(length(trips.ID)-1)] , trip2 = trips.ID[2:length(trips.ID)], overlap =NA)
       for(x in 1:(length(trips.ID)-1)){
         df$overlap[x] <- overlap[x+1,x]
         }
       df$overlap <- as.numeric(df$overlap)
       df$seal_ID <- rep(seal.Ids[y], length(df$trip1))
       
       seal.trips <- trip.summaries[which(trip.summaries$ID==seal.Ids[y]),]
  
       for(i in 1:length(df$trip1)){
         t1 <- which(seal.trips$Trip_Code == df$trip1[i])
         t2 <- which(seal.trips$Trip_Code == df$trip2[i])
    
         df$trip1.duration[i] <- seal.trips$Trip_Duration[t1]
         df$trip1.start[i] <- seal.trips$Trip_Start[t1]
         df$trip1.end[i] <- seal.trips$Trip_End[t1]
         df$trip1.max.distance.haulout[i] <- max(seal.trips$MaxDist_Start_Haulout[t1], seal.trips$MaxDist_End_Haulout[t1]) 
         df$trip1.max.depth[i] <- seal.trips$Max_Depth[t1]
    
    
         df$trip2.duration[i] <- seal.trips$Trip_Duration[t2]
         df$trip2.start[i] <- seal.trips$Trip_Start[t2]
         df$trip2.end[i] <- seal.trips$Trip_End[t1]
         df$trip2.max.distance.haulout[i] <- max(seal.trips$MaxDist_Start_Haulout[t2], seal.trips$MaxDist_End_Haulout[t2])
         df$trip2.max.depth[i] <- seal.trips$Max_Depth[t2]
         }
  
       df <- df[,c(4,1:3,5:14)]
       write.csv(df, here::here("Output","ARS trips overlap within individuals",paste0("ARS_", utilizations.list[z] ,"UD"), paste0("ARS_",utilizations.list[z],"UD_overlap_seal",seal.Ids[y],"_2017.csv")), row.names = FALSE)
       print(paste0("Seal ", seal.Ids[y]," ",utilizations.list[z],"% Utilization distribution done"))
       }
  print(y)
}


