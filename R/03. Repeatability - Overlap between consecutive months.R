###Utilization distributions
##Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
##Purpose: The code calculates kernel densities curves and utilization distributions of seals foraging patches
##Output: For each seal UD maps (95%, 75%, 50% and 25%) and overlap in ARS UD between consecutive trips (95%, 75%, 50% and 25%)

#Created on: 19/01/2021
#Updated on: 13/05/2021

library(pacman)
p_load(ggplot2, sf, rgdal, sp, adehabitatHR, tidyverse, magrittr, lubridate, ggspatial)

## Data preparation ----------------------------------------------------------------------------
#Read in trip summaries
trip.summaries <- read.table(here::here("Dryad", "pv64-2017_trip_summaries.txt"), sep="\t",header=TRUE)
trip.summaries$Trip_Code <- format(trip.summaries$Trip_Code, nsmall=3)

dat <- read.csv(here::here("Output", "Model 1 - HMM dive batches classified.csv"), header=TRUE)
dat <- dat %>% select(ID,seal_ID,PTT,start.time,end.time,batch.duration,batch.start.lon,batch.start.lat,
                      batch.end.long,batch.end.lat,x,y,HMMstate,state)

dat$ID <- format(dat$ID, nsmall=3)

#For each seall extract month 1 (= April) and month 2 (= May)
trip.summaries$month <- month(as.Date(trip.summaries$Trip_Start))
#months <- c(4,5)
# trips <- trip.summaries[which(trip.summaries$month %in% months),]
# trips.list <- trips$Trip_Code
month1 <- trip.summaries$Trip_Code[which(trip.summaries$month==4)]
month2 <- trip.summaries$Trip_Code[which(trip.summaries$month==5)]

dat$month <- ifelse(dat$ID %in% month1, "1", NA)
dat$month <- ifelse(dat$ID %in% month2, "2", dat$month)

dat <- dat[which(dat$month>0),]

#read in the coastaline shapefile as you will need it for plotting
coastline <- st_read(here::here("Datasets","Coastline_UTM30","Coastline_UTM30.shp"))

#Create grid for which you are going to calculate the kernel UD. 
#This grid is roughly 250mx250m and fits within the original 1x1km MF grid
#Try with three different grid sizes: 100x100 m, 500x500 m and 1x1 km
grid <- readOGR(here::here("Datasets","Moray_Firth_1km_grid_shapefile","MF_grid_1km_UTM30.shp"))
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

df <- data.frame(seal_ID = unique(dat$seal_ID), month1_n_trips = NA,
                 month2_n_trips = NA, BA_overlap = NA)

for(y in 27:length(df$seal_ID)){
  tmp <- dat[which(dat$seal_ID==df$seal_ID[y]),]
  
  #Extract all locations classified as ARS
  tmpfr <- tmp[which(tmp$state=="ARS"),]
  
  #Create list of Trip ID that don't have enough observations
  # splst <- tmpfr %>% group_by(ID) %>% summarize(n = n()) %>% arrange(n) %>%
  #   filter(n > 9) %>% dplyr::select(ID) %>% extract2(1) %>% as.character()
  # #Remove trips that don't have eough observations
  #tmpfr <- tmpfr %>% filter(ID %in% splst) %>% droplevels()
  
  #Make ARS location a spatial point dataframe to maintain the 
  spfor <- SpatialPointsDataFrame(coordinates(cbind(tmpfr$x, tmpfr$y)),  data=tmpfr, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
  
  df$month1_n_trips[y] <- length(unique(tmpfr$ID[which(tmpfr$month==1)]))
  df$month2_n_trips[y] <- length(unique(tmpfr$ID[which(tmpfr$month==2)]))
  
  #Calculate kernel density curve
  #h <- c(50,250,500,1000,1500,2000)
  kd <- kernelUD(spfor[,15], h=2000, grid=xy, kern=c("bivnorm"))
  
  #For each trip extract the 95% utilization distribution
  utilizations <- list()
   for(x in 1:length(kd)){
     ud1 <- getverticeshr.estUD(kd[[x]], 95)
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
      geom_polygon(data= udsp, aes(x=x,y=y, fill=as.factor(Id), group=as.factor(group)), alpha=.6)+
      scale_fill_manual(values=c("#FFB1A1", "#70F5D6"), labels=c("April","May"), name="Months")+
      geom_point(data= tmpfr, aes(x=x, y=y), col="black")+
      xlim(min(udsp$x)-1000,max(udsp$x)+1000)+ylim(min(udsp$y)-1000,max(udsp$y)+1000)+
      xlab("Longitude")+ylab("Latitude")+
      theme_classic()+
      theme(legend.position = "top")
    ggsave(plot=ARS.UD.plot , filename=here::here("Figures","Repeatability", paste0("ARS_95_UD_seal",seal.Ids[y],"_2017.png")),
           device="png", height=6, width = 8, unit="in")
    
    #Calculate overlap between consecutive trips
    overlap <- kerneloverlaphr(kd, method="BA", percent=95, conditional = TRUE)
    overlap <- as.data.frame(overlap)
    df$BA_overlap[y] <- overlap[1,2]
  print(y)
}

write.csv(df, here::here("Output","Repeatability overlap between April and May.csv"), row.names=FALSE)

## Overlap null distribution ---------------------------------------------------------------------------------------
tmp <- unique(dat$seal_ID)
tmp <- tmp[-which(tmp==59 | tmp==280)]
df <- data.frame(seala = NA, sealb= NA)
for(x in 1:length(tmp)){
  seala <- rep(tmp[x], length(tmp)-1)
  sealb <- tmp[which(tmp != tmp[x])]
  cols <- cbind(seala,sealb)
  df <- rbind(df, cols)
}
df <- df[-1,]

df$BA_overlap <- NA

#Calculate the overlap between seal's 1 May and every other seal's April
for(y in 1:length(df$seala)){
  tmp1 <- dat[which(dat$seal_ID==df$sealb[y] & dat$month==1 & dat$state=="ARS"),]
  tmp2 <- dat[which(dat$seal_ID==df$seala[y] & dat$month==2 & dat$state=="ARS"),]
  
  #Extract all locations classified as ARS
  tmpfr <- rbind(tmp1, tmp2)
  
  #Make ARS location a spatial point dataframe to maintain the 
  spfor <- SpatialPointsDataFrame(coordinates(cbind(tmpfr$x, tmpfr$y)),  data=tmpfr, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
  
  #Calculate kernel density curve
  kd <- kernelUD(spfor[,15], h=2000, grid=xy, kern=c("bivnorm"))
  #kd1 <- kernelUD(spfor[,1], h=5000, grid=xy3, kern=c("bivnorm"))
  
  #Calculate overlap between consecutive trips
  overlap <- kerneloverlaphr(kd, method="BA", percent=95, conditional = TRUE)
  overlap <- as.data.frame(overlap)
  df$BA_overlap[y] <- overlap[1,2]
  print(y)
}

write.csv(df, here::here("Output","Repeatability - null distribution overlap.csv"), row.names=FALSE)
