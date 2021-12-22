## Create dive batches and filter foraging trips for Model 2
## Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
## Purpose: Select foraging trips and prepare dive bacthes for the foraging trips with accelerometer data
## Output: -csv with dataset
#Created on: 04/12/2020
#Updated on: 21/12/2021

#Load necessary packages
library(pacman)
p_load(sp, sf, magrittr, lubridate, tidyverse, geosphere, ggplot2,ggspatial)

## Create buffers around haul-out sites ----------------------------------------------------
#Load coastline needed for plotting
coastline <- st_read(here::here("Dryad","Coastline_UTM30","Coastline_UTM30.shp"))

#Create the two haul-out sites buffers to remove locations within 2km from the haul-out site
#Create Loch Fleet buffer based on two points
loch.fleet <- data.frame(location = rep("Loch Fleet",1), lat = 57.938128, lon = -4.035541) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326) %>% st_transform(., crs=32630)
loch.fleet1 <- data.frame(location = rep("Loch Fleet",1), lat = 57.941909, lon = -4.065843) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326) %>% st_transform(., crs=32630) 

# create a 2km buffer around loch fleet
LF_2k_buffer <- st_buffer(loch.fleet, dist=2000)
LF_2k_buffer1 <- st_buffer(loch.fleet1, dist=2000)

#Join the two buffers
buffer <- st_union (LF_2k_buffer, LF_2k_buffer1$geometry)

#South coast buffer
south.coast <- data.frame(location = rep("South coast",1), lat = 57.635, lon = -3.74) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326) %>% st_transform(., crs=32630) 
south.coast2 <- data.frame(location = rep("South coast",1), lat = 57.65, lon = -3.71) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326) %>%  st_transform(SC_sf1, crs=32630)

# create a 2km buffer around south coast haul-out site
SC_2k_buffer <- st_buffer(south.coast, dist=2000)
SC_2k_buffer1 <- st_buffer(south.coast2, dist=2000)

buffer2 <- st_union (SC_2k_buffer, SC_2k_buffer1$geometry)

## Data preparation ---------------------------------------------------------------------------------------------------
#select only the variables used the analysis
path <- here::here("Dryad","Outputs","Processed accelerometer parameters")
files <- paste0(path,"/",list.files(paste0(path)))
tables <- lapply(files, read.csv)
seal <- do.call(rbind, tables) 
seal <- seal %>% dplyr::select(ID,PTT,DS_DATE,DE_DATE,DIVE_DUR,SURF_DUR,MAX_DEP,START_LAT,START_LON,END_LAT,END_LON,PERCENT_AREA,
         Trip_No,Posn_in_Trip,Trip_Code,date,TAD,descent.speed,ascent.speed,B_PrCA_arch,B_pitch.diff20,
         B_overlap_A) %>%
  mutate(seal_ID = ID,
         trip_code = format(Trip_Code, nsmall=3),
         ID = trip_code,
         time = as.POSIXct(DS_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         time.end = as.POSIXct(DE_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"))

llcord <- SpatialPoints(seal[,c("START_LON","START_LAT")],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
seal$x <- attr(utmcoord, "coords")[,1]
seal$y <- attr(utmcoord, "coords")[,2]

## Filter trips and locations ----------------------------------------------------------------
##Remove locations within the haul-out site buffer
#Loch Fleet
seal_sf <- st_as_sf(seal, coords=c(x="x", y="y"), crs=32630)
# to remove data that are inside the buffer
seal_at_sea <- st_difference(seal_sf, buffer)
# convert sf object to dataframe 
seal <- seal_at_sea 

#Repeat for the South coast
seal_sf2 <- st_as_sf(seal, coords=c(x="x", y="y"), crs=32630) # already in UTM 30N sp no need to project it 
seal_at_sea2 <- st_difference(seal_sf2, buffer2)
seal <- seal_at_sea2 
st_geometry(seal) <- NULL
seal <- seal[,-c(28,29)]


##Filter the dataset for trips longer than 12 hours
trip <- read.table(here::here("Dryad","pv64-2017_trip_summaries.txt"),sep="\t", head=TRUE)
trip <- trip[-c(which(trip$PTT=="99999")),]
trip$Trip_Code <- format(trip$Trip_Code, nsmall=3)
long_trips <- trip$Trip_Code[which(trip$Trip_Duration>12)]
seal <- seal[seal$ID %in% long_trips,]


##Filter the dataset for trips occurring in the first week post tagging
#Read seal summary
seal.summary <- read.table(here::here("Dryad","pv64-2017_seal_summary.txt"),sep="\t", header=TRUE) %>%
  mutate(date.cap = as.Date(Tag_start, format="%d/%m/%Y"),
         first.week = date.cap + 7)

#Find the trips starting in the first week
trip$date <- as.Date(trip$Trip_Start)
trip$first.week <- NA
for(x in 1:length(trip$ID)){
  n <- which(seal.summary$ID == trip$ID[x])
  trip$first.week[x] <- paste(seal.summary$first.week[n])
}
trip$first.week <- as.Date(trip$first.week)
trip$to.keep <- ifelse(trip$date <= trip$first.week, 0, 1)
after.first.week.trips <- trip$Trip_Code[which(trip$to.keep==1 & trip$Trip_Duration >12)]
seal <- seal[which(seal$ID %in% after.first.week.trips),]



##Filter for round.trips at the same haul out location
trip$haulout.dist <- NA
for(x in 1:length(trip$Trip_Code)){
  trip$haulout.dist[x] <- distm(c(trip$Start_Haulout_Lon[x], trip$Start_Haulout_Lat[x]), c(trip$End_Haulout_Lon[x], trip$End_Haulout_Lat[x]), fun=distHaversine)
}
round_trips <- trip$Trip_Code[which(trip$to.keep==1 & trip$Trip_Duration >12 & trip$haulout.dist<=5000)]
seal <- seal[which(seal$ID %in% round_trips),]


#Remove trips with missing accelerometer data
forNA <- seal[which(is.na(seal$B_PrCA_arch) | is.na(seal$B_pitch.diff20)),]
tripsNA <- as.factor(unique(forNA$ID)) 

'%!in%' <- function(x,y)!('%in%'(x,y))
seal <- seal[which(seal$ID %!in% tripsNA),] 


## Create dive batches --------------------------------------------------------------------------------------------------
trips <- unique(seal$ID)
batch_series <- NA

i <- 4
x <- 0
for(i in 1:length(trips)){
  batch <- subset(seal, seal$ID==trips[i])
  Ndives <- ceiling(length(batch$ID)/5)
  Ndives <- Ndives+x
  batch$batch <- rep(x:Ndives, each= 5,length.out = length(batch$ID))
  batch_series <- c(batch_series, batch$batch)
  x <- length(batch$ID)+1
}
batch_series <- batch_series[-1]
seal$batch <- batch_series

seal$foraging.index <- seal$B_PrCA_arch + seal$B_pitch.diff20 - seal$B_overlap_A 

seal.batch <- data.frame(seal_ID = NA, PTT = NA, trip.N = NA, ID= NA, start.time = NA, 
                         end.time = NA, batch.duration_mins = NA,
                         batch.dive.dur_secs = NA, batch.surface.dur_secs = NA, batch.max.depth = NA,
                         batch.pos.trip.start = NA, batch.pos.trip.end = NA, 
                         batch.PrCA = NA, batch.benthic = NA, batch.foraging.index = NA, 
                         batch.start.lon = NA, batch.start.lat = NA, batch.end.long = NA, batch.end.lat = NA)
batch.tmp <- seal.batch

x <- 1
n <- 5
options(scipen = 999)
for(n in 1:(length(seal$ID)-1)){
  if(seal$batch[n]!=seal$batch[n+1]){
    y <- n 
    batch.tmp$seal_ID <- seal$seal_ID[x]
    batch.tmp$PTT <- seal$PTT[x]
    batch.tmp$trip.N <- seal$Trip_No[x]
    batch.tmp$ID <- seal$ID[x]
    
    batch.tmp$start.time <- paste(seal$time[x])
    batch.tmp$end.time <- paste(seal$time[y]+seal$SURF_DUR[y])
    options(digits=2)
    batch.tmp$batch.duration_mins <- difftime(seal$time[y], seal$time[x], units="mins")
    batch.tmp$batch.pos.trip.start <- seal$Posn_in_Trip[x]
    batch.tmp$batch.pos.trip.end <- seal$Posn_in_Trip[y]
    
    #Batch covariates
    batch.tmp$batch.dive.dur_secs <- mean(seal$DIVE_DUR[x:y])
    batch.tmp$batch.surface.dur_secs <- mean(seal$SURF_DUR[x:y])
    batch.tmp$batch.max.depth <- mean(seal$MAX_DEP[x:y])

    batch.tmp$batch.PrCA <- as.integer(mean(seal$B_PrCA_arch[x:y], na.rm=TRUE))
    batch.tmp$batch.benthic <- as.integer(mean(seal$B_pitch.diff20[x:y], na.rm =TRUE))
    batch.tmp$batch.foraging.index <- as.integer(mean(seal$foraging.index[x:y], na.rm =TRUE))
    
    #Batch Lat and Long
    options(digits=8)
    batch.tmp$batch.start.lon <- seal$START_LON[x]
    batch.tmp$batch.start.lat <- seal$START_LAT[x]
    batch.tmp$batch.end.lon <- seal$END_LON[y]
    batch.tmp$batch.end.lat <- seal$END_LAT[y]
    
    seal.batch <- rbind(seal.batch, batch.tmp)
    x <- y+1
  } 
}
seal.batch <- seal.batch[-c(1),]

llcord <- SpatialPoints(seal.batch[,c("batch.start.lon","batch.start.lat")],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
seal.batch$x <- attr(utmcoord, "coords")[,1]
seal.batch$y <- attr(utmcoord, "coords")[,2]

ggplot()+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_point(data=seal.batch, aes(x=x, y=y))+
  theme_classic()

write.table(seal.batch, here::here("Dryad","Outputs","Dive batches dataset - Accelerometer.txt"), 
            sep="\t", row.names=FALSE)
