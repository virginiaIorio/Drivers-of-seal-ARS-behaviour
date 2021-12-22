###Utilization distributions and overlap for all the individuals between two consecutive months
##Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
##Purpose: The code calculates kernel densities curves and utilization distributions of seals foraging patches
##Output: We first used the methods by Lascelles et al. 2016 to estimate the most appropriate 
#         smoothing paramter h.
#         For each seal we first calculate the overlap in utilization distribution between
#         two consecutive months (here April and May).
#         The we calculate a null distribution overlap between a randomised pair-wise comparison
#         with another individual in the population. 

#Created on: 19/01/2021
#Updated on: 16/12/2021

library(pacman)
p_load(ggplot2, sf, rgdal, sp, adehabitatHR, tidyverse, magrittr, lubridate, ggspatial, adehabitatLT)

## Data preparation ----------------------------------------------------------------------------
#Read in trip summaries
trip.summaries <- read.table(here::here("Dryad", "pv64-2017_trip_summaries.txt"), sep="\t",header=TRUE)
trip.summaries$Trip_Code <- format(trip.summaries$Trip_Code, nsmall=3)

dat <- read.table(here::here("Dryad","Outputs", "Model 1 - HMM dive batches classified.txt"), sep="\t", header=TRUE)
dat <- dat %>%dplyr::select(ID,seal_ID,PTT,start.time,end.time,batch.duration,batch.start.lon,batch.start.lat,
                      batch.end.lon,batch.end.lat,x,y,HMMstate,state)

dat$ID <- format(dat$ID, nsmall=3)

#For each seall extract month 1 (= April) and month 2 (= May)
trip.summaries$month <- month(as.Date(trip.summaries$Trip_Start))
month1 <- trip.summaries$Trip_Code[which(trip.summaries$month==4)]
month2 <- trip.summaries$Trip_Code[which(trip.summaries$month==5)]

dat$month <- ifelse(dat$ID %in% month1, "1", NA)
dat$month <- ifelse(dat$ID %in% month2, "2", dat$month)

dat <- dat[which(dat$month>0),]

#read in the coastaline shapefile as you will need it for plotting
coastline <- st_read(here::here("Dryad","Coastline_UTM30","Coastline_UTM30.shp"))

#Create grid for which you are going to calculate the kernel UD. 
#This grid is roughly 500mx500m and fits within the original 1x1km MF grid
grid <- readOGR(here::here("Dryad","Moray_Firth_1km_grid_shapefile","MF_grid_1km_UTM30.shp"))
x <- seq(grid@bbox[1]+250, grid@bbox[3]-250, by=500)
y <- seq(grid@bbox[2]+250, 6600733, by=500)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE


## Determine h ======================================
# Use the methods of Lascelles et al. 2016 in Diversity and distribution to define h
# This next section of the code was copied from the appendix S5 of Lascelles et al. 2016
# https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12411
time <- as.POSIXct(dat$start.time, format="%Y-%m-%d %H:%M:%S", tz="UTC")
traj <- as.ltraj(data.frame(dat$x, dat$y),date=time,id=dat$seal_ID, typeII = TRUE)

Scales <- c(seq(0.5, 15, 0.5))
Scales <- Scales*1000
fpt.out <- fpt(traj, radii = Scales, units = "seconds")
fpt.scales <- varlogfpt(fpt.out, graph = FALSE)
Temp <- as.double(fpt.scales[1,])
plot(Scales, Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)))
Peak <- "Flexible"

ars.scales <- NULL
UIDs <- unique(dat$seal_ID)
for(i in 1:length(UIDs)){
  if(length(Scales) == length(which(is.na(fpt.scales[i,])))) {print(paste("Warning: ID", UIDs[i], "is smaller than smallest scale and will be ignored")); next}
  Temp <- as.double(fpt.scales[i,])
  #lines(Scales,Temp)
  plot(Scales, Temp, type="l")
  
  q <- which(!is.na(Temp))
  p <- 2
  while(!is.na(Temp[q[p]]) & Temp[q[p]] < Temp[q[p-1]] & q[p] != length(Temp)) {p <- p + 1}
  while(!is.na(Temp[q[p]]) & Temp[q[p]] > Temp[q[p-1]]) {p <- p + 1}
  
  rfpt <- Scales[q[p-1]]
  if(suppressWarnings(min(which(is.na(Temp))) == p)) {print(paste("ID", UIDs[i], "has no peak")); next}
  FirstPeak <- Scales[q[p-1]]
  MaxPeak <- Scales[which(Temp == max(Temp[q[p-1]:length(Temp)], na.rm=T))]
  if(Peak == "Flexible") {
    if(FirstPeak < MaxPeak[1])
    {
      MaxPeak <- MaxPeak[MaxPeak >= FirstPeak]
      ifelse(MaxPeak[1] < FirstPeak + (max(Scales)/3), ars.sc <- MaxPeak[1], ars.sc <- FirstPeak)
    }  else  {ars.sc <- FirstPeak}
  }
  if(Peak == "Max") {ars.sc <- MaxPeak}
  if(Peak == "First")  {ars.sc <- FirstPeak}
  if(Peak == "User")
  {
    print("Select Peak on Graph")
    N <- identify(Scales, Temp, n=1)
    ars.sc <- Scales[N]
  }
  abline(v=ars.sc, col="red", lty=2)
  ars.scales <- c(ars.scales, ars.sc)
  #print(ars.sc)
  #readline("proceed?")
}

AprScale <- mean(ars.scales)
AprScale <- round(AprScale/1000,3)


jpeg(filename = here::here("Figures","Repeatability - smoothing parameter selection.jpeg"), 
      width = 1800, height = 1580, res=300)
plot((Scales/1000), Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)), xlab="Scales (km)", ylab="")
for(i in 1:length(UIDs)){
  Temp <- as.double(fpt.scales[i,])
  lines((Scales/1000),Temp)
}
abline(v=ars.scales/1000, col="red", lty=2)
abline(v=AprScale, col="darkred", lty=1, lwd=3)
#print(ars.scales)
#print(AprScale)
text(max(Scales/1000)/2, 5, paste(AprScale, "km"), col="darkred", cex=3)
dev.off()



## Utilization distribution ============================================
#Create utilization distribution and calculate foraging trips overlap wihtin individuals
seal.Ids <- unique(dat$seal_ID)

#Remove seal 59 and 280 because there are not trips in May
dat <- dat[-which(dat$seal_ID==59),]
dat <- dat[-which(dat$seal_ID==280),]

df <- data.frame(seal_ID = unique(dat$seal_ID), month1_n_trips = NA,
                 month2_n_trips = NA, BA_overlap = NA)

#The code calculates the overlap between the 50% UD for each individuals between May and April, and
# saves a map of the two UDs.
for(y in 1:length(df$seal_ID)){
  tmp <- dat[which(dat$seal_ID==df$seal_ID[y]),]
  
  #Extract all locations classified as ARS
  tmpfr <- tmp[which(tmp$state=="ARS"),]
  
  #Make ARS location a spatial point dataframe to maintain the 
  spfor <- SpatialPointsDataFrame(coordinates(cbind(tmpfr$x, tmpfr$y)),  data=tmpfr, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
  
  df$month1_n_trips[y] <- length(unique(tmpfr$ID[which(tmpfr$month==1)]))
  df$month2_n_trips[y] <- length(unique(tmpfr$ID[which(tmpfr$month==2)]))
  
  #Calculate kernel density curve
  kd <- kernelUD(spfor[,15], h=5859, grid=xy, kern=c("bivnorm"))
  
  #For each trip extract the 95% utilization distribution
  utilizations <- list()
   for(x in 1:length(kd)){
     ud1 <- getverticeshr.estUD(kd[[x]], 50)
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
    ggsave(plot=ARS.UD.plot , filename=here::here("Figures","Repeatability","New plots", paste0("ARS_95_UD_seal",seal.Ids[y],"_2017.png")),
           device="png", height=6, width = 8, unit="in")
    
    #Calculate overlap between consecutive trips
    overlap <- kerneloverlaphr(kd, method="BA", percent=50, conditional = TRUE)
    overlap <- as.data.frame(overlap)
    df$BA_overlap[y] <- overlap[1,2]
  print(y)
}
hist(df$BA_overlap)

write.table(df, here::here("Dryad","Outputs","Repeatability - overlap between April and May.txt"),  sep="\t", row.names=FALSE)

## Overlap null distribution ---------------------------------------------------------------------------------------
#Randomized pair-wise comparison between each individual distribution in May with a randomly
# selected individual in April
tmp <- unique(dat$seal_ID)
df <- data.frame(seala = tmp, sealb= NA)

for(x in 1:length(tmp)){
  tmp2 <- tmp[-which(tmp==tmp[x])]
  sealb <- sample(tmp2, 1)
  df$sealb[x] <- sealb
}

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
  kd <- kernelUD(spfor[,15], h=5859, grid=xy, kern=c("bivnorm"))
  
  #Calculate overlap between consecutive trips
  overlap <- kerneloverlaphr(kd, method="BA", percent=50, conditional = TRUE)
  overlap <- as.data.frame(overlap)
  df$BA_overlap[y] <- overlap[1,2]
  print(y)
}
hist(df$BA_overlap)
write.table(df, here::here("Dryad","Outputs","Repeatability - null distribution overlap.txt"), sep="\t", row.names=FALSE)

