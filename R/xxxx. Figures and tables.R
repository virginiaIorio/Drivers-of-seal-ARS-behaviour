### Model 1 Figures ------------------------------------------------------------------------
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose: 
# Output: 
# Created on: 12/01/2021
# Updated on: 13/01/2021

## Load packages --------------------------------------------------------------------------------------
library(pacman)
p_load(ggplot2, ggspatial, ggsn, ggpubr, sf, sp, magrittr, tidyverse, geosphere, flextable, officer, rgdal,
       adehabitatHR)


## Map of trips of all seals and focus with acc seals ------------------------------------------------------------------------
#Code requires access to the list of trips included in model 2, so need to wait to finalise 
gpsx <- read.csv(here::here("Datasets","pv64-2017_gps_data_with_haulout_&_trip_info.csv"), header=TRUE)
gpsx <- gpsx[,c(1:5,11,12)]
gpsx$time <- as.POSIXct(gpsx$D_DATE, tz="UTC")

#Calculate mean interbal between GPS locations
gpsx$gps_interval <- NA
for(i in 2:length(gpsx$ID)){
  if(gpsx$time[i]>gpsx$time[i-1]){
    gpsx$gps_interval[i] <- difftime(gpsx$time[i], gpsx$time[i-1], unit="mins")
  }
}
mean(gpsx$gps_interval, na.rm=TRUE)
sd(gpsx$gps_interval, na.rm=TRUE)
min(gpsx$gps_interval, na.rm=TRUE)

#Load coastline map
coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")


# A) plot with all the seals
llcord <- SpatialPoints(gpsx[,c(5,4)],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
gpsx$x <- attr(utmcoord, "coords")[,1]
gpsx$y <- attr(utmcoord, "coords")[,2]

#Specify the seals for which we have accelerometer data
IDs <- c("90","158","242","283","285")
gpsx_all <- gpsx
gpsx_seals <- gpsx[which(gpsx$ID %in% IDs),]

map_all <-  ggplot()+
  geom_path(data=gpsx_all, aes(x=x, y=y, group=ID), col="#0072B2", lwd=0.5)+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_path(data=gpsx_seals, aes(x=x, y=y, group=ID),col=c("#EDB90E"), lwd=0.5)+
  xlab("Longitude")+ylab("Latitude")+
  theme_classic()+
  #north(land, location="topright")
  scalebar(y.min =  6375392, y.max= 6581648,
           x.min = 416466.9 , x.max=530568.8,
           dist=20, dist_unit = "km",
           transform = FALSE, model = "WGS84", location = "bottomright",
           st.color="#757575", box.fill=c("#757575","white"), height = 0.009,
           box.color="#757575", st.size=3, border.size=0.1)+
  theme(legend.position="none",
        panel.grid.major = element_line(colour="transparent"),
        text = element_text(size=10))

gpsx$trip_code <- format(gpsx$trip_code, nsmall=3)
#Need access to list of trips in model 2!!!!!
#trips_analysis <- unique(df.HMM$ID)
gpsx_seals$trip_in <- ifelse(gpsx_seals$trip_code %in% trips_analysis, 1, 0)
gpsx_seals$seal_ID <- as.factor(gpsx_seals$ID)

seal_maps <- ggplot()+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_path(data=gpsx_seals, aes(x=x, y=y, group=ID, col=as.factor(trip_in)), lwd=0.5)+
  scale_color_manual(values=c("#0072B2", "red"))+
  facet_wrap(~ seal_ID, nrow=2)+
  ylab("Latitude")+
  scalebar(y.min =  6374392, y.max= 6435955,
           x.min = 416466.9 , x.max=470884.8,
           dist=10, dist_unit = "km",
           transform = FALSE, model = "WGS84", location = "bottomright",
           st.color="#757575", box.fill=c("#757575","white"), box.color="#757575",
           st.size=2,height = 0.02, st.dist=0.03, border.size=0.1)+
  theme_bw()+
  theme(legend.position="none",
        panel.grid.major = element_line(colour="transparent"),
        #text = element_text(size=7.5),
        strip.text = element_text(size=10),
        axis.text= element_text(size=6),
        axis.title = element_text(size=10))

arrange <- ggarrange(map_all, seal_maps, ncol=2, widths = c(0.3,0.7), labels=c("A", "B"))
ggsave(plot=arrange, filename=here::here("Figures","Figure 1 - seal tracks.tiff"),
       device="tiff", width = 10, height=5, units="in", dpi=300)

## Model 1 F - PR plots and parameters density distributions ---------------------------------------------------------------------------
#This code requires the momentuHMM extra functions script and for model 1 to be present in the Environment
source(here::here("R","03. Figures - MomentuHMM specific functions.R"))

# Evaluation plots model 1
tiff(here::here("Figures","Model 1 - evaluation plots.tiff"), width = 10, height = 6, units="in", res=300)
plotPR(m1)
dev.off()

#Need to debug to obtain the distribution densities
debug(momentuHMM:::plot.momentuHMM)
plot(m1)
#Debug untill line 620 then 
assign("genDensities.steps", genDensities, globalenv())
#Keep going until you plot the first graph
#then repeat until line 620
assign("genDensities.angles", genDensities, globalenv())
#now you can stop
undebug(momentuHMM:::plot.momentuHMM)

step1 <- as.data.frame(genDensities.steps[[1]])
step2 <- as.data.frame(genDensities.steps[[2]])
angle1 <- as.data.frame(genDensities.angles[[1]])
angle2 <- as.data.frame(genDensities.angles[[2]])

df.HMM <- read.csv(here::here("Output", "Model 1 - HMM dive batches classified.csv"), header=TRUE)

step <- ggplot()+
  geom_histogram(data=df.HMM, aes(x=step, y=..density..), binwidth=500 , fill="light grey",  boundary=0)+
  geom_path(data=step1, aes(x=grid , y=V2, colour="Transit"), lwd=1)+
  geom_path(data=step2, aes(x=grid, y=V2, colour="ARS"), lwd=1)+
  scale_colour_manual(name="HMM State", values=c("Transit"="#22B7F2","ARS"="#F58318"))+
  ylab("Frequency density\n")+xlab("\nStep length (m)")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="transparent"),
        panel.grid.minor = element_line(colour="transparent"))

angle <- ggplot()+
  geom_histogram(data=df.HMM, aes(x=angle, y=..density..), binwidth=0.4 , fill="lightgrey",  boundary=0)+
  geom_path(data=angle1, aes(x=grid, y=V2, colour="Transit"), lwd=1)+
  geom_path(data=angle2, aes(x=grid, y=V2, colour="ARS"), lwd=1)+
  scale_colour_manual(name="HMM State", values=c("Transit"="#22B7F2","ARS"="#F58318"))+
  ylab("Frequency density\n")+xlab("\nTurning angle (radians)")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="transparent"),
        panel.grid.minor = element_line(colour="transparent"))

m1density <- ggarrange(step, angle, common.legend=TRUE)  

ggsave(plot= m1density, filename=here::here("Figures","Model 1 - parameters density functions.tiff"), 
       device="tiff", width = 10, height=5, units="in", dpi=300)    


## Model 1 T - seals data available --------------------------------------------------------------------------------------
trips <- read.csv(here::here("Datasets", "pv64-2017_trip_summaries.csv"), header=TRUE) %>%
  mutate(Trip_Start = as.POSIXct(Trip_Start, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         Trip_End = as.POSIXct(Trip_End, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         Trip_Code = format(Trip_Code, nsmall=3),
         date = as.Date(Trip_Start))

#Filter for more than 12 hours trips
trip <- trips[-c(which(trips$PTT=="99999")),]
long_trips <- trip$Trip_Code[which(trip$Trip_Duration>12)]
trip <- trip[trip$Trip_Code %in% long_trips,]

seal.summary <- read.csv(here::here("Datasets","Seal_summary.csv"), header=TRUE) %>%
  mutate(date.cap = as.Date(Date_Captured, format="%d/%m/%Y"),
         first.week = date.cap + 7)

#Find the trips starting in the first week
trip$first.week <- NA
for(x in 1:length(trip$ID)){
  n <- which(seal.summary$Tag_Number == trip$PTT[x])
  trip$first.week[x] <- paste(seal.summary$first.week[n])
}
trip$first.week <- as.Date(trip$first.week)
trip$to.keep <- ifelse(trip$date <= trip$first.week, 0, 1)
trip <- trip[which(trip$to.keep==1),]

##Filter for round.trips at the same haul out location
trip$haulout.dist <- NA
for(x in 1:length(trip$Trip_Code)){
  trip$haulout.dist[x] <- distm(c(trip$Start_Haulout_Lon[x], trip$Start_Haulout_Lat[x]), c(trip$End_Haulout_Lon[x], trip$End_Haulout_Lat[x]), fun=distHaversine)
}
trip <- trip[which(trip$haulout.dist<=2000),]

tableS <- trip %>% group_by(ID) %>% summarise(
  `First trip start` = paste(first(Trip_Start)),
  `Last trip end`= paste(last(Trip_End)),
  `Number of foraging trips` = n(),
  mean_duration = mean(Trip_Duration),
  sd = sd(Trip_Duration)
  ) %>% mutate(
  `Mean trip duration (hours)` = paste0(format(mean_duration, digits=4)," (± ", format(sd, digits=4),")"),
  ID = as.character(ID)
  )
tableS <- tableS[,-c(5,6)]

IDs <- c("90","158","242","283","285")
tableS1 <- tableS[which(tableS$ID %in% IDs),]
tableS2 <- tableS[which(!tableS$ID %in% IDs),]

tb1 <- flextable(tableS1) %>% autofit(.) %>% align_nottext_col(., align="center")
tb2 <- flextable(tableS2) %>% autofit(.) %>% align_nottext_col(., align="center")

doc <- read_docx()
body_add_flextable(doc, tb1)
body_add_par(doc, value = "")
body_add_flextable(doc, tb2)
print(doc, target = here::here("Output", "Table - Seals data available.docx"))

## Repeatability F - Map of UD all individuals --------------------------------------------------------------------------- 
dat <- read.csv(here::here("Output", "Model 1 - HMM dive batches classified.csv"), header=TRUE)
dat <- dat[,c(2,5:15,17,18)] #for dive batches
dat$ID <- format(dat$ID, nsmall=3)
coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")
grid <- readOGR("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/MF_grid/MF_grid_1km_UTM30.shp")
x <- seq(grid@bbox[1]+250, grid@bbox[3]-250, by=500)
y <- seq(grid@bbox[2]+250, 6600733, by=500)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

ud <- dat[which(dat$state=="ARS"),]
splst <- ud %>% group_by(ID) %>% summarize(n = n()) %>% arrange(n) %>%
  filter(n > 9) %>% dplyr::select(ID) %>% extract2(1) %>% as.character()
#Remove trips that don't have eough observations
ud <- ud %>% filter(ID %in% splst) %>% droplevels()
spfor <- SpatialPointsDataFrame(coordinates(cbind(ud$x, ud$y)), data=ud, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
kd <- kernelUD(spfor[,13], h=2000, grid=xy, kern=c("bivnorm")) #Using the HMMstate as column as I want them altogether
ud.list <- list(ud95 = getverticeshr.estUD(kd[[1]], 95),
                ud75 = getverticeshr.estUD(kd[[1]], 75),
                ud50 = getverticeshr.estUD(kd[[1]], 50),
                ud25 = getverticeshr.estUD(kd[[1]], 25))

if(length(ud.list[[1]]@polygons[[1]]@Polygons)==1){
  udsp <- ud.list[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
    as.data.frame(.) %>%
    mutate(Id = rep("95", length(V1)), poly = rep(1, length(V1)))
} else {
  n <- length(ud.list[[1]]@polygons[[1]]@Polygons)
  udsp <- ud.list[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
    as.data.frame(.) %>%
    mutate(Id = rep("95", length(V1)), poly = rep(1, length(V1)))
  for(i in 2:n){
    udsp1 <- ud.list[[1]]@polygons[[1]]@Polygons[[n]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep("95", length(V1)), poly = rep(n, length(V1)))
    udsp <- rbind(udsp, udsp1)
  }
}
uds <- c(95, 75,50,25)
for(x in 2:length(uds)){
  if(length(ud.list[[x]]@polygons[[1]]@Polygons)==1){
    udtmp <- ud.list[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(1, length(V1)))
  }  else {
    n <- length(ud.list[[x]]@polygons[[1]]@Polygons)
    udtmp <- ud.list[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(1, length(V1)))
    for(i in 2:n){
      udtmp1 <- ud.list[[x]]@polygons[[1]]@Polygons[[i]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(i, length(V1)))
      udtmp <- rbind(udtmp, udtmp1)
    }
  }
  udsp <- rbind(udsp, udtmp)
}
colnames(udsp) <- c("x","y","Id","poly")
udsp$group <- paste(udsp$Id, udsp$poly)

UD_map <- ggplot()+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_polygon(data= udsp[which(udsp$Id==95),], aes(x=x,y=y, group=as.factor(group), fill= "95%"), alpha=0.5)+
  geom_polygon(data= udsp[which(udsp$Id==75),], aes(x=x,y=y, group=as.factor(group), fill= "75%"), alpha=0.5)+
  geom_polygon(data= udsp[which(udsp$Id==50),], aes(x=x,y=y, group=as.factor(group), fill= "50%"), alpha=0.5)+
  geom_polygon(data= udsp[which(udsp$Id==25),], aes(x=x,y=y, group=as.factor(group), fill= "25%"), alpha=0.5)+
  scale_fill_manual(values=c(c("25%"="#F26E2C", "50%"="#F5BE27", "75%"="#F5EE27","95%"="#91CF14")), name="Utilization distribution")+
  scale_color_manual(values=c("#000000"), lab=c("ARS"), name="HMM state")+
  xlim(412635.3,529903.1)+ylim(6373710, 6501614)+
  xlab("Longitude")+ylab("Latitude")+
  theme_classic()+
  theme(legend.position = "top")
ggsave(plot= UD_map, filename=here::here("Figures","Repeatability - IMF UDs.tiff"), 
       device="tiff", width = 10, height=5, units="in", dpi=300)   

# Repetability F - Comparison of consecutive trip -----------------------------------------------------------
#Read in dataset and other useful things
dat <- read.csv(here::here("Output", "Model 1 - HMM dive batches classified.csv"), header=TRUE)
dat <- dat[,c(2,5:15,17,18)] #for dive batches
dat$ID <- format(dat$ID, nsmall=3)
coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")
grid <- readOGR("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/MF_grid/MF_grid_1km_UTM30.shp")
x <- seq(grid@bbox[1]+250, grid@bbox[3]-250, by=500)
y <- seq(grid@bbox[2]+250, 6600733, by=500)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

##Map an individual trip utilization distributions
trip1 <- dat[which(dat$ID=="14468.027"),]
trip2 <- dat[which(dat$ID=="14468.031"),]

xmax <- max(max(trip1$x), max(trip2$x))
xmin <- min(min(trip1$x), min(trip2$x))
ymax <- max(max(trip1$y), max(trip2$y))
ymin <- min(min(trip1$y), min(trip2$y))

UD_and_ARS_map1 <- {
  trip1ud <- trip1[which(trip1$state=="ARS"),]
  spfor <- SpatialPointsDataFrame(coordinates(cbind(trip1ud$x, trip1ud$y)),  data=trip1ud, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
  kd <- kernelUD(spfor[,1], h=2000, grid=xy, kern=c("bivnorm"))
  ud.list <- list(ud95 = getverticeshr.estUD(kd[[1]], 95),
                  ud75 = getverticeshr.estUD(kd[[1]], 75),
                  ud50 = getverticeshr.estUD(kd[[1]], 50),
                  ud25 = getverticeshr.estUD(kd[[1]], 25))
  
  if(length(ud.list[[1]]@polygons[[1]]@Polygons)==1){
    udsp <- ud.list[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep("95", length(V1)), poly = rep(1, length(V1)))
  } else {
    n <- length(ud.list[[1]]@polygons[[1]]@Polygons)
    udsp <- ud.list[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep("95", length(V1)), poly = rep(1, length(V1)))
    for(i in 2:n){
      udsp1 <- ud.list[[1]]@polygons[[1]]@Polygons[[n]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep("95", length(V1)), poly = rep(n, length(V1)))
      udsp <- rbind(udsp, udsp1)
    }
  }
  uds <- c(95, 75,50,25)
  for(x in 2:length(uds)){
    if(length(ud.list[[x]]@polygons[[1]]@Polygons)==1){
      udtmp <- ud.list[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(1, length(V1)))
    }  else {
      n <- length(ud.list[[x]]@polygons[[1]]@Polygons)
      udtmp <- ud.list[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(1, length(V1)))
      for(i in 2:n){
        udtmp1 <- ud.list[[x]]@polygons[[1]]@Polygons[[i]]@coords %>% 
          as.data.frame(.) %>%
          mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(i, length(V1)))
        udtmp <- rbind(udtmp, udtmp1)
      }
    }
    udsp <- rbind(udsp, udtmp)
  }
  colnames(udsp) <- c("x","y","Id","poly")
  udsp$group <- paste(udsp$Id, udsp$poly)
  
  ggplot()+
    annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
    geom_polygon(data= udsp[which(udsp$Id==95),], aes(x=x,y=y, group=as.factor(group), fill= "95%"), alpha=0.5)+
    geom_polygon(data= udsp[which(udsp$Id==75),], aes(x=x,y=y, group=as.factor(group), fill= "75%"), alpha=0.5)+
    geom_polygon(data= udsp[which(udsp$Id==50),], aes(x=x,y=y, group=as.factor(group), fill= "50%"), alpha=0.5)+
    geom_polygon(data= udsp[which(udsp$Id==25),], aes(x=x,y=y, group=as.factor(group), fill= "25%"), alpha=0.5)+
    scale_fill_manual(values=c(c("25%"="#F26E2C", "50%"="#F5BE27", "75%"="#F5EE27","95%"="#91CF14")), name="Utilization\ndistribution")+
    #geom_point(data=tmp, aes(x=x,y=y,col=as.factor(HMMstate)))+
    #scale_color_manual(values=c("#2EB3F5", "#000000", "#FFFFFF00"), lab=c("Transit", "ARS", " "), name="HMM state")+
    geom_path(data=trip1, aes(x=x,y=y), lty="dashed")+
    geom_point(data=trip1[which(trip1$state=="ARS"),], aes(x=x,y=y, col="ARS"))+
    scale_color_manual(values=c("#000000"), lab=c("ARS"), name="HMM state")+
    xlim(xmin-5000, xmax+5000)+ylim(ymin-5000, ymax+5000)+
    xlab("Longitude")+ylab("Latitude")+
    theme_classic()
}
UD_and_ARS_map1
UD_and_ARS_map2 <- {
  trip2ud <- trip2[which(trip2$state=="ARS"),]
  spfor <- SpatialPointsDataFrame(coordinates(cbind(trip2ud$x, trip2ud$y)),  data=trip2ud, proj4string = CRS("+proj=utm +zone=30 ellps=WGS84 +datum=WGS84"))
  kd <- kernelUD(spfor[,1], h=2000, grid=xy, kern=c("bivnorm"))
  ud.list <- list(ud95 = getverticeshr.estUD(kd[[1]], 95),
                  ud75 = getverticeshr.estUD(kd[[1]], 75),
                  ud50 = getverticeshr.estUD(kd[[1]], 50),
                  ud25 = getverticeshr.estUD(kd[[1]], 25))
  
  if(length(ud.list[[1]]@polygons[[1]]@Polygons)==1){
    udsp <- ud.list[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep("95", length(V1)), poly = rep(1, length(V1)))
  } else {
    n <- length(ud.list[[1]]@polygons[[1]]@Polygons)
    udsp <- ud.list[[1]]@polygons[[1]]@Polygons[[1]]@coords %>% 
      as.data.frame(.) %>%
      mutate(Id = rep("95", length(V1)), poly = rep(1, length(V1)))
    for(i in 2:n){
      udsp1 <- ud.list[[1]]@polygons[[1]]@Polygons[[n]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep("95", length(V1)), poly = rep(n, length(V1)))
      udsp <- rbind(udsp, udsp1)
    }
  }
  uds <- c(95, 75,50,25)
  for(x in 2:length(uds)){
    if(length(ud.list[[x]]@polygons[[1]]@Polygons)==1){
      udtmp <- ud.list[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(1, length(V1)))
    }  else {
      n <- length(ud.list[[x]]@polygons[[1]]@Polygons)
      udtmp <- ud.list[[x]]@polygons[[1]]@Polygons[[1]]@coords %>% 
        as.data.frame(.) %>%
        mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(1, length(V1)))
      for(i in 2:n){
        udtmp1 <- ud.list[[x]]@polygons[[1]]@Polygons[[i]]@coords %>% 
          as.data.frame(.) %>%
          mutate(Id = rep(as.character(uds[x]), length(V1)), poly = rep(i, length(V1)))
        udtmp <- rbind(udtmp, udtmp1)
      }
    }
    udsp <- rbind(udsp, udtmp)
  }
  colnames(udsp) <- c("x","y","Id","poly")
  udsp$group <- paste(udsp$Id, udsp$poly)
  
  ggplot()+
    annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
    geom_polygon(data= udsp[which(udsp$Id==95),], aes(x=x,y=y, group=as.factor(group), fill= "95%"), alpha=0.5)+
    geom_polygon(data= udsp[which(udsp$Id==75),], aes(x=x,y=y, group=as.factor(group), fill= "75%"), alpha=0.5)+
    geom_polygon(data= udsp[which(udsp$Id==50),], aes(x=x,y=y, group=as.factor(group), fill= "50%"), alpha=0.5)+
    geom_polygon(data= udsp[which(udsp$Id==25),], aes(x=x,y=y, group=as.factor(group), fill= "25%"), alpha=0.5)+
    scale_fill_manual(values=c(c("25%"="#F26E2C", "50%"="#F5BE27", "75%"="#F5EE27","95%"="#91CF14")), name="Utilization\ndistribution")+
    #geom_point(data=tmp, aes(x=x,y=y,col=as.factor(HMMstate)))+
    #scale_color_manual(values=c("#2EB3F5", "#000000", "#FFFFFF00"), lab=c("Transit", "ARS", " "), name="HMM state")+
    geom_path(data=trip2, aes(x=x,y=y), lty="dashed")+
    geom_point(data=trip2[which(trip2$state=="ARS"),], aes(x=x,y=y, col="ARS"))+
    scale_color_manual(values=c("#000000"), lab=c("ARS"), name="HMM state")+
    xlim(xmin-5000, xmax+5000)+ylim(ymin-5000, ymax+5000)+
    xlab("Longitude")+ylab("Latitude")+
    theme_classic()
}
UD_and_ARS_map2
ggarrange(UD_and_ARS_map1, UD_and_ARS_map2, ncol = 1,common.legend = TRUE)

overlap.trips <- data.frame(Trips = "Bhattacharyya’s affinity",
                            UD95 = c(0.88),
                            UD75 = c(0.67),
                            UD50 = c(0.43),
                            UD25 = c(0.19))
colnames(overlap.trips) <- c("Trips overlap","95% UD","75% UD","50% UD","25% UD")
tbf <- ggtexttable(overlap.trips, rows=NULL, theme = ttheme(
  colnames.style = colnames_style(fill="white", linecolor = "black"),
  rownames.style = rownames_style(color = "white"),
  tbody.style = tbody_style(fill=c("white", "white"), linecolor = "black"))) 
tbf <- table_cell_font(tbf, row =2, column = 1, face="bold", size=11)

comparison_arrange <- ggarrange(UD_and_ARS_map1, UD_and_ARS_map2, tbf, ncol = 1,common.legend = TRUE, heights = c(1,1,0.3))
ggsave(plot= comparison_arrange, filename=here::here("Figures","Repeatability - Consecutive trips comparison.tiff"), 
       device="tiff", width = 6, height=8, units="in", dpi=300)   

## Repeatability T - Table summary data -------------------------------------------------------
path <- "C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal repeatability/Output/ARS trips overlap within individuals/"
files <- paste0(path,"ARS_95UD","/",list.files(paste0(path,"ARS_95UD")))
tables <- lapply(files, read.csv)
UD95 <- do.call(rbind, tables) %>% mutate(UD= rep("UD 95", length(seal_ID)))
files <- paste0(path,"ARS_75UD","/",list.files(paste0(path,"ARS_75UD")))
tables <- lapply(files, read.csv)
UD75 <- do.call(rbind, tables) %>% mutate(UD= rep("UD 75", length(seal_ID)))
files <- paste0(path,"ARS_50UD","/",list.files(paste0(path,"ARS_50UD")))
tables <- lapply(files, read.csv)
UD50 <- do.call(rbind, tables) %>% mutate(UD= rep("UD 50", length(seal_ID)))
files <- paste0(path,"ARS_25UD","/",list.files(paste0(path,"ARS_25UD")))
tables <- lapply(files, read.csv)
UD25 <- do.call(rbind, tables) %>% mutate(UD= rep("UD 25", length(seal_ID)))

UD <- rbind(UD95,UD75,UD50,UD25)

#Summary table of individuals repeatability 
summary <- UD %>% group_by(seal_ID, UD) %>% summarise(
  N_trips = n(),
  overlap_mean = mean(overlap),
  overlap_sd = sd(overlap)
)
tmp_mean <- summary[,c(4,2)] %>% unstack()
tmp_sd <- summary[,c(5,2)] %>% unstack()
summary2 <- summary %>% group_by(seal_ID) %>% summarise(N_trips = unique(N_trips)) %>% 
  mutate(mean.95 = tmp_mean[,4], sd.95 = tmp_sd[,4],
         mean.75 = tmp_mean[,3], sd.75 = tmp_sd[,3],
         mean.50 = tmp_mean[,2], sd.50 = tmp_sd[,2],
         mean.25 = tmp_mean[,1], sd.25 = tmp_sd[,1]) %>%
  mutate(overlap95 = paste0(format(mean.95, digits=2)," (± ", format(sd.95, digits=2),")"),
         overlap75 = paste0(format(mean.75, digits=2)," (± ", format(sd.75, digits=2),")"),
         overlap50 = paste0(format(mean.50, digits=2)," (± ", format(sd.50, digits=2),")"),
         overlap25 = paste0(format(mean.25, digits=2)," (± ", format(sd.25, digits=2),")"))
summary2 <- summary2[,c(1,2,11:14)]
colnames(summary2) <- c("ID","Number of foraging trips","overlap.95","overlap.75",
                        "overlap.50", "overlap.25")

ft <- flextable(summary2) 
typology <- data.frame(
  col_keys = c("ID","Number of foraging trips","overlap.95","overlap.75","overlap.50", "overlap.25"),
  what =c("ID","Number of foraging trips","Consecutive trips overlap","Consecutive trips overlap","Consecutive trips overlap","Consecutive trips overlap"),
  measure=c("ID","Number of foraging trips", "95% UD", "75% UD", "50% UD", "25% UD")
)
ft <- set_header_df(ft, mapping = typology, key="col_keys") %>%
 merge_h(., part = "header") %>% merge_v(., part = "header") %>%
 theme_booktabs(.) %>% autofit(.) %>% fix_border_issues(.) %>%
 align_nottext_col(., align="center") %>% align_text_col(., align="center")
ft

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = here::here("Output", "Table - UD analysis overlap.docx"))

## Accelerometry - Foraging tactics --------------------------------------------------------------------
path <- "C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/Processed accelerometer data/"
files <- paste0(path,list.files(paste0(path)))
tables <- lapply(files, read.csv)
dat <- do.call(rbind, tables) 

dat <- dat [,c(1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,25,26,27,33,34,37,47)]
sealdf$dive_type <- NA
n <- 1
for(n in 1:length(sealdf$ID)){
  if(is.na(sealdf$B_PrCA_arch[n])){
    sealdf$dive_type[n] <- NA
  } else {
    if(sealdf$B_PrCA_arch[n]>0 & sealdf$B_pitch.diff20[n]==0){sealdf$dive_type[n] <- "1"} #Only prca
    if(sealdf$B_PrCA_arch[n]==0 & sealdf$B_pitch.diff20[n]>0){sealdf$dive_type[n] <- "2"} #only benthic
    if(sealdf$B_PrCA_arch[n]>0 & sealdf$B_pitch.diff20[n]>0){sealdf$dive_type[n] <- "3"} #both
    if(sealdf$B_PrCA_arch[n]==0 & sealdf$B_pitch.diff20[n]==0){sealdf$dive_type[n] <- "4"}
  }
  print(n)
} 

graph.plot <- sealdf %>% group_by(seal_ID, dive_type) %>% summarise (count=n())
graph.plot <- graph.plot[which(graph.plot$dive_type>0),]
tot.dives <- graph.plot %>% group_by(seal_ID) %>% summarise(n.dives = sum(count))
graph.plot$tot.dives <- rep(tot.dives$n.dives, each=4)
graph.plot$p.dives <- graph.plot$count/graph.plot$tot.dives

ggplot(graph.plot, aes(x=as.factor(seal_ID), y=p.dives, fill=dive_type))+
  geom_bar(position="dodge", stat = "identity")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                    name= "Dive type",
                    labels=c("Peaks in dynamic acceleration", "Changes in head pitch angle", "Both catch strategies", "No attempts"))+
  ylab("Proportion of dives")+
  xlab("Seal ID")+
  theme_classic()
theme(legend.position="bottom", legend.box="vertical")



## Model 2 - Panel covariates and output ---------------------------------------------------------------
#Memory panel
pixel <- as(rr, "SpatialPixelsDataFrame")
pixel_df <- as.data.frame(pixel)
colnames(pixel_df) <- c("value", "x", "y")
pixelvalues <- pixel_df[which(pixel_df$value!=0),]

ggplot()+
  geom_tile(data=pixel_df, aes(x=x, y=y, fill=value))+
  annotation_spatial(land, fill = "lightgrey", lwd = 0)+
  scale_fill_distiller(name="Proportion of dives\n spent searching", palette="YlGn", trans="reverse")+
  xlab("Longitude")+ylab("Latitude")+
  theme_bw()+
  guides(fill=guide_colourbar(title.position = "right", title.vjust=1))+
  theme(legend.position="top",
        panel.grid.major = element_line(colour="transparent"), text=element_text(family="serif"))
