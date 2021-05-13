### Figures and tables
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose: Create all the figures and tables in Iorio-Merlo et al. 

# Created on: 12/01/2021
# Updated on: 13/05/2021

## Load packages --------------------------------------------------------------------------------------
library(pacman)
p_load(ggplot2, ggspatial, ggsn, ggpubr, sf, sp, magrittr, tidyverse, geosphere, flextable, officer, rgdal,
       adehabitatHR, REdaS, roll, mefa, lubridate, raster)

## Table 1 - Summary table----------------------------------------------------------------------------------
ids <- c(242,158,90,285,283)

seal.info <- read.csv(here::here("Datasets","Seal_summary.csv"), header=TRUE) 
seal.info <- seal.info[which(seal.info$ID..Purple...recaptured.another.year. %in% ids & seal.info$Year.Captured==2017),]
seal.info <- seal.info[,c(2,5,10,11,12,14,20)]
colnames(seal.info) <- c("Tag_start","ID","Sex","Weight","Length","Tag_end","Pregnancy status")

#Trips in model 1
# m1data <- read.csv(here::here("Output","Model 1 - HMM dive batches classified.csv"), header=TRUE) %>%
#   mutate(ID = format(ID, nsmal=3))
# m1data <- m1data[which(m1data$seal_ID %in% ids),]
# m1trips <- m1data %>% group_by(seal_ID) %>% summarise(m1_trips = length(unique(ID)))
# seal.info$m1_trips <- m1trips$m1_trips

#Manually insert the UD overlap from the summaty table you have already created
# seal.info$UD_mean_over <- c("0.59 (± 0.330)", "0.69 (± 0.303)","0.76 (± 0.233)","0.84 (± 0.159)","0.71 (± 0.340)")

#Spatial memory
mem <- read.csv(here::here("Output","Spatial memory - dataset for rasters.csv")) %>%
  mutate(ID = format(ID, nsmall=3))
mem.trips <- mem %>% group_by(seal_ID) %>% summarise(
  start_date = first(as.Date(start.time, format="%d/%m/%Y")),
  mem_n_trips = length(unique(ID)))
seal.info$mem.start <- mem.trips$start_date
seal.info$mem_n_trips <- mem.trips$mem_n_trips

#Prey encounters
path <- here::here("Output","Processed accelerometer parameters")
files <- paste0(path,"/",list.files(paste0(path)))
tables <- lapply(files, read.csv)
seal <- do.call(rbind, tables) 
seal <- seal[,c(1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,25,26,27,33,34,37,47)] %>%
  mutate(seal_ID = ID,
         trip_code = format(Trip_Code, nsmall=3),
         ID = trip_code,
         time = as.POSIXct(DS_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         time.end = as.POSIXct(DE_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"))
acc.summ <- seal %>% group_by(seal_ID) %>% summarise(
  n.dives.pe = length(which(B_PrCA_arch>0 | B_pitch.diff20>0)),
  Cox_met = sum(B_PrCA_arch, na.rm=TRUE),
  Brass_met = sum(B_pitch.diff20, na.rm=TRUE),
  tot.attempts = (sum(B_PrCA_arch, na.rm=TRUE) + sum(B_pitch.diff20, na.rm=TRUE)) - sum(B_overlap_A, na.rm=TRUE)
  ) %>% mutate (
    Cox_percent = Cox_met/tot.attempts*100,
    Brass_percent = Brass_met/tot.attempts*100
  )
seal.info$n.dive.pe <- acc.summ$n.dives.pe
seal.info$Cox_met <- paste0(acc.summ$Cox_met, " (",round(acc.summ$Cox_percent),"%)")
seal.info$Brass_met <- paste0(acc.summ$Brass_met, " (",round(acc.summ$Brass_percent),"%)")
seal.info$tot.attempts <- acc.summ$tot.attempts

#Trips in model 2
m2data <- read.csv(here::here("Output","Model 2 - HMM dive batches classified.csv"), header=TRUE) %>%
  mutate(ID = format(ID, nsmal=3))
m2data <- m2data[which(m2data$seal_ID %in% ids),]
m2trips <- m2data %>% group_by(seal_ID) %>% summarise(
  m2_start = first(as.Date(start.time)),
  m2_end = last(as.Date(end.time)),
  m2_trips = length(unique(ID)))
seal.info$m2_start <- m2trips$m2_start
seal.info$m2_end <- m2trips$m2_end
seal.info$m2_trips <- m2trips$m2_trips

info.short <- seal.info[,c(2:5,7,8,9,10,13:16)]
colnames(info.short) <- c("Seal ID", "Sex", "Weight (Kg)","Length (cm)", "Pregnancy status",
                          "Start date for memory map","# trips in memory map","Dives with prey encounter",
                          "Total prey encounters", "Data start date","Data end date","# trips in model 2")
ft <- flextable(info.short)  %>% autofit(.) %>% align_nottext_col(., align="center")
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = here::here("Output", "Table - Table 1 data and res summary.docx"))

info.ext <- seal.info[,c(2,3:5,1,6:17)]
colnames(info.ext) <- c("Seal ID", "Sex", "Weight (Kg)","Length (cm)", 
                          "Tag start date", "Tag end date", "# trips in Model 1", "Mean 95% UD overlap",
                          "Trips for memory map start","# trips in memory map",
                          "Dives with prey encounter","Encounters by Cox et al. method",
                          "Encounters by Brasseur et al. method", "Total prey encounters",
                          "Trips with accelerometer data start", "Trips with accelerometer data end", 
                          "# trips in model 2")
ft2 <- flextable(info.ext)  %>% autofit(.) %>% align_nottext_col(., align="center")
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft2)
print(doc, target = here::here("Output", "Table - Table S1 data and res extended.docx"))


## Table 2 - Model selection table ---------------------------------------------------------------------------
selection <- read.csv(here::here("Output","Model 2 - Covariates model selection.csv"), header=TRUE)

colnames(selection) <- c(" ","Log-Likelihood","AIC","BIC","Delta AIC","Delta BIC")
ft <- flextable(selection)  %>% autofit(.) %>% align_nottext_col(., align="center")

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = here::here("Output", "Table - Model 2 covariates selection.docx"))

## Figure 1 - Map of trips of all seals and focus with acc seals ------------------------------------------------------------------------
#Code requires access to the list of trips included in model 2, so need to wait to finalise 
gpsx <- read.csv(here::here("Datasets","pv64-2017_gps_data_with_haulout_&_trip_info.csv"), header=TRUE)
gpsx <- gpsx[,c(1:5,11,12)]
gpsx$time <- as.POSIXct(gpsx$D_DATE, tz="UTC")

seal.HMM <- read.csv(here::here("Output","Model 2 - HMM dive batches classified.csv"), header=TRUE) %>%
  mutate(ID= format(ID, nsmall=3))

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

c("#2998FF")

map_all <-  ggplot()+
  geom_path(data=gpsx_all, aes(x=x, y=y, group=ID), col="#2998FF", lwd=0.5)+
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
        text = element_text(size=12))

gpsx$trip_code <- format(gpsx$trip_code, nsmall=3)
#Need access to list of trips in model 2!!!!!
trips_analysis <- unique(seal.HMM$ID)
gpsx_seals$trip_in <- ifelse(gpsx_seals$trip_code %in% trips_analysis, 1, 0)
gpsx_seals$seal_ID <- as.factor(gpsx_seals$ID)

seal_maps <- ggplot()+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_path(data=gpsx_seals[which(gpsx_seals$trip_in==0),], aes(x=x, y=y, group=ID), lwd=0.5, col="#2998FF")+
  geom_path(data=gpsx_seals[which(gpsx_seals$trip_in==1),], aes(x=x, y=y, group=ID), lwd=0.5, col="red")+
  #scale_color_manual(values=c("#0072B2", "red"))
  facet_wrap(~ seal_ID, nrow=2)+
  xlab("Longitude")+ylab("Latitude")+
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
        strip.text = element_text(size=12),
        axis.text= element_text(size=6),
        axis.title = element_text(size=12))

arrange <- ggarrange(map_all, seal_maps, ncol=2, widths = c(0.3,0.7), labels=c("A", "B"))
ggsave(plot=arrange, filename=here::here("Figures","Figure 1 - seal tracks.tiff"),
       device="tiff", width = 10, height=5, units="in", dpi=300)

## Figure 2 - Repeatability: Map of mean ARS locations and density curve of observed and null distribution of overlap------------------------------------------------------------
#Density curve
obs <- read.csv(here::here("Output", "Repeatability - overlap between April and May.csv"), header=TRUE)
obs <- obs %>%  mutate(curve = rep("obs", length(obs$seal_ID)))
null <- read.csv(here::here("Output", "Repeatability - null distribution overlap.csv"), header=TRUE)
null <- null %>% mutate(curve = rep("null", length(null$seala)))

dat <- data.frame(curve= c(obs$curve, null$curve), overlap = c(obs$BA_overlap, null$BA_overlap)) 

density.curve <- ggplot(dat, aes(x=overlap, fill=curve, col=curve))+
  geom_density(mapping = (aes(y = ..scaled..)),alpha=0.5)+
  scale_color_manual(values = c("#D6D6D6", "#69D66E"), labels=c("null","Observed"), name=" ")+
  scale_fill_manual(values = c("#D6D6D6", "#69D66E"),labels=c("null","Observed"), name=" ")+
  xlim(0,1)+
  ylab("Frequency density")+xlab("Bhattacharyya's affinity")+
  theme_classic()+
  theme(legend.position = "right",
        plot.margin = margin(0,12,0,12,"cm"),
        text=element_text(size=25))

#Map 
dat <- read.csv(here::here("Figures","Data for figures","Repeatability map data.csv"), header=TRUE)

df <- data.frame(seal.ID = unique(dat$seal_ID), lon1 = NA, lat1 = NA, lon2 = NA, lat2 = NA)
df <- df[-which(df$seal.ID==59 | df$seal.ID==280),]
months <- c(1,2)
for(i in 1:length(df$seal.ID)){
  for(y in 1:length(months)){
    tmp <- dat[which(dat$seal_ID==df$seal.ID[i] & dat$month == months[y]),]
    tmp2 <- data.frame(cord.lat = tmp$batch.start.lat, cord.long = tmp$batch.start.lon, lat = NA, long = NA,
                       a= NA, b= NA, c= NA)
    X = 0.0
    Y = 0.0
    Z = 0.0
    for(x in 1:length(tmp2$cord.lat)){
      tmp2$lat[x] <- tmp2$cord.lat[x] *pi/180
      tmp2$long[x] <- tmp2$cord.long[x]* pi/180
      
      tmp2$a[x] <- cos(tmp2$lat[x])*cos(tmp2$long[x])
      tmp2$b[x] <- cos(tmp2$lat[x])*sin(tmp2$long[x])
      tmp2$c[x] <- sin(tmp2$lat[x])
      
      X = X + tmp2$a[x]
      Y = Y + tmp2$b[x]
      Z = Z + tmp2$c[x] 
    }
    X = X/length(tmp2$cord.lat)
    Y = Y/length(tmp2$cord.lat)
    Z = Z/length(tmp2$cord.lat)
    
    lon = atan2(Y, X)*180/pi
    hyp = sqrt(X*X+Y*Y)
    lat = atan2(Z, hyp)*180/pi
    
    if(y == 1){
      df$lon1[i] <- lon
      df$lat1[i] <- lat
    }
    if(y == 2){
      df$lon2[i] <- lon
      df$lat2[i] <- lat
    } 
  }
}

llcord <- SpatialPoints(df[,c(2,3)],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
df$x1 <- attr(utmcoord, "coords")[,1]
df$y1 <- attr(utmcoord, "coords")[,2]

llcord <- SpatialPoints(df[,c(4,5)],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
df$x2 <- attr(utmcoord, "coords")[,1]
df$y2 <- attr(utmcoord, "coords")[,2]

coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")

c("#E06636", "#F2DF53")
color1 = "#E06636"
color2 = "#F2DF53"

mean.loc.map <- ggplot()+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_segment(data=df,aes(x=x1,y=y1,xend=x2,yend=y2))+
  geom_point(data=df, aes(x=x1, y=y1, fill="month1"),shape=21, size=5)+
  geom_point(data=df, aes(x=x2, y=y2, fill="month2"), shape=21, size=5)+
  xlab("Longitude")+ylab("Latitude")+
  xlim(418635.3,489903.1)+ylim(6383710, 6451614)+  
  scale_fill_manual(values=c("month1" = color1,"month2" = color2), labels = c("April","May"), name=" ")+ 
  theme_classic()+
  theme(legend.position = "right",
        plot.margin = margin(0,0,0,0,"cm"),
        text=element_text(size=25))

ggarrange(mean.loc.map, density.curve, ncol=1, heights=c(1,0.3), align="v", labels = c("A", "B"), font.label = list(size=25))
ggsave(here::here("Figures","Repeatability - map and density distribution.tiff"),
       device="tiff", height = 15, width = 23, unit="in")

## Figure 3 - Panel covariates and output ---------------------------------------------------------------
#Memory panel
seal.HMM <- read.csv(here::here("Output","Model 2 - HMM dive batches classified.csv"), header=TRUE)
seal <- seal.HMM[which(seal.HMM$seal_ID == 242),]

coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")

raster_path <- here::here("Output","Memory_grids")
all_rasters <- list.files(raster_path,
                          full.names = TRUE,
                          pattern = ".tif$")
stack <- stack(all_rasters)
mem.raster <- stack$seals_memory_grids_1

pixel <- as(mem.raster, "SpatialPixelsDataFrame")
pixel_df <- as.data.frame(pixel)
colnames(pixel_df) <- c("value", "x", "y")
pixelvalues <- pixel_df[which(pixel_df$value!=0),]

xmin <- min(seal$x)-2000
xmax <- max(seal$x)+2000
ymin <- min(seal$y)-2000
ymax <- max(seal$y)+1000

memory.panel <- ggplot()+
  geom_tile(data=pixel_df, aes(x=x, y=y, fill=value))+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  scale_fill_distiller(name="Proportion of dive \nbatches spent searching", palette="YlGn", trans="reverse")+
  xlim(xmin, xmax)+ ylim(ymin, ymax)+
  xlab("Longitude")+ylab("Latitude")+
  theme_bw()+
  guides(fill=guide_colourbar(title.position = "right", title.vjust=1))+
  theme(legend.position="top",
        panel.grid.major = element_line(colour="transparent"))
#text=element_text(family="serif") for times new roman
#Attempts track
library(viridis)
foraging.graph.df <- seal[which(seal$HMMstate>0),]

foraging.panel <- ggplot()+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_point(data=foraging.graph.df, aes(x=x, y=y, col=batch.foraging.index), size=1.5)+
  scale_color_viridis(name="Mean batch prey encounters",trans="reverse")+
  xlim(xmin, xmax)+ ylim(ymin, ymax)+
  xlab("Longitude")+ylab("Latitude")+
  theme_bw()+
  guides(col=guide_colourbar(title.position = "right", title.vjust=0.8))+
  theme(legend.position="top",
        panel.grid.major = element_line(colour="transparent")) 

#HMM track panel
HMM.panel <- ggplot()+ 
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_point(data=seal, aes(x=x, y=y, col=as.factor(HMMstate), group=ID), size=1.5)+
  geom_path(data=seal, aes(x=x, y=y, col=as.factor(HMMstate), group=ID), lwd=0.5)+
  scale_color_manual(values=c("#06368F","#1EB4E6","#B3B1B1E2"), name="HMM state",
                     label=c("ARS", "Transit", " "))+
  xlim(xmin, xmax)+ ylim(ymin, ymax)+
  xlab("Longitude")+ylab("Latitude")+
  theme_bw()+
  guides(fill=guide_colourbar(title.position = "right", title.vjust=1))+
  theme(legend.position="top",
        panel.grid.major = element_line(colour="transparent"))
c("#1EB4E6")

panel <- ggarrange(foraging.panel,memory.panel, HMM.panel, ncol=1, nrow=3, labels = c("A","B","C"))

ggsave(plot = panel, filename=here::here("Figures","Model 2 - covariates example seal 242.tiff"), 
       device="tiff", width = 130, height=260, units="mm", dpi=300)    

## Figure 4 - Covairates stationary and transition probability -----------------------------------
#This code sources another script where there are copied functions from the momentuHMM pakcage (McClintock & Michelot 2018)
#These were required to calculate the predicted curves and the confidence intervals.
#McClintock, B. T. and T. Michelot (2018). "momentuHMM : R package for generalized hidden Markov models of animal movement." Methods in Ecology and Evolution 9(6): 1518-1530.

#This code requires the m2 model to be in the environment
source(here::here("R","Figures - MomentuHMM specific functions.R"))
seal.HMM <- read.csv(here::here("Output","Model 2 - HMM dive batches classified.csv"), header=TRUE)

{#Run to prepare the data and all the objects needed for plots
  m <- m2
  data <- m$data
  beta <- m$mle$beta
  nbStates <- length(m$stateNames)
  ref<- c(1,2)
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$data,par=m$mle)
  recharge <- reForm$recharge
  hierRecharge <- reForm$hierRecharge
  newformula <- reForm$newformula
  nbCovs <- reForm$nbCovs
  aInd <- reForm$aInd
  nbG0covs <- reForm$nbG0covs
  nbRecovs <- reForm$nbRecovs
  g0covs <- reForm$g0covs
  recovs <- reForm$recovs
  mixture <- m$conditions$mixtures
  alpha <- 0.95 
  Sigma <- m$mod$Sigma
  rawCovs <- m$rawCovs
  covs <-  data.frame(batch.foraging.index=mean(seal.HMM$batch.foraging.index, na.rm=TRUE), 
                      memory=mean(seal.HMM$memory, na.rm=TRUE)) 
  #Create gamIndex
  gamInd<-(length(m$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixture+1):(length(m$mod$estimate))-(ncol(m$covsPi)*(mixture-1))-ifelse(nbRecovs,nbRecovs+1+nbG0covs+1,0)-ncol(m$covsDelta)*(nbStates-1)*(!m$conditions$stationary)*mixture
  quantSup<-qnorm(1-(1-alpha)/2)
  
  # Stationary probability foraging index
  cov = 1
  othercova = 2
  gridLength <- 101
  hGridLength <- gridLength
  inf <- min(rawCovs[,cov],na.rm=TRUE)
  sup <- max(rawCovs[,cov],na.rm=TRUE)
  tempCovs <- NA
  tempCovs <- data.frame(
    cov1 =  rep(seq(inf,sup,length=gridLength),each=hGridLength/gridLength),
    cov2  = rep(covs[1,othercova], gridLength))
  colnames(tempCovs) <-
    c(names(m$rawCovs[cov]),names(m$rawCovs[othercova]))
  tmpcovs<-covs[names(rawCovs)]
  splineList=NULL
  tmpSplineInputs<-getSplineFormula(newformula,data,tempCovs)
  desMat <- model.matrix(tmpSplineInputs$formula, data=tmpSplineInputs$covs)
  probs <- stationary(m, covs=desMat)  
  probs2 <- as.data.frame(probs)
  
  lci <- matrix(NA,gridLength,length(ref))
  uci <- matrix(NA,gridLength,length(ref))
  
  state = c(1,2)
  for(state in 1:length(state)){
    dN <- t(apply(desMat, 1, function(x)
      numDeriv::grad(get_stat,m$mod$estimate[gamInd[unique(c(m$conditions$betaCons))]],covs=matrix(x,1),nbStates=nbStates,i=state,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=m$conditions$workBounds$beta,mixture=mixture,ref=ref)))
    tmpSig <- Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]
    se <- t(apply(dN, 1, function(x)
      suppressWarnings(sqrt(x%*%tmpSig%*%x))))
    lci[,state] <- 1/(1 + exp(-(log(probs2[,state]/(1-probs2[,state])) -
                                  qnorm(1-(1-alpha)/2) * (1/(probs2[,state]-probs2[,state]^2)) * se)))
    uci[,state] <- 1/(1 + exp(-(log(probs2[,state]/(1-probs2[,state])) +
                                  qnorm(1-(1-alpha)/2) * (1/(probs2[,state]-probs2[,state]^2)) * se)))
  }
  
  stat.pro.for <- tempCovs %>%
    mutate( 
      statprob1 = probs2$state1,statprob2 = probs2$state2,
      statprob_lci1 =  as.vector(lci[,1]),statprob_uci1 =  as.vector(uci[,1]),
      statprob_lci2 =  as.vector(lci[,2]),statprob_uci2 =  as.vector(uci[,2]),
    )
  
  # Stationary probability memory
  cov = 2
  othercova = 1
  gridLength <- 101
  hGridLength <- gridLength
  inf <- min(rawCovs[,cov],na.rm=TRUE)
  sup <- max(rawCovs[,cov],na.rm=TRUE)
  tempCovs <- NA
  tempCovs <- data.frame(
    cov1 =  rep(seq(inf,sup,length=gridLength),each=hGridLength/gridLength),
    cov2  = rep(covs[1,othercova], gridLength))
  colnames(tempCovs) <-
    c(names(m$rawCovs[cov]),names(m$rawCovs[othercova]))
  tmpcovs<-covs[names(rawCovs)]
  splineList=NULL
  tmpSplineInputs<-getSplineFormula(newformula,data,tempCovs)
  desMat <- model.matrix(tmpSplineInputs$formula, data=tmpSplineInputs$covs)
  probs <- stationary(m, covs=desMat)  
  probs2 <- as.data.frame(probs)
  
  lci <- matrix(NA,gridLength,length(ref))
  uci <- matrix(NA,gridLength,length(ref))
  
  state = c(1,2)
  for(state in 1:length(state)){
    dN <- t(apply(desMat, 1, function(x)
      numDeriv::grad(get_stat,m$mod$estimate[gamInd[unique(c(m$conditions$betaCons))]],covs=matrix(x,1),nbStates=nbStates,i=state,betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=m$conditions$workBounds$beta,mixture=mixture,ref=ref)))
    tmpSig <- Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]
    se <- t(apply(dN, 1, function(x)
      suppressWarnings(sqrt(x%*%tmpSig%*%x))))
    lci[,state] <- 1/(1 + exp(-(log(probs2[,state]/(1-probs2[,state])) -
                                  qnorm(1-(1-alpha)/2) * (1/(probs2[,state]-probs2[,state]^2)) * se)))
    uci[,state] <- 1/(1 + exp(-(log(probs2[,state]/(1-probs2[,state])) +
                                  qnorm(1-(1-alpha)/2) * (1/(probs2[,state]-probs2[,state]^2)) * se)))
  }
  
  stat.pro.mem <- tempCovs %>%
    mutate( 
      statprob1 = probs2$state1,statprob2 = probs2$state2,
      statprob_lci1 =  as.vector(lci[,1]), statprob_uci1 =  as.vector(uci[,1]),
      statprob_lci2 =  as.vector(lci[,2]), statprob_uci2 =  as.vector(uci[,2]),
    )
  
  #Transition probabilities foraging index
  covIndex <- 1:ncol(m$rawCovs)
  cov = 1
  othercova = 2
  gridLength <- 101
  hGridLength <- gridLength*ifelse(inherits(m,"hierarchical"),nlevels(m$data$level),1)
  inf <- min(m$rawCovs[,cov],na.rm=TRUE) #Minimum value of the range
  sup <- max(m$rawCovs[,cov],na.rm=TRUE) #Maximum value of the range
  tempCovs <- NA
  tempCovs <- data.frame(
    cov1 =  rep(seq(inf,sup,length=gridLength),each=hGridLength/gridLength),
    cov2  = rep(covs[1,othercova], gridLength))
  colnames(tempCovs) <- c(names(m$rawCovs[cov]),names(m$rawCovs[othercova]))
  splineList=NULL
  tmpSplineInputs<-getSplineFormula(newformula,m$data,tempCovs)
  desMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
  #desMat <- model.matrix(newformula, tempCovs)
  #Create matrix with all the transition probabilities 
  beta <- list(beta=m$mle$beta)
  trMat <- trMatrix_rcpp(nbStates,beta$beta[(mixture-1)*(nbCovs+1)+1:(nbCovs+1),,drop=FALSE],desMat,m$conditions$betaRef)
  
  i <- 2 #Assumes state 2 represents Transit
  second_state <- c(2,1) #State to transition to
  lci <- matrix(NA,gridLength,length(ref))
  uci <- matrix(NA,gridLength,length(ref))
  for(h in 1:length(second_state)){
    j = second_state[h]
    Sigma <- m$mod$Sigma
    tmpSig <- Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]
    dN<-t(apply(desMat,1,function(x) tryCatch(numDeriv::grad(get_gamma,m$mod$estimate[gamInd[unique(c(m$conditions$betaCons))]],covs=matrix(x,1,dimnames=list(NULL,names(x))),nbStates=nbStates,i=ref[i],j=ref[j],betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=m$conditions$workBounds$beta,mixture=mixture),error=function(e) NA)))
    se<-t(apply(dN,1,function(x) tryCatch(suppressWarnings(sqrt(x%*%tmpSig%*%x)),error=function(e) NA)))
    lci[,h]<-1/(1+exp(-(log(trMat[ref[i],ref[j],]/(1-trMat[ref[i],ref[j],]))-quantSup*(1/(trMat[ref[i],ref[j],]-trMat[ref[i],ref[j],]^2))*se)))
    uci[,h]<-1/(1+exp(-(log(trMat[ref[i],ref[j],]/(1-trMat[ref[i],ref[j],]))+quantSup*(1/(trMat[ref[i],ref[j],]-trMat[ref[i],ref[j],]^2))*se)))
  }
  
  tran.pro.for <- tempCovs %>%
    mutate( 
      prob22 = trMat[ref[2],ref[2],], prob21 = trMat[ref[2],ref[1],],
      tranprob22_lci =  as.vector(lci[,1]), tranprob22_uci =  as.vector(uci[,1]),
      tranprob21_lci =  as.vector(lci[,2]), tranprob21_uci =  as.vector(uci[,2]),
    )
  
  #Transition probabilities memory
  covIndex <- 1:ncol(m$rawCovs)
  cov = 2
  othercova = 1
  gridLength <- 101
  hGridLength <- gridLength*ifelse(inherits(m,"hierarchical"),nlevels(m$data$level),1)
  inf <- min(m$rawCovs[,cov],na.rm=TRUE) #Minimum value of the range
  sup <- max(m$rawCovs[,cov],na.rm=TRUE) #Maximum value of the range
  tempCovs <- NA
  tempCovs <- data.frame(
    cov1 =  rep(seq(inf,sup,length=gridLength),each=hGridLength/gridLength),
    cov2  = rep(covs[1,othercova], gridLength))
  colnames(tempCovs) <- c(names(m$rawCovs[cov]),names(m$rawCovs[othercova]))
  splineList=NULL
  tmpSplineInputs<-getSplineFormula(newformula,m$data,tempCovs)
  desMat <- model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
  #desMat <- model.matrix(newformula, tempCovs)
  #Create matrix with all the transition probabilities 
  beta <- list(beta=m$mle$beta)
  trMat <- trMatrix_rcpp(nbStates,beta$beta[(mixture-1)*(nbCovs+1)+1:(nbCovs+1),,drop=FALSE],desMat,m$conditions$betaRef)
  
  i <- 2 #Assumes state 2 represents Transit
  second_state <- c(2,1) #State to transition to
  lci <- matrix(NA,gridLength,length(ref))
  uci <- matrix(NA,gridLength,length(ref))
  for(h in 1:length(second_state)){
    j = second_state[h]
    Sigma <- m$mod$Sigma
    tmpSig <- Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]
    dN<-t(apply(desMat,1,function(x) tryCatch(numDeriv::grad(get_gamma,m$mod$estimate[gamInd[unique(c(m$conditions$betaCons))]],covs=matrix(x,1,dimnames=list(NULL,names(x))),nbStates=nbStates,i=ref[i],j=ref[j],betaRef=m$conditions$betaRef,betaCons=m$conditions$betaCons,workBounds=m$conditions$workBounds$beta,mixture=mixture),error=function(e) NA)))
    se<-t(apply(dN,1,function(x) tryCatch(suppressWarnings(sqrt(x%*%tmpSig%*%x)),error=function(e) NA)))
    lci[,h]<-1/(1+exp(-(log(trMat[ref[i],ref[j],]/(1-trMat[ref[i],ref[j],]))-quantSup*(1/(trMat[ref[i],ref[j],]-trMat[ref[i],ref[j],]^2))*se)))
    uci[,h]<-1/(1+exp(-(log(trMat[ref[i],ref[j],]/(1-trMat[ref[i],ref[j],]))+quantSup*(1/(trMat[ref[i],ref[j],]-trMat[ref[i],ref[j],]^2))*se)))
  }
  
  tran.pro.mem <- tempCovs %>%
    mutate( 
      prob22 = trMat[ref[2],ref[2],], prob21 = trMat[ref[2],ref[1],],
      tranprob22_lci =  as.vector(lci[,1]), tranprob22_uci =  as.vector(uci[,1]),
      tranprob21_lci =  as.vector(lci[,2]), tranprob21_uci =  as.vector(uci[,2]),
    )
}
#Now you can create the plots
color1 = "#1EB4E6"
color2 = "#06368F"

stationary.foraging <- ggplot(stat.pro.for)+
  geom_line(aes(x=batch.foraging.index, y = statprob1, colour="ARS"), lwd=1,lty=5)+
  geom_ribbon(aes(ymin=statprob_lci1, ymax=statprob_uci1, x=batch.foraging.index), alpha=0.3, fill=color2)+
  geom_line(aes(x=batch.foraging.index, y = statprob2, col="Transit"), lwd=1)+
  geom_ribbon(aes(ymin= statprob_lci2, ymax= statprob_uci2, x=batch.foraging.index), alpha=0.3, fill=color1)+
  ylim(0,1)+xlim(min(stat.pro.for$batch.foraging.index),max(stat.pro.for$batch.foraging.index))+
  xlab("Mean batch prey encounters")+ ylab("Stationary probability")+
  scale_color_manual(values=c(
    "Transit" = color1,
    "ARS" = color2))+ 
  labs(color="State")+
  theme_classic()+theme(legend.position="top", text=element_text(size=10))

stationary.memory <- ggplot(stat.pro.mem)+
  geom_line(aes(x=memory, y = statprob1, colour="ARS"), lwd=1,lty=5)+
  geom_ribbon(aes(ymin=statprob_lci1, ymax=statprob_uci1, x=memory), alpha=0.3, fill=color2)+
  geom_line(aes(x=memory, y = statprob2, col="Transit"), lwd=1)+
  geom_ribbon(aes(ymin= statprob_lci2, ymax= statprob_uci2, x=memory), alpha=0.3, fill=color1)+
  ylim(0,1)+xlim(min(stat.pro.mem$memory),max(stat.pro.mem$memory))+
  xlab("Proportion of dive batches spent searching")+ ylab("Stationary probability")+
  scale_color_manual(values=c(
    "Transit" = color1,
    "ARS" = color2))+ 
  labs(color="State")+
  theme_classic()+theme(legend.position="top", text=element_text(size=10))

color3 = "#9CCC62"
color4 = "#006B15"

transition.foraging <- ggplot(tran.pro.for)+
  geom_line(aes(x=batch.foraging.index, y = prob22, colour="Remain in Transit"), lwd=1)+
  geom_ribbon(aes(ymin=tranprob22_lci, ymax=tranprob22_uci, x=batch.foraging.index), alpha=0.3, fill=color3)+
  geom_line(aes(x=batch.foraging.index, y = prob21, col="From Transit to ARS"), lwd=1, lty=5)+
  geom_ribbon(aes(ymin=tranprob21_lci, ymax=tranprob21_uci, x=batch.foraging.index), alpha=0.3, fill=color4)+
  ylim(0,1)+xlim(min(tran.pro.for$batch.foraging.index),max(tran.pro.for$batch.foraging.index))+
  xlab("Mean batch prey encounters")+ ylab("Transition probability")+
  scale_color_manual(values=c(
    "Remain in Transit" = color3,
    "From Transit to ARS" = color4))+
  labs(color="State transition")+
  theme_classic()+theme(legend.position="top", text=element_text(size=10))

transition.memory <- ggplot(tran.pro.mem)+
  geom_line(aes(x=memory, y = prob22, colour="Remain in Transit"), lwd=1)+
  geom_ribbon(aes(ymin=tranprob22_lci, ymax=tranprob22_uci, x=memory), alpha=0.3, fill=color3)+
  geom_line(aes(x=memory, y = prob21, col="From Transit to ARS"), lwd=1, lty=5)+
  geom_ribbon(aes(ymin=tranprob21_lci, ymax=tranprob21_uci, x=memory), alpha=0.3, fill=color4)+
  ylim(0,1)+xlim(min(tran.pro.mem$memory),max(tran.pro.mem$memory))+
  xlab("Proportion of dive batches spent searching")+ ylab("Transition probability")+
  scale_color_manual(values=c(
    "Remain in Transit" = color3,
    "From Transit to ARS" = color4))+
  labs(color="State transition")+
  theme_classic()+theme(legend.position="top", text=element_text(size=10))


stationary <- ggarrange(stationary.foraging, stationary.memory, ncol=2, nrow=1, common.legend = TRUE)
#transitionbw <- ggarrange(stationary.foragingbw, membw,  ncol=2, nrow=1, common.legend = TRUE)
transition <- ggarrange(transition.foraging, transition.memory, ncol=2, nrow=1, common.legend = TRUE)
#stationarybw <- ggarrange(catchstatbw, memstatbw,  ncol=2, nrow=1, common.legend = TRUE)

fig3 <- ggarrange(stationary, transition, nrow=2, labels = c("A", "B"))
ggsave(plot = fig3, filename=here::here("Figures","Model 2 - stat and tran probabilities.tiff"),
       device="tiff", width = 220, height=140, units="mm", dpi=300)  

## Table S1 - seals data available --------------------------------------------------------------------------------------
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


## Table S2- Table summary data -------------------------------------------------------
obs <- read.csv(here::here("Output", "Repeatability - overlap between April and May.csv"), header=TRUE)
obs <- obs %>%  mutate(curve = rep("obs", length(obs$seal_ID)))
obs <- obs[,-5]
obs <- obs[-which(obs$seal_ID==59 | obs$seal_ID==280),]
obs <- obs[order(obs$seal_ID),]
colnames(obs) <- c("Seal ID", "# Trips in April", "# Trips in May", "BA overlap")


ft <- flextable(obs) %>% autofit(.) %>% align_nottext_col(., align="center")

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = here::here("Output", "Table - UD analysis overlap.docx"))



## Figure S1 - PR plots and parameters density distributions ---------------------------------------------------------------------------
#This code requires the momentuHMM extra functions script and for model 1 to be present in the Environment
source(here::here("R","Figures - MomentuHMM specific functions.R"))

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
  scale_colour_manual(name="HMM State", values=c("Transit"="#1EB4E6","ARS"="#06368F"))+
  ylab("Frequency density\n")+xlab("\nStep length (m)")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="transparent"),
        panel.grid.minor = element_line(colour="transparent"))

angle <- ggplot()+
  geom_histogram(data=df.HMM, aes(x=angle, y=..density..), binwidth=0.4 , fill="lightgrey",  boundary=0)+
  geom_path(data=angle1, aes(x=grid, y=V2, colour="Transit"), lwd=1)+
  geom_path(data=angle2, aes(x=grid, y=V2, colour="ARS"), lwd=1)+
  scale_colour_manual(name="HMM State", values=c("Transit"="#1EB4E6","ARS"="#06368F"))+
  ylab("Frequency density\n")+xlab("\nTurning angle (radians)")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="transparent"),
        panel.grid.minor = element_line(colour="transparent"))

m1density <- ggarrange(step, angle, common.legend=TRUE)  

ggsave(plot= m1density, filename=here::here("Figures","Model 1 - parameters density functions.tiff"), 
       device="tiff", width = 10, height=5, units="in", dpi=300)    


## Figure S2 - raw data during dive -----------------------------------------------------------------
#Raw accelerometer dive data
source(here::here("R","06. Prey encounter - Accelerometer processing functions.R"))
acc.data <- read.csv(here::here("Figures","Data for figures","pv14464_125_foraging_trip_data.csv"), header=TRUE) %>%
  mutate(datetime = ddate,
         time = as.POSIXct(seconds, format="%Y-%m-%d %H:%M:%OS", tz="UTC"),
         pitch_rad = asin(Xs),
         pitch_deg = rad2deg(pitch_rad),
         roll_rad = asin(Ys),
         roll_deg = rad2deg(roll_rad))

acc.depth <- read.table(here::here("Figures","Data for figures","14464_accel_depth.txt"), sep="\t", header=TRUE) %>%
  mutate(time = as.POSIXct(DATE..UTC., format="%Y/%m/%d %H:%M:%OS", tz="UTC"))

dive.data <- read.csv(here::here("Figures","Data for figures","pv14464_125_foraging_trip_dive_data.csv"), header=TRUE) %>%
  mutate(time = as.POSIXct(DS_DATE, format="%Y-%m-%d %H:%M:%OS", tz="UTC"))

#Specify the times of the dive you want to plot
dive.start <- as.POSIXct("2017-06-01 03:36:08",  format="%Y-%m-%d %H:%M:%S", tz="UTC")
dive.end <- as.POSIXct("2017-06-01 03:39:28",  format="%Y-%m-%d %H:%M:%S", tz="UTC")

{#Prepare data for plot
  acc.dive <- acc.data[which(acc.data$time >= dive.start & acc.data$time <= dive.end),]
  
  #Calculate standard deviation of dynamic acceleration
  Xsd <- as.matrix(acc.dive$Xd)
  Ysd <- as.matrix(acc.dive$Yd)
  Zsd <- as.matrix(acc.dive$Zd)
  
  acc.dive <- acc.dive %>% mutate(
    Xsd = roll_sd(Xsd, 19),
    Ysd = roll_sd(Ysd, 19),
    Zsd = roll_sd(Zsd, 19),
    highX = ifelse(Xsd>=0.2030459, 1, 0),
    highY = ifelse(Ysd>=0.170182, 1, 0),
    highZ = ifelse(Zsd>=0.1373276, 1, 0),
    highstate = highX+highY+highZ
  )
  
  #Calculate prey catch attempts events
  attempts <- acc.dive[which(acc.dive$highstate==3),]
  attempts$posix_sec <- as.POSIXct(attempts$seconds, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
  for(x in 1:length(attempts$ddate)){
    attempts$diff.time[x] <- as.numeric(difftime(attempts$posix_sec[x+1],attempts$posix_sec[x], unit="sec"))
  }
  attempts2 <- attempts[which(attempts$diff.time>=1),]
  df_Prca2 <- data.frame("attempts" = attempts2)
  
  row_end <- which(attempts$diff.time>=1)
  row_start <- row_end+1
  row_start <- c(1, row_start)
  row_end <- c(row_end, length(attempts$ddate))
  time_PrCA <- data.frame(start = attempts$posix_sec[row_start], end = attempts$posix_sec[row_end])
  time_PrCA$y <- rep(0.8, length(time_PrCA$start))
  
  #Calculate benthic attempts
  ben <- BENTHIC_prca("acc.dive", acc.dive$pitch_deg, 20,50,30) 
  ben$peak_time <- as.POSIXct(ben$peak_time, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  ben$valley_time <- as.POSIXct(ben$valley_time, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  bt_preydf <- ben[which(ben$drop_duration>0),]
  bt_preydf$y <- rep(0, length(bt_preydf$peak_time))
  
  
  peaks <- as.data.frame(find_peaks(acc.dive$pitch_deg, m=50))
  valley <- as.data.frame(find_peaks(-acc.dive$pitch_deg, m= 30))
  colnames(peaks) <- c("peak")
  colnames(valley) <- c("peak")
  peaks$lab <- rep("peak", length(peaks))
  valley$lab <- rep("valley", length(valley))
  pv <- rbind(peaks, valley)
  pv <- pv[order(pv$peak),]
  LEN <- length(acc.dive$X)
  acc.dive$row.num <- seq(1, LEN, 1)
  dftmp <- merge(acc.dive,pv, by.x="row.num", by.y="peak", all.x=TRUE)
  dftmp <- dftmp[complete.cases(dftmp$lab),]
  dftmp <- spread(dftmp, lab, pitch_deg)
  if(is.na(dftmp$peak[1])){dftmp$peak[1] <- 0}else{}
  dftmp <- dftmp %>% mutate(
    start_peak = ifelse(is.na(peak), NA, paste(seconds)),
    row_peak = ifelse(is.na(peak), NA, paste(row.num)),
    peak = fill.na(peak),
    start_peak = fill.na(start_peak),
    row_peak = fill.na(row_peak)
  )
  dftmp <- dftmp[complete.cases(dftmp$valley),]
  dftmp <- dftmp %>% mutate(
    diff = (peak - valley),
    time_diff = as.numeric(difftime(seconds, start_peak, units="secs")),
    attempts = ifelse(diff>=20 & time_diff<=5 & time_diff>0, 2, 0),
    dattempts = ifelse(diff>=20 & time_diff>5, 1, attempts),
    start_peak = as.POSIXct(start_peak, format="%Y-%m-%d %H:%M:%OS", tz="UTC"),
    time = as.POSIXct(seconds, format="%Y-%m-%d %H:%M:%OS", tz="UTC"),
    adj.diff = -80 + diff,
    mid.drop = start_peak+time_diff/2
  )
  dftmp <- dftmp[-which(dftmp$row.num==1533),]
  dftmp <- dftmp[-which(dftmp$row.num==1805),]
  ben_attempts <- dftmp[which(dftmp$attempts==2),]
  ben_attempts <- ben_attempts[!duplicated(attempts$row_peak),]
  
  #Arhcived depth profile
  adepth.dive <- acc.depth[which(acc.depth$time >= dive.start & acc.depth$time <= dive.end),]
  adepth.dive$DEPTH..m. <- - adepth.dive$DEPTH..m.
  
  #Transmitted depth profile
  bdepth.dive <- dive.data[which(dive.data$time == dive.start ),]
  bdepth <- stack(bdepth.dive[,38:60])
  bdepth.time <- stack(bdepth.dive[,15:37])
  dive.duration <- as.numeric(difftime(dive.end, dive.start, unit="secs"))
  bdepth.time$change <- bdepth.time$values*dive.duration/100
  bdepth$change <- bdepth.time$change
  bdepth$time <- dive.start + bdepth$change
  bdepth$values <- - bdepth$values
}

depth_plot <- ggplot()+
  geom_path(data=adepth.dive, aes(x=time, y=DEPTH..m.,color="a"),lwd=1)+
  geom_point(data=bdepth, aes(x=time, y=values, color="b"))+
  geom_path(data=bdepth, aes(x=time, y=values, color="b"))+
  scale_color_manual(
    values = c("black", "red"),
    labels = c("Archived", "Transmitted"),
    name = c(" ")
  )+
  xlab(" ")+ylab("Depth (m)")+
  theme_classic()+
  theme(legend.position = c(0.7, 0.8))

dynamic_plot <- ggplot()+
  geom_path(data=acc.dive, aes(x=time, y=Xsd, color="a"))+
  geom_path(data=acc.dive, aes(x=time, y=Ysd, color="b"))+
  geom_path(data=acc.dive, aes(x=time, y=Zsd, color="c"))+
  geom_point(data=time_PrCA, aes(x=start, y=y, color="d"))+
  scale_color_manual(values = c("red", "blue", "dark green", "black"),
                     labels = c("St. dev x-axis","St. dev y-axis","St. dev z-axis", "PrCA behaviour"),
                     name = c(" "))+
  xlab(" ")+ylab("Standard deviation\n dynamic acc. (g)")+
  theme_classic()+
  theme(legend.position = "top",
        legend.key = element_rect(colour = NA, fill = NA))

y_plot <- ggplot(acc.dive, aes(x=time, y=Yd))+
  geom_path()+
  xlab("Time")+ylab("Dynamic\n acceleration\n Y axis (g)")+
  theme_classic()

pitch_plot <- ggplot()+
  geom_path(data=acc.dive, aes(x=time, y=pitch_deg))+
  geom_point(data=bt_preydf[2:4,], aes(x=peak_time, y=y, col="a"))+
  geom_segment(data=dftmp[6:(length(dftmp$row.num)-3),], aes(x=mid.drop, y=-80, xend=mid.drop, yend=adj.diff), col="orange")+
  geom_hline(yintercept=-60, linetype="dashed", color="orange")+
  scale_color_manual(values = c("red"),
                     labels = c("Benthic attempts"),
                     name= c(" "))+
  xlab(" ")+ylab("Body pitch\n angle (degrees)")+
  theme_classic()+
  theme(legend.position="top")

ggarrange(depth_plot, dynamic_plot, pitch_plot, y_plot, ncol=1, heights = c(1, .8, .8, .6),labels = c("A","B","C","D"))
ggsave(here::here("Figures","Accelerometer - dive example.jpg"),
       device="jpg", width = 6, height = 12, units="in")





  


## Figure S3 - Foraging tactics --------------------------------------------------------------------
path <- here::here("Output","Processed accelerometer parameters")
files <- paste0(path,"/",list.files(paste0(path)))
tables <- lapply(files, read.csv)
dat <- do.call(rbind, tables) 

dat <- dat [,c(1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,25,26,27,33,34,37,47)]
dat$dive_type <- NA
n <- 1
for(n in 1:length(dat$ID)){
  if(is.na(dat$B_PrCA_arch[n])){
    dat$dive_type[n] <- NA
  } else {
    if(dat$B_PrCA_arch[n]>0 & dat$B_pitch.diff20[n]==0){dat$dive_type[n] <- "1"} #Only prca
    if(dat$B_PrCA_arch[n]==0 & dat$B_pitch.diff20[n]>0){dat$dive_type[n] <- "2"} #only benthic
    if(dat$B_PrCA_arch[n]>0 & dat$B_pitch.diff20[n]>0){dat$dive_type[n] <- "3"} #both
    if(dat$B_PrCA_arch[n]==0 & dat$B_pitch.diff20[n]==0){dat$dive_type[n] <- "4"}
  }
} 

graph.plot <- dat %>% group_by(ID, dive_type) %>% summarise (count=n())
graph.plot <- graph.plot[which(graph.plot$dive_type>0),]
tot.dives <- graph.plot %>% group_by(ID) %>% summarise(n.dives = sum(count))
graph.plot$tot.dives <- rep(tot.dives$n.dives, each=4)
graph.plot$p.dives <- graph.plot$count/graph.plot$tot.dives

catching.strategy <- ggplot(graph.plot, aes(x=as.factor(ID), y=p.dives, fill=dive_type))+
  geom_bar(position="dodge", stat = "identity")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                    name= "Dive type",
                    labels=c("Peaks in dynamic acceleration", "Changes in body pitch angle", "Both catch strategies", "No attempts"))+
  ylab("Proportion of dives")+
  xlab("Seal ID")+
  theme_classic()

ggsave(here::here("Figures","Accelerometer - catching strategies.tiff"),
       device="tiff", width = 10, height = 6, units="in")


## Figure S4 - PR plots and parameters density distributions ---------------------------------------------------------------------------
#This code requires the momentuHMM extra functions script and for model 1 to be present in the Environment
source(here::here("R","Figures - MomentuHMM specific functions.R"))

# Evaluation plots model 1
tiff(here::here("Figures","Model 2 - evaluation plots.tiff"), width = 10, height = 6, units="in", res=300)
plotPR(m2)
dev.off()

#Need to debug to obtain the distribution densities
debug(momentuHMM:::plot.momentuHMM)
plot(m2)
#Debug untill line 620, check that genDensities is a list of 2, or the number of your states
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

df.HMM <- read.csv(here::here("Output", "Model 2 - HMM dive batches classified.csv"), header=TRUE)

step <- ggplot()+
  geom_histogram(data=df.HMM, aes(x=step, y=..density..), binwidth=500 , fill="light grey",  boundary=0)+
  geom_path(data=step1, aes(x=grid , y=V2, colour="ARS"), lwd=1)+
  geom_path(data=step2, aes(x=grid, y=V2, colour="Transit"), lwd=1)+
  scale_colour_manual(name="HMM State", values=c("Transit"="#1EB4E6","ARS"="#06368F"))+
  ylab("Frequency density\n")+xlab("\nStep length (m)")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="transparent"),
        panel.grid.minor = element_line(colour="transparent"))

angle <- ggplot()+
  geom_histogram(data=df.HMM, aes(x=angle, y=..density..), binwidth=0.4 , fill="lightgrey",  boundary=0)+
  geom_path(data=angle1, aes(x=grid, y=V2, colour="ARS"), lwd=1)+
  geom_path(data=angle2, aes(x=grid, y=V2, colour="Transit"), lwd=1)+
  scale_colour_manual(name="HMM State", values=c("Transit"="#1EB4E6","ARS"="#06368F"))+
  ylab("Frequency density\n")+xlab("\nTurning angle (radians)")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour="transparent"),
        panel.grid.minor = element_line(colour="transparent"))

m2density <- ggarrange(step, angle, common.legend=TRUE)  

ggsave(plot= m2density, filename=here::here("Figures","Model 2 - parameters density functions.tiff"), 
       device="tiff", width = 10, height=5, units="in", dpi=300)  

