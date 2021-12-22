### Processing accelerometer data
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk) and Matt Carter (midc@st-andrews.ac.uk)
# Purpose: This code uses the functions specified in script 6 and other calculations to calculate various 
#          dive parameters and summarise the acceleormeter data in meaningful ways. 
#          Dive metrics: - descent/bottom/ascent phase start time
#                        - TAD (Time-at-depth index)
#                        - descent/ascent speed
#          Accelerometer: - Prey catch attempts using archived method from Cox et al. 2018 and benthic attempts from Brasseur et al. 2012 only for the bottom phase
#                         - overlap attempts between the benthic attempts and prey catch attempts


#References:
#Brasseur, S., et al. (2012). "Habitat Preferences of Harbour Seals in the Dutch Coastal Area: Analysis and Estimate of Effects of Offshore Wind Farms (Report No. OWEZ R 252 T1 20120130 C043-10), IMARES - Wageningen UR, Noordzeewind: 58".
#Cox, S. L., et al. (2018). "Processing of acceleration and dive data on board satellite relay tags to investigate diving and foraging behaviour in free ranging marine predators." Methods in Ecology and Evolution 9(1): 64-77".

# Output: A final dataset with rows representing individual dives and the parameters for that dive.
# Created on: 12/01/2021
# Updated on: 21/12/2021

## Load packages --------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, magrittr, ggplot2, roll, pracma, sf, rgdal, mefa, ggsn)

#Load functions
source(here::here("R","08. Prey encounter - Accelerometer processing functions.R"))

## Data preparation --------------------------------------------------------------------------------------
#Dive data
all_seal <- read.table(here::here("Dryad", "pv64-2017_dive.txt"), sep="\t",header=TRUE) %>%
  mutate(DS_DATE = as.POSIXct(DS_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         DE_DATE = as.POSIXct(DE_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"))

all_trips <- read.table(here::here("Dryad", "pv64-2017_trip_summaries.txt"),sep="\t" ,header=TRUE) %>%
  mutate(Trip_Start = as.POSIXct(Trip_Start , format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         Trip_End = as.POSIXct(Trip_End , format="%Y-%m-%d %H:%M:%S", tz="UTC"))

#Run one seal at the time
#These values need to be manually modified for each seal
PTTacc = 14478
IDN = 158
seal <- all_seal[which(all_seal$ID == IDN),]
trip_seal <- all_trips[which(all_trips$ID == IDN),]

#Make list of the trips for which you accelerometer data
files <- list.files(paste0(here::here("Dryad","Outputs","Raw accelerometer data per trip"),"/",PTTacc,"/"))
split <- do.call(rbind, strsplit(files, "_"))
numbers <- as.numeric(split[,2])
trip_seal <- trip_seal[which(trip_seal$Trip_No %in% numbers),]
trips <- unique(trip_seal$Trip_No)

#For the PrCA analysis need to load the previously calculates thresholds
thresholds <- read.csv(here::here("Dryad","Outputs","Prey encounters - cluster analysis axis thresholds.txt"), 
                       header=TRUE, sep="\t")

Sys.setenv(TZ='GMT')
sealdf <- data.frame(ID= NA, REF = NA, PTT = NA, DS_DATE = as.POSIXct(NA), JUL =NA, DE_DATE = as.POSIXct(NA),
                     DIVE_DUR =NA, SURF_DUR = NA, MAX_DEP = NA,START_LAT = NA, START_LON =NA,
                     END_LAT = NA, END_LON = NA, PERCENT_AREA =NA, Trip_No =NA, Posn_in_Trip = NA, 
                     Trip_Code =NA, date = NA, sunrise =NA, sunset= NA, daynight = NA,
                     descent_start = as.POSIXct(NA), bottom_start = as.POSIXct(NA),
                     ascent_start = as.POSIXct(NA), PrCA = NA, pitch.diff20 = NA, overlap = NA)

sealdfempty <- sealdf
tripsdf <- sealdf
skipped <- NA

for(l in 1:length(trip_seal$ID)){
    if(trip_seal$Trip_Duration[l]>1){
      trip <- seal[which(seal$Trip_No==trip_seal$Trip_No[l]),]
    if(length(trip$ID)>0){
    ## Define bottom phase of the dive -------------------------------------------------------------------------
    #Extract a dive threshold such as the 80% of the maximum depth of the dive
    trip$DThresh<-trip$MAX_DEP/100*80
    trip$DThresh <- round(trip$DThresh, digits=1)
    #Calculate the percentage of the dive which was at the bottom
    trip$BTcols<-rowSums(trip[,which(colnames(trip)=="D1"):which(colnames(trip)=="D23")] > trip$DThresh) 
    trip$BTprop<-trip$BTcols/23 #ignore T1:T2 and T20:T21 because they are not equally spaced
    #Calculate bottom time
    trip$BT<- trip$BTprop*trip$DIVE_DUR
    
    #Find the infl_prop (inflection proportion) at which you can split the dive between descent and bottom, and bottom and ascent. 
    trip$infl_prop <- NA
    trip$infl_prop2 <- NA
    for(y in 1:length(trip$ID)){
    tmp <- trip[y,]
    tmpT <- gather(tmp[which(colnames(tmp)=="T1"):which(colnames(tmp)=="T23")], "Ts", "percentage")
    tmpD <- gather(tmp[which(colnames(tmp)=="D1"):which(colnames(tmp)=="D23")], "Ds", "depth")
    tmp <- cbind(tmpT, tmpD)
    thrsh <- trip$DThresh[y]
    
    if(tmp$depth[1]>=thrsh){
      trip$infl_prop[y] <- 0
    } else{
      for(r in 1:23){
      if(tmp$depth[r]>=thrsh){
        trip$infl_prop[y] <- tmp$percentage[r]; break}  #infl_prop beginning of bottom phase
      }
    }
    if(tmp$depth[23]>=thrsh){
      trip$infl_prop2[y] <- 100
    } else{
      for(r in 23:1){
        if(tmp$depth[r]>=thrsh){
          trip$infl_prop2[y] <- tmp$percentage[r]; break} #infl_prop2 end of bottom phase
      }
     }
    }
    
    #define beginng of diving phases
    trip$bphase_sec_start <- ((trip$infl_prop/100)*trip$DIVE_DUR)
    trip$bphase_start <- trip$DS_DATE + trip$bphase_sec_start
    trip$bphase_sec_end <- ((trip$infl_prop2/100)*trip$DIVE_DUR)
    trip$bphase_end<- trip$DS_DATE + trip$bphase_sec_end
    
    
    #Get rid of sll the time depth columns
    seal_tmp <- trip %>% select(ID, REF, PTT, DS_DATE, JUL, DE_DATE, DIVE_DUR, SURF_DUR, MAX_DEP, 
                                START_LAT, START_LON, END_LAT, END_LON, PERCENT_AREA, Trip_No,
                                Posn_in_Trip, Trip_Code, date, sunrise, sunset, daynight, bphase_start,
                                bphase_end)
    
    seal_tmp$descent_start <- as.POSIXct(seal_tmp$DS_DATE, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
    colnames(seal_tmp)[which(colnames(seal_tmp)=="bphase_start")] <- c("bottom_start")
    colnames(seal_tmp)[which(colnames(seal_tmp)=="bphase_end")] <- c("ascent_start")
    
    #re-order columns to match sealdf
    seal_tmp <- seal_tmp[,c(1:21,24,22,23)]
    
    sealdf[1:length(seal_tmp$ID),1:length(seal_tmp)] <- seal_tmp

    ##Dive-by-dive accelerometer  ------------------------------------------------------------------------
    acc <- read.table(paste0(here::here("Dryad","Outputs","Raw accelerometer data per trip"),"/",PTTacc,"/pv",PTTacc,"_",trip_seal$Trip_No[l],"_trip_data.txt"), 
                      sep="\t", header=TRUE)
    acc$posix <- as.POSIXct(acc$seconds, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    acc$posix_sec <- as.POSIXct(acc$seconds, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    
    #Calculate pitch angle
    b<-(acc$Ys^2)+(acc$Zs^2)
    a<-acc$Xs
    acc$pitch <-atan((a/sqrt(b)))
    acc$pitch <- acc$pitch*(180/pi)
    acc$pitch_deg <- acc$pitch+90
    
    #Standard deviation of dynamic acceleration
    Xsd <- as.matrix(acc$Xd)
    acc$Xsd  <- roll_sd(Xsd, 19)
    Ysd <- as.matrix(acc$Yd)
    acc$Ysd  <- roll_sd(Ysd, 19)
    Zsd <- as.matrix(acc$Zd)
    acc$Zsd  <- roll_sd(Zsd, 19)
  
    thresholdX <- thresholds$ThresholdX[which(thresholds$PTT==PTTacc)]
    thresholdY <- thresholds$ThresholdY[which(thresholds$PTT==PTTacc)]
    thresholdZ <- thresholds$ThresholdZ[which(thresholds$PTT==PTTacc)]

    for(x in 1:length(sealdf$ID)){
      #subset accelerometer data for individual dive
      dive_acc <- subset(acc, acc$posix>= sealdf$DS_DATE[x] &
                           acc$posix<= sealdf$DE_DATE[x])
      if(length(dive_acc$ddate)>1){
        diff <-  as.numeric(difftime(dive_acc$ddate[length(dive_acc$ddate)], dive_acc$ddate[1], unit = "sec"))
        #Check the length of the dive
        if(diff>0) {
          ## Bottom phase - accelerometer ---------------------------------------------------------------------
          bt <- subset(dive_acc, dive_acc$posix>= sealdf$bottom_start[x] &
                              dive_acc$posix < sealdf$ascent_start[x])
          if(length(bt$ddate)>0){
            #Archived PrCA
            prca2 <- PRCA_arch("bt", ColXd = bt$Xd, ColYd = bt$Yd, ColZd = bt$Zd, 17, thresholdX, thresholdY,
                               thresholdZ)
            df_prcaA <- as.data.frame(prca2[1])
            sealdf$PrCA[x] <- df_prcaA$attempts
            time_PrCAA <- as.data.frame(prca2[2])
            
            #Benthic attempts
            ben <- BENTHIC_prca("bt", bt$pitch_deg, 20,50,30)
            if(length(ben)>1){
              sealdf$pitch.diff20[x] <- length(ben$peak_time)
              ben$peak_time <- as.POSIXct(ben$peak_time, format="%Y-%m-%d %H:%M:%S", tz="GMT")
              ben$valley_time <- as.POSIXct(ben$valley_time, format="%Y-%m-%d %H:%M:%S", tz="GMT")
            } else{
              sealdf$pitch.diff20[x] <- 0
            }
            
            if(df_prcaA$attempts > 0 & length(ben)>1){
              time_PrCAA$overlap <- NA
              for(e in 1:length(time_PrCAA$start)){
                for(t in 1:length(ben$peak_time)){
                  if(time_PrCAA$start[e]<=ben$peak_time[t] & time_PrCAA$end[e]>= ben$valley_time[t]){
                    time_PrCAA$overlap[e] <- 1
                  }
                }
              }
              sealdf$overlap[x] <- sum(time_PrCAA$overlap, na.rm=TRUE)
            } else{
              sealdf$overlap[x] <- 0
            }
          
          } else{
            sealdf[x,which(colnames(sealdf)=="PrCA"):which(colnames(sealdf)=="overlap")] <- NA
          }

        }
      }
    }
    tripsdf <- rbind(tripsdf, sealdf)
    sealdf <- sealdfempty
    print(paste0(l," Trip done out of ", length(trip_seal$ID)))
    } else{
      sealdf <- sealdfempty
      skipped <- c(skipped, trip_seal$Trip_No[l])
      print(paste0(l," Trip done out of ", length(trip_seal$ID)))
      }
    } else{
      sealdf <- sealdfempty
      skipped <- c(skipped, trip_seal$Trip_No[l])
      print(paste0(l," Trip done out of ", length(trip_seal$ID)))
    }
}

tripsdf <- tripsdf[-1,]

#Check that you have analysed all the trips
print(paste0("Total trips available for analysis ",length(which(trip_seal$Trip_Duration>1))))
print(paste0("Total trips analysed ",length(unique(tripsdf$Trip_No))))
skipped <- skipped[-1]
print(paste0("Total trips skipped because of missing data ",length(skipped)))
trip_processed <- length(unique(tripsdf$Trip_No))+length(skipped)
print(paste0("Total trips processed ",trip_processed))

write.table(tripsdf, paste0(here::here("Dryad","Outputs", "Processed accelerometer parameters"),"/",PTTacc,"_processed_accelerometer_parameters.txt"),
            sep="\t" ,row.names=FALSE)
