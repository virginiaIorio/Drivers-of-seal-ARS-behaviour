### Processing accelerometer data
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose: This code uses the functions specified in script 6 and other calculations to calculate various 
#          dive parameters and summarise the acceleormeter data in meaningful ways. 
#          Dive metrics: - descent/bottom/ascent phase start time
#                        - TAD (Time-at-depth index)
#                        - descent/ascent speed
#          Accelerometer: - stroke rate and swimming effort for all dive phases
#                         - Mean body pitch angle for all dive phases
#                         - Prey catch attempts using "transmitted" and "archived" method from Cox et al. 2018 and benthic attempts from Brasseur et al. 2012 only for the bottom phase
#                         - overlap attempts between the benthic attempts and either "transmitted" or "archived" methods
#                         - Mean roll angle during the bottom phase of te dive
#                         - Standard deviation of the dynamic acceleration during the bottom phase of the dive

#References:
#Brasseur, S., et al. (2012). Habitat Preferences of Harbour Seals in the Dutch Coastal Area: Analysis and Estimate of Effects of Offshore Wind Farms (Report No. OWEZ R 252 T1 20120130 C043-10), IMARES - Wageningen UR, Noordzeewind: 58.
#Cox, S. L., et al. (2018). "Processing of acceleration and dive data on‐board satellite relay tags to investigate diving and foraging behaviour in free‐ranging marine predators." Methods in Ecology and Evolution 9(1): 64-77.

# Output: A final dataset with rows representing individual dives and the parameters for that dive.
# Created on: 12/01/2021
# Updated on: 13/05/2021

## Load packages --------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, magrittr, ggplot2, roll, pracma, sf, rgdal, mefa, ggsn)

#Load functions
source(here::here("R","06. Prey encounter - Accelerometer processing functions.R"))

## Data preparation --------------------------------------------------------------------------------------
#Dive data
all_seal <- read.csv(here::here("Datasets", "pv64-2017_dive.csv"), header=TRUE) %>%
  mutate(DS_DATE = as.POSIXct(DS_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         DE_DATE = as.POSIXct(DE_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC"))

all_trips <- read.csv(here::here("Datasets", "pv64-2017_trip_summaries.csv"), header=TRUE) %>%
  mutate(Trip_Start = as.POSIXct(Trip_Start , format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         Trip_End = as.POSIXct(Trip_End , format="%Y-%m-%d %H:%M:%S", tz="UTC"))

#Run one seal at the time
#These values need to be manually modified for each seal
PTTacc = 14438
IDN = 90
seal <- all_seal[which(all_seal$ID == IDN),]
trip_seal <- all_trips[which(all_trips$ID == IDN),]

#Make list of the trips for which you accelerometer data
files <- list.files(paste0(here::here("Output","Raw accelerometer data per trip"),"/",PTTacc,"/"))
split <- do.call(rbind, strsplit(files, "_"))
numbers <- as.numeric(split[,2])
trip_seal <- trip_seal[which(trip_seal$Trip_No %in% numbers),]
trips <- unique(trip_seal$Trip_No)

#For the PrCA analysis need to load the previously calculates thresholds
thresholds <- read.csv(here::here("Output","Prey encounters - cluster analysis axis thresholds.csv"), header=TRUE)

Sys.setenv(TZ='GMT')
sealdf <- data.frame(ID= NA, REF = NA, PTT = NA, DS_DATE = as.POSIXct(NA), JUL =NA, DE_DATE = as.POSIXct(NA),
                     DIVE_DUR =NA, SURF_DUR = NA, MAX_DEP = NA,START_LAT = NA, START_LON =NA,
                     END_LAT = NA, END_LON = NA, PERCENT_AREA =NA, Trip_No =NA, Posn_in_Trip = NA, 
                     Trip_Code =NA, date = NA, sunrise =NA, sunset= NA, daynight = NA,
                     descent_start = as.POSIXct(NA), bottom_start = as.POSIXct(NA),
                     ascent_start = as.POSIXct(NA),  TAD=NA,  descent.speed = NA, 
                     ascent.speed = NA, D_stroke_rate = NA, D_swim_eff.g = NA,
                     D_swim_eff.m_s = NA, D_mean_pitch = NA, B_PrCA_trans = NA, B_PrCA_arch = NA,
                     B_pitch.diff20 = NA, B_overlap_T = NA, B_overlap_A = NA,
                     B_stroke_rate = NA, B_swim_eff.g = NA, B_swim_eff.m_s = NA, 
                     B_mean_roll = NA, B_dyn_sd = NA, B_mean_pitch = NA, 
                     A_stroke_rate = NA, A_swim_eff.g = NA, A_swim_eff.m_s = NA, A_mean_pitch = NA)

sealdfempty <- sealdf
tripsdf <- sealdf
skipped <- NA

for(l in 1:length(trip_seal$ID)){
    if(trip_seal$Trip_Duration[l]>1){
      trip1 <- seal[which(seal$Trip_No==trip_seal$Trip_No[l]),]
    if(length(trip$ID)>0){
    ## Define bottom phase of the dive -------------------------------------------------------------------------
    #Extract a dive threshold such as the 80% of the maximum depth of the dive
    trip$DThresh<-trip$MAX_DEP/100*80
    trip$DThresh <- round(trip$DThresh, digits=1)
    #Calculate the percentage of the dive which was at the bottom
    trip$BTcols<-rowSums(trip[,38:60] > trip$DThresh) 
    trip$BTprop<-trip$BTcols/23 #ignore T1:T2 and T20:T21 because they are not equally spaced
    #Calculate bottom time
    trip$BT<- trip$BTprop*trip$DIVE_DUR
    
    #Find the infl_prop (inflection proportion) at which you can split the dive between descent and bottom, and bottom and ascent. 
    trip$infl_prop <- NA
    trip$infl_prop2 <- NA
    for(y in 1:length(trip$ID)){
    tmp <- trip[y,]
    tmpT <- gather(tmp[15:37], "Ts", "percentage")
    tmpD <- gather(tmp[38:60], "Ds", "depth")
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
    
    ## TAD and descent-ascent speed --------------------------------------------------------------
    #For details look at the Matt Carter code
    DThresh<-1.5 #depth threshold defined by device
    trip$CorrDv<-trip$MAX_DEP-DThresh	
    x=seq(1,3,0.1)
    y=0
    dat1<-data.frame("x"=x,"y"=y)
    for(i in 1: dim(dat1)[1]){
      S <- x[i]
      TD <-(trip$PERCENT_AREA/100*trip$CorrD*trip$DIVE_DUR-trip$CorrD^2/S)/(trip$CorrD*trip$DIVE_DUR-2*trip$CorrD^2/S)
      temp<-subset(TD, TD<=1 & TD>=0.5) #subset to remove erroneous values
      temp2<-length(temp)
      if(i==1){TAD <- temp2}else
      {TAD <- rbind(TAD, temp2)} 
    }  
    dat2<-data.frame("x"=x,"y"=TAD)
    
    x<-dat2$x
    y<-dat2$y
    for(i in 1: dim(dat2)[1]){
      d1 <- (y[i+1]-y[i])/(x[i+1]-x[i])  #first derivative
      if(i==1){D <- d1}else
      {D <- rbind(D, d1)} 
    }  
    
    Der<-data.frame("x"=seq(1,3,0.1))
    Der$D<-D
    
    for(i in 1: dim(Der)[1]){
      d2 <- (D[i+1]-D[i])/((x[i+2]-x[i])/2)  
      if(i==1){D2 <- d2}else
      {D2 <- rbind(D2, d2)} 
    }  
    Der2<-data.frame("x"=seq(1.1,3.1,0.1))
    Der2$D2<-D2
    
    for (k in 1:dim(Der2)[1]){
      if(Der2$D2[1]>=0){
        S <-  min(which(Der2$D2>=0))} else{
      if (Der2$D2[k]<0 & Der2$D2[k+1]>=0){
        S <- Der2$x[k+1]; break} 
        }
      }

    DThresh<-1.5
    trip$CorrD<-trip$MAX_DEP-DThresh	
    trip$TAD<-(trip$PERCENT_AREA/100*trip$CorrD*trip$DIVE_DUR-trip$CorrD^2/S)/(trip$CorrD*trip$DIVE_DUR-2*trip$CorrD^2/S)
    
    trip$descent.speed <- abs(trip$D2-trip$D1)/((2.5-1)*(1/100)*(trip$DIVE_DUR-4))
    trip$ascent.speed <- (trip$D22-trip$D23)/((99-97.5)*(1/100)*(trip$DIVE_DUR-4))
    trip$ascent_dur <- trip$DIVE_DUR- trip$bphase_sec_end
    
    #Get rid of sll the time depth columns
    seal_tmp <- trip[,c(1:14,62:68,76,78,81:83)] #check the columns that you are extracting
    seal_tmp$descent_start <- as.POSIXct(seal_tmp$DS_DATE, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
    colnames(seal_tmp)[22] <- c("bottom_start")
    colnames(seal_tmp)[23] <- c("ascent_start")
    
    seal_tmp <- seal_tmp[,c(1:21,27,22:26)] #Re-organise columns order
    
    sealdf[1:length(seal_tmp$ID),1:27] <- seal_tmp

    ##Dive-by-dive accelerometer  ------------------------------------------------------------------------
    acc <- read.csv(paste0(here::here("Output","Raw accelerometer data per trip"),"/",PTTacc,"/pv",PTTacc,"_",trip_seal$Trip_No[l],"_trip_data.csv"), header=TRUE)
    acc$posix <- as.POSIXct(acc$seconds, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    acc$posix_sec <- as.POSIXct(acc$seconds, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    
    #Calculate pitch angle
    b<-(acc$Ys^2)+(acc$Zs^2)
    a<-acc$Xs
    acc$pitch <-atan((a/sqrt(b)))
    acc$pitch <- acc$pitch*(180/pi)
    acc$pitch_deg <- acc$pitch+90
    
    #Calculate roll angle
    acc$roll<-atan((-acc$Ys/acc$Zs))
    acc$roll <- acc$roll*(180/pi)
      
    acc$YDabs <- abs(acc$Yd)
    acc$G_tot <- sqrt((acc$X)^2+(acc$Y)^2+(acc$Z)^2)
    
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
      
          ## Descent phase - Accelerometer ---------------------------------------------------------
      descent <- subset(dive_acc, dive_acc$posix>= sealdf$descent_start[x] &
                          dive_acc$posix < sealdf$bottom_start[x])
      descent <- descent[complete.cases(descent$Xs),]
      if(length(descent$ddate)>0){
      n <- as.numeric(length(descent$ddate))
      
      descYDabs <- descent$YDabs
      swimm_eff_df <- swim_eff("descYDabs", dive1 = descent ,n, 0.03, 4)
      sealdf$D_stroke_rate[x] <- swimm_eff_df$stroke_rate
      sealdf$D_swim_eff.g[x] <- swimm_eff_df$swim_eff.g
      sealdf$D_swim_eff.m_s[x] <- swimm_eff_df$swim_eff.m_s
      sealdf$D_mean_pitch[x] <- mean(descent$pitch)
      } else{
        sealdf[x,28:31] <- NA
      }
    
      ## Bottom phase - accelerometer ---------------------------------------------------------------------
      bt <- subset(dive_acc, dive_acc$posix>= sealdf$bottom_start[x] &
                          dive_acc$posix < sealdf$ascent_start[x])
      if(length(bt$ddate)>0){
      #Transmitted PrCA
      prca <- PRCA_transmit("bt", "posix", "G_tot" ,10, 3, 1.5)
      df_prcaT <- as.data.frame(prca[1])
      sealdf$B_PrCA_trans[x] <- df_prcaT$attempts
      time_PrCAT <- as.data.frame(prca[2])
      
      #Archived PrCA
      prca2 <- PRCA_arch("bt", ColXd = bt$Xd, ColYd = bt$Yd, ColZd = bt$Zd, 17, thresholdX, thresholdY,
                         thresholdZ)
      df_prcaA <- as.data.frame(prca2[1])
      sealdf$B_PrCA_arch[x] <- df_prcaA$attempts
      time_PrCAA <- as.data.frame(prca2[2])
      
      #Benthic attempts
      ben <- BENTHIC_prca("bt", bt$pitch_deg, 20,50,30)
      if(length(ben)>1){
        sealdf$B_pitch.diff20[x] <- length(ben$peak_time)
        ben$peak_time <- as.POSIXct(ben$peak_time, format="%Y-%m-%d %H:%M:%S", tz="GMT")
        ben$valley_time <- as.POSIXct(ben$valley_time, format="%Y-%m-%d %H:%M:%S", tz="GMT")
      } else{
        sealdf$B_pitch.diff20[x] <- 0
      }
      
      #Foraging overlap
      if(length(time_PrCAT$get.colNamestr.) > 0 & length(ben)>1){
        time_PrCAT$overlap <- NA
        for(e in 1:length(time_PrCAT$get.colNamestr.)){
          for(t in 1:length(ben$peak_time)){
            if(is.between(time_PrCAT$get.colNamestr.[e], ben$peak_time[t], ben$valley_time[t])){
              time_PrCAT$overlap[e] <- 1
            }
          }
        }
        sealdf$B_overlap_T[x] <- sum(time_PrCAT$overlap, na.rm=TRUE)
      } else {
        sealdf$B_overlap_T[x] <- 0
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
        sealdf$B_overlap_A[x] <- sum(time_PrCAA$overlap, na.rm=TRUE)
      } else{
        sealdf$B_overlap_A[x] <- 0
      }
      
      #Swimming effort
      btn <- as.numeric(length(bt$ddate))
      btYDabs <- bt$YDabs 
      if(is.na(btYDabs)){
        btYDabs <- btYDabs[-which(is.na(btYDabs))]
      }
      swimm_eff_df <- swim_eff("btYDabs",dive1 = bt, btn, 0.03, 4)
      sealdf$B_stroke_rate[x] <- swimm_eff_df$stroke_rate
      sealdf$B_swim_eff.g[x] <- swimm_eff_df$swim_eff.g
      sealdf$B_swim_eff.m_s[x] <- swimm_eff_df$swim_eff.m_s
        
      #Other accelerometer summary values
      sealdf$B_mean_roll[x] <- mean(bt$roll)
      sealdf$B_dyn_sd[x] <- sd(bt$Gd)
      sealdf$B_mean_pitch[x] <- mean(bt$pitch)
      } else{
        sealdf[x,32:43] <- NA
      }

      ## Ascent phase - Accelerometer -----------------------------------------------------------------------
      ascent <- subset(dive_acc, dive_acc$posix>= sealdf$ascent_start[x] &
                           dive_acc$posix <= sealdf$DE_DATE[x])
      if(length(ascent$ddate)>0){
        n <- as.numeric(length(ascent$ddate))
        ascYDabs <- ascent$YDabs
        swimm_eff_df <- swim_eff("ascYDabs", dive1 = ascent ,n, 0.03, 4)
        
        sealdf$A_stroke_rate[x] <- swimm_eff_df$stroke_rate
        sealdf$A_swim_eff.g[x] <- swimm_eff_df$swim_eff.g
        sealdf$A_swim_eff.m_s[x] <- swimm_eff_df$swim_eff.m_s
        sealdf$A_mean_pitch[x] <- mean(ascent$pitch)
      } else{
        sealdf[x,44:47] <- NA
        }
      }else{
        sealdf[x,27:47] <- NA
        }
      }else{
        sealdf[x,27:47] <- NA
      }
    }
    tripsdf <- rbind(tripsdf, sealdf)
    sealdf <- sealdfempty
    print(l)
    } else{
      sealdf <- sealdfempty
      skipped <- c(skipped, trip_seal$Trip_No[l])
      print(l)
      }
    } else{
      sealdf <- sealdfempty
      skipped <- c(skipped, trip_seal$Trip_No[l])
      print(l)
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

write.csv(tripsdf, paste0(here::here("Output", "Processed accelerometer parameters"),"/",PTTacc,"_processed_accelerometer_parameters.csv"), row.names=FALSE)
