### Fuctions created to process the accelerometer data
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose: Functions created to calculate various accelerometer outputs:
#          1. Swimming effort as the absolute acceleration on the Y axes
#          2. Prey catch attempts using the "transmitted" method described in Cox et al. 2018
#          3. Prey catch attempts using the "archived" method described in Cox et al. 2018
#          4. Find peaks function used the pitch angle analysis for the benthic attempts
#          5. Prey catch attempts using body pitch angle, based on Brasseur et al 2012

#References:
#Brasseur, S., et al. (2012). Habitat Preferences of Harbour Seals in the Dutch Coastal Area: Analysis and Estimate of Effects of Offshore Wind Farms (Report No. OWEZ R 252 T1 20120130 C043-10), IMARES - Wageningen UR, Noordzeewind: 58.
#Cox, S. L., et al. (2018). "Processing of acceleration and dive data on‐board satellite relay tags to investigate diving and foraging behaviour in free‐ranging marine predators." Methods in Ecology and Evolution 9(1): 64-77.

# Output: 
# Created on: 12/01/2021
# Updated on: 13/05/2021

## Empty vector --------------------------------------------------------------------------------------------------
# Find out if a vector is empty 
vector.is.empty <- function(x) return(length(x) ==0)


# Calculate swimming effort--------------------------------------------------------------------------------------------------
#a=column with Yabsolute values, n= length of the column, e= minpeakheight, f= minpeakdistance
swim_eff = function(colNamestr,dive1,n,peakheight,peakdistance){
  peaks <- findpeaks(get(colNamestr), minpeakheight = peakheight, minpeakdistance = peakdistance) 
  if(vector.is.empty(peaks)) {
    stroke_rate <- NA
    start_dive <- dive1$posix[1]
    end_dive <- dive1$posix[length(dive1$DATE)]
    t <- as.numeric(difftime(end_dive, start_dive, unit="secs"))
    dive1$m_s <- dive1$YDabs*9.80665
    swim_eff.g <- sum(dive1$YDabs)/t
    swim_eff.m_s <- sum(dive1$m_s)/t
    output <- data.frame("stroke_rate" = stroke_rate, "swim_eff.g" = swim_eff.g, 
                         "swim_eff.m_s" = swim_eff.m_s)
  } else{
    peaks <- as.data.frame(peaks)
    peaks_index <- peaks$V2
    dive1$row.num <- seq(1, n, 1)
    for(i in 1:n){
      if(dive1$row.num[i] %in% peaks_index){
        dive1$peak[i] <- 1
      } else{dive1$peak[i] <- 0}
    }
    start_dive <- dive1$posix[1]
    end_dive <- dive1$posix[length(dive1$DATE)]
    t <- as.numeric(difftime(end_dive, start_dive, unit="secs"))
    eff <- filter(dive1, peak> 0) #Threhsold from Heerah 0.2 m/s = 0.02 g
    stroke_rate <- length(eff$DATE)/t
    dive1$m_s <- dive1$YDabs*9.80665
    swim_eff.g <- sum(dive1$YDabs)/t
    swim_eff.m_s <- sum(dive1$m_s)/t
    output <- data.frame("stroke_rate" = stroke_rate, "swim_eff.g" = swim_eff.g, 
                         "swim_eff.m_s" = swim_eff.m_s)
    return(output)
  }
}



## Transmitted PrCA --------------------------------------------------------------------------------------------------
# Calculate Prey catch attempts using the simplified method of Cox et al.
varS_fun = function(a,i){
  sum(abs(a[i]-a[i-1]), na.rm=TRUE)
}

PRCA_transmit = function(dfNameStr, colNamestr, magA, window, observation, threshold_g){
  require(dplyr)
  require(magrittr)
  df <- get(dfNameStr)
  dat <- df %>% group_by(get(colNamestr)) %>% 
    summarise(varS= varS_fun(get(magA),(2:n())))
  dat1 <- as.matrix(dat$varS)
  dat$varA <- roll_mean(dat1, width = window, min_obs = observation)
  dat$diff <- dat$varS-dat$varA
  dat$PrCA <- ifelse(dat$diff >= threshold_g, 1, 0)
  dat <- dat[!is.na(dat$varS), ]
  
  attempts = sum(dat$PrCA, na.rm = TRUE)
  df_Prca <- data.frame("attempts" = attempts)
  time_PrCA <- subset(dat, dat$PrCA==1)
  return(list(df_Prca, time_PrCA))
  return(df_Prca)
}

## Archived PrCA --------------------------------------------------------------------------------------------------
# Calculate Prey catch attempts using the archival method of Cox et al.
PRCA_arch = function(dfNameStr, ColXd, ColYd, ColZd, window, thresholdX, thresholdY, thresholdZ){
  require(roll)
  df <- get(dfNameStr)
  df$highX <- ifelse(df$Xsd>=thresholdX, 1, 0)
  df$highY <- ifelse(df$Ysd>=thresholdY, 1, 0)
  df$highZ <- ifelse(df$Zsd>=thresholdZ, 1, 0)
  df$highstate <-  df$highX+df$highY+df$highZ
  
  HIGH <- df[which(df$highstate==3),]
  LEN <- length(which(df$highstate==3))
  if(LEN>0){
    HIGH$diff <- NA
    for(w in 1:LEN){
      HIGH$diff[w] <- as.numeric(difftime(HIGH$posix_sec[w+1],HIGH$posix_sec[w], unit="sec"))
    }
    attempts = length(which(HIGH$diff>=1))+1
    df_Prca2 <- data.frame("attempts" = attempts)
    
    row_end <- which(HIGH$diff>=1)
    row_start <- row_end+1
    row_start <- c(1, row_start)
    row_end <- c(row_end, length(HIGH$ddate))
    time_PrCA <- data.frame(start = HIGH$posix_sec[row_start], end = HIGH$posix_sec[row_end])
    
    return(list(df_Prca2, time_PrCA))
    return(df_Prca2)
  } else{
    attempts = 0
    df_Prca2 <- data.frame("attempts" = attempts)
    time_PrCA <- NA
    return(list(df_Prca2, time_PrCA))
    return(df_Prca2)
  }
}

## Find peaks function --------------------------------------------------------------------------------------------------
#https://github.com/stas-g/findPeaks
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
} 

## Benthic attempts --------------------------------------------------------------------------------------------------
#Function to calculate benthic attempts
BENTHIC_prca = function(dfNameStr, col_pitch_deg, degree_thresh, m_peaks , m_valley){
  df <- get(dfNameStr)
  peaks <- as.data.frame(find_peaks(col_pitch_deg, m=  m_peaks))
  valley <- as.data.frame(find_peaks(-col_pitch_deg, m= m_valley))
  if(is.na(peaks[1,1]) || is.na(valley[1,1]) || vector.is.empty(peaks) || vector.is.empty(valley)){
    bt_preydf <- NA
    return(bt_preydf)
  } else {
    colnames(peaks) <- c("peak")
    colnames(valley) <- c("peak")
    peaks$lab <- rep("peak", length(peaks))
    valley$lab <- rep("valley", length(valley))
    pv <- rbind(peaks, valley)
    pv <- pv[order(pv$peak),]
    LEN<- length(df$X)
    df$row.num <- seq(1, LEN, 1)
    dftmp <- merge(df,pv, by.x="row.num", by.y="peak", all.x=TRUE)
    dftmp <- dftmp[complete.cases(dftmp$lab),]
    dftmp <- spread(dftmp, lab, pitch_deg)
    if(is.na(dftmp$peak[1])){dftmp$peak[1] <- 0}else{}
    dftmp$start_peak <- ifelse(is.na(dftmp$peak), NA, paste(dftmp$seconds))
    dftmp$row_peak <- ifelse(is.na(dftmp$peak), NA, paste(dftmp$row.num))
    dftmp$peak <- fill.na(dftmp$peak)
    dftmp$start_peak <- fill.na(dftmp$start_peak)
    dftmp$row_peak <- fill.na(dftmp$row_peak)
    dftmp <- dftmp[complete.cases(dftmp$valley),]
    dftmp$diff <- dftmp$peak - dftmp$valley
    dftmp$time_diff <- difftime(dftmp$seconds, dftmp$start_peak, units="secs")
    
    dftmp$attempts <- ifelse(dftmp$diff>=degree_thresh & dftmp$time_diff<=5, 1, 0)
    attempts <- dftmp[which(dftmp$attempts==1),]
    attempts <- attempts[!duplicated(attempts$row_peak),]
    if(length(attempts$row.num)>0){
      bt_preydf <- data.frame(peak_time = attempts$start_peak, 
                              valley_time = attempts$seconds, drop_duration = attempts$time_diff)
      return(bt_preydf)
    } else{
      bt_preydf <- NA
      return(bt_preydf)
    }
  }
}

## is.between --------------------------------------------------------------------------------------------------
# Function to find out if a time is between a specific window
is.between<-function(x, a, b) {
  (x >= a) & (b >= x)
}

## is.integer0 --------------------------------------------------------------------------------------------------
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}