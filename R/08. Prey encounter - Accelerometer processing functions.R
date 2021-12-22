### Fuctions created to process the accelerometer data
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose: Functions created to calculate various accelerometer outputs:
#          1. Prey catch attempts using the "archived" method described in Cox et al. 2018
#          2. Find peaks function used the pitch angle analysis for the benthic attempts
#          3. Prey catch attempts using body pitch angle, based on Brasseur et al 2012

#References:
#Brasseur, S., et al. (2012). Habitat Preferences of Harbour Seals in the Dutch Coastal Area: Analysis and Estimate of Effects of Offshore Wind Farms (Report No. OWEZ R 252 T1 20120130 C043-10), IMARES - Wageningen UR, Noordzeewind: 58.
#Cox, S. L., et al. (2018). "Processing of acceleration and dive data on‐board satellite relay tags to investigate diving and foraging behaviour in free‐ranging marine predators." Methods in Ecology and Evolution 9(1): 64-77.

# Output: 
# Created on: 12/01/2021
# Updated on: 13/05/2021

## Empty vector --------------------------------------------------------------------------------------------------
# Find out if a vector is empty 
vector.is.empty <- function(x) return(length(x) ==0)

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