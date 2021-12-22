### Acceleration PrCA threshold eStimation
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose:This code is used to calculate acceleration axes thresholds for all seals using the accelerometer. data 
# Output: One summary table with a line for each seal and the respective threshold for each axes
# Created on: 13/01/2021
# Updated on: 13/05/2021

## Load packages --------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, magrittr, roll)

# Calculate rolling standard deviation of dynamic acceleration -------------------------------------------------
#Due to the huge amount of data it is better to load one trip at the time to calculate the rolling standard deviation

df <- data.frame(ID = c(90,242,283,158,283), PTT = c(14438,14464,14477,14478,14479),
                 ThresholdX = NA, ThresholdY = NA, ThresholdZ = NA)
for (x in 1:length(df$ID)){
  PTT = df$PTT[x]
  files <- list.files(paste0(here::here("Dryad","Outputs","Raw accelerometer data per trip"),"/",PTT))

  clusterdf <- data.frame(Xsd = NA, Ysd = NA, Zsd = NA)

  for(l in 1:length(files)){
    file <- read.table(paste0(here::here("Dryad","Outputs","Raw accelerometer data per trip"),"/",PTT,"/",files[l]), 
                       sep="\t", header=TRUE)
    Xsd <- as.matrix(file$Xd)
    file$Xsd  <- roll_sd(Xsd, 19)
    Ysd <- as.matrix(file$Yd)
    file$Ysd  <- roll_sd(Ysd, 19)
    Zsd <- as.matrix(file$Zd)
    file$Zsd  <- roll_sd(Zsd, 19)
    file2 <- file[!is.na(file$Xsd), ]
    tmp <- select(file2, c("Xsd","Ysd","Zsd"))
    clusterdf <- rbind(clusterdf, tmp)
    print(l)
    }
  clusterdf <- clusterdf[-1,]

  # Cluster analysis ------------------------------------------------------------------------
  clusterdf <- clusterdf[!is.na(clusterdf$Xsd), ]
  acc_clusterX <- kmeans(clusterdf$Xsd, 2)
  acc_clusterY <- kmeans(clusterdf$Ysd, 2)
  acc_clusterZ <- kmeans(clusterdf$Zsd, 2)
  

  #Select threshold
  df$ThresholdX[x] <- max(acc_clusterX$centers)
  df$ThresholdY[x] <- max(acc_clusterY$centers)
  df$ThresholdZ[x] <- max(acc_clusterZ$centers)
}

write.table(df, here::here("Dryad","Outputs","Prey encounters - cluster analysis axis thresholds.csv"), 
            sep="\t" ,row.names = FALSE)
