### Extracting accelerometer data
# Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
# Purpose: This code extract from the accelerometer dataset only the times that sel was at sea.
# Output: Individual trips output files with the accelerometer data.
# Created on: 14/01/2021
# Updated on: 13/05/2021

## Load packages --------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, magrittr)


# Store individual trip subsets
#Seals with accelerometer data
IDs <- c(90,242,285,158,283)
PTT <- c(14438,14464,14477,14478,14479)

#CSV with trip data
all_trips <- read.csv(here::here("Datasets", "pv64-2017_trip_summaries.csv"), header=TRUE) %>%
  mutate(Trip_Start = as.POSIXct(Trip_Start , format="%Y-%m-%d %H:%M:%S", tz="UTC"),
         Trip_End = as.POSIXct(Trip_End , format="%Y-%m-%d %H:%M:%S", tz="UTC"))

#Due to the dimension of the data they need to be loaded a seal at the time
#This code used the accelerometer data stored as rdata, needs to be moodified with any other format
for(i in 1:length(PTT)){
  PTTseal = PTT[i]
  load(paste0("F:/pv64_tag_archived_accelerometer_and_temperature_data/",PTTseal,"/",PTTseal,"_all_processed_StatDyn.rdata"))
  
  trip_seal <- all_trips[which(all_trips$PTT %in% PTTseal),]
  
  for(x in 1:length(trip_seal$ID)){
    dates <- seq(as.Date(trip_seal$Trip_Start[x]), as.Date(trip_seal$Trip_End[x]), by="days")
    dates <- format(dates, "%Y/%m/%d")
    tmp <- subset(pvac, DATE %in% dates)
    if(length(tmp$ddate > 0)){
      tmp<- subset(tmp, tmp$ddate>= trip_seal$Trip_Start[x] &
                     tmp$ddate<= trip_seal$Trip_End[x])
      tmp$seconds <- format(tmp$ddate, "%Y-%m-%d %H:%M:%OS6")
      write.csv(tmp,paste0(here::here("Output","Raw accelerometer data per trip"),"/",PTTseal,"/pv",PTTseal,"_",trip_seal$Trip_No[x],"_trip_data.csv"), row.names=FALSE)
    } 
    tmp <- NA
    print(x)
  }
}

