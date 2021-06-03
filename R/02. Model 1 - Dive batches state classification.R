## Assign HMM state to dive batches ------------------------------------------------------------
##Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
##Purpose: The code uses an Hidden Markov Model to classify harbour seal dive batches into transit and ARS state.
##Output: 1. Model 1 - dataset.csv
#         2. Model 1 - initial parameters selection output.csv
#         3. Model 1 - HMM dive batches classified.csv
#Created on: 06/01/2021
#Updated on: 12/01/2021

#Load necessary package and create these functions
library(pacman)
p_load(momentuHMM, ggplot2, magrittr, tidyverse)

is.between<-function(x, a, b) {
  (x >= a) & (b >= x)
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

# Data preparation -----------------------------------------------------------------------------------------
##Read in dataset
df <- read.csv(here::here("Output", "Dive batches dataset - 2017.csv"), header=TRUE) 

df <- df %>% mutate(
    ID = format(ID, nsmall=3),
    start.time = as.POSIXct(start.time, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    end.time = as.POSIXct(end.time, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    batch.duration = as.numeric(difftime(end.time, start.time, units = "mins")))

tmp <- df %>% group_by(df$ID) %>% summarise(n=n())
short.trips <- tmp$`df$ID`[which(tmp$n<3)]
df <- df[which(!df$ID %in% short.trips),]

##Prepare data for HMM
df.HMM <- prepData(df, type="UTM", coordNames=c("x", "y"))

##Flag uncertain locations
#Flag interpolated locations that don't have at least 1 GPS locations within a 30 minutes window
gps <- read.table(here::here("Dryad", "pv64-2017_gps_data_with_haulout_&_trip_info.txt"), sep="\t",header=TRUE)
gps$time <- as.POSIXct(gps$D_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC")
gps$trip_code <- format(gps$trip_code, nsmal=3)

#Select only the trips in the analysis
trips <- unique(df.HMM$ID)
gps <- gps[gps$trip_code %in% trips,]

#Create a column by which merge them
flag.tmp <- df.HMM[1,]
flag.tmp[1,] <- NA
flag.tmp$flag <- NA

#Flagging for batch data
i <-  1
x <-1 
for(i in 1:length(trips)){
  HMM_tmp <- subset(df.HMM, df.HMM$ID==trips[i])
  gps_tmp <- subset(gps, gps$trip_code==trips[i])
  
  for(n in 1:length(HMM_tmp$ID)){
    m <- as.integer(which(is.between(gps_tmp$time, HMM_tmp$start.time[n], HMM_tmp$end.time[n])))
    if(is.integer0(m)){
      if(HMM_tmp$batch.duration[n]<=25){
        p <- as.integer(which(is.between(gps_tmp$time, HMM_tmp$start.time[n-1], HMM_tmp$end.time[n-1])))
        o <- as.integer(which(is.between(gps_tmp$time, HMM_tmp$start.time[n+1], HMM_tmp$end.time[n+1])))
        HMM_tmp$flag[n] <- ifelse(is.integer0(p) & is.integer0(o), 1 ,0)
      } else {
        HMM_tmp$flag[n] <- 1
      }
    } else {
      HMM_tmp$flag[n] <- 0
    }
  }
  flag.tmp <- rbind(flag.tmp, HMM_tmp)
}

flag.tmp <- flag.tmp[-1,]
df.HMM$flag <- flag.tmp$flag
df.HMM$step <- ifelse(df.HMM$flag==1, NA, df.HMM$step)
df.HMM$angle <- ifelse(df.HMM$flag==1, NA, df.HMM$angle)

#Remove ipossible step lenghts 
#Moving at 1.5 m/s in 30 minutes a marine mammal cannot move further than 2700 m
quant99.duration <- as.numeric(quantile(df.HMM$batch.duration, probs=.99))
max.distance <- quant99.duration*60*1.5

df.HMM$step <- ifelse(df.HMM$step>max.distance, NA, df.HMM$step)
df.HMM$angle <- ifelse(df.HMM$step>max.distance, NA, df.HMM$angle)

hist(df.HMM$step)
hist(df.HMM$angle)

#Save dataset for model
write.csv(df.HMM, here::here("Output","Model 1 - dataset.csv"), row.names=TRUE)

# Select model initial values ===================================================================================
data = df.HMM
m_list<-list()
n_its<-50
output<-data.frame(iter<-seq(1,n_its), s1_mean = NA, s2_mean = NA,
                   s1_sd = NA, s2_sd = NA, s1_zero= NA, s2_zero =NA,
                   s1_angle= NA, s2_angle = NA, AIC = NA, loglik = NA)
stateNames <- c("state1","state2")
dist <- list(step="gamma", angle="wrpcauchy")
i <- 1
for(i in 1:n_its){
  output$s1_mean[i]<-runif(1,700,1500)
  output$s2_mean[i]<-runif(1,50,500)
  output$s1_sd[i]<-runif(1,100,200)
  output$s2_sd[i]<-runif(1,50,100)
  output$s1_zero[i]<-runif(1,0,1)
  output$s2_zero[i]<-runif(1,0,1)
  output$s1_angle[i]<-runif(1,0,0.5)
  output$s2_angle[i]<-runif(1,0.5,1)
  
  stepPar0<-c(output$s1_mean[i],output$s2_mean[i],
              output$s1_sd[i],output$s2_sd[i],
              output$s1_zero[i],output$s2_zero[i])
  anglePar0<-c(output$s1_angle[i],output$s2_angle[i])
  
  m_list[[i]]<-try(fitHMM(data=data, nbStates=2, dist=dist, 
                          Par0=list(step=stepPar0[1:4], angle=anglePar0),
                          stateNames=stateNames),silent=TRUE)
  
  try(output$AIC[i]<-AIC(m_list[[i]]),silent=TRUE)
  try(output$loglik[i]<-m_list[[i]]$mod$minimum,silent=TRUE)
  print(i)
}
plot(output$iter....seq.1..n_its., output$AIC)

#Include the saving output
write.csv(output, here::here("Output","Model 1 - initial parameters selection output.csv"), row.names=TRUE)

# Running the model ###########################################################################################
data=df.HMM

nbstates <- 2
stateNames <- c("state1","state2")
dist <- list(step="gamma", angle="wrpcauchy")
b <- as.numeric(which.min(output$AIC))

stepPar0<-c(output$s1_mean[b],output$s2_mean[b],
            output$s1_sd[b],output$s2_sd[b],
            output$s1_zero[b],output$s2_zero[b])
anglePar0<-c(output$s1_angle[b],output$s2_angle[b])

#Initial parameters of the best model in the iterations
# stepPar0 <- c(1390.76332208,131.2383175,190.18469572,91.55341004,0.44187353,0.08603832)
# anglePar0 <- c(0.17120304,0.79857152)

m1 <- fitHMM(data=data, nbStates= nbstates, dist=dist, 
             Par0=list(step=stepPar0[1:4], angle=anglePar0),
             stateNames=stateNames)
m1
#Copy m1 output into a text file
-m1$mod$minimum
AIC(m1)
-2*(-m1$mod$minimum)+length(m1$mod$wpar)*log(length(data$ID))
plot(m1, plotTracks = FALSE)
plotPR(m1)
#state2 <- ARS
df.HMM$HMMstate <- viterbi(m1) 

df.HMM$HMMstate <- ifelse(is.na(df.HMM$step), NA, df.HMM$HMMstate)
df.HMM$state <- ifelse(df.HMM$HMMstate==2, "ARS", "Transit")

write.csv(df.HMM, here::here("Output", "Model 1 - HMM dive batches classified.csv"), row.names = TRUE)
