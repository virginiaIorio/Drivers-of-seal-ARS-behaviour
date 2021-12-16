## Assign HMM state to dive batches and test covariates effect 
##Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
##Purpose: The code uses an Hidden Markov Model to classify harbour seal dive batches into transit and ARS state.
#           I also assess the influence of spatial memory and prey encounter on the transition probability.
##Output: Dataset with step length and turning angle, output of the iteration with the selection of 
#         initial parameters, output of the covariates selection process, final dataset with the HMM state classification
#Created on: 04/12/2020
#Updated on: 13/01/2021

library(pacman)
p_load(sp, sf, magrittr, lubridate, tidyverse, geosphere, momentuHMM, ggplot2, raster)

is.between<-function(x, a, b) {
  (x >= a) & (b >= x)
}
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

## Load memory raster maps ------------------------------------------------------------------------
raster_path <- here::here("Output","Memory_grids")
all_rasters <- list.files(raster_path,
                          full.names = TRUE,
                          pattern = ".tif$")
stack <- stack(all_rasters)
Brick <- brick(stack)

#IDs <- c(242,90,285,158,283)
#Brick@data@names <- as.character(IDs)
Brick <- setZ(Brick, Brick@data@names, name='seal_ID')
plot(Brick$seal_158)

Brick[is.na(Brick[])] <- 0


## Data preparation -----------------------------------------------------------------------------------------
##Read in dataset
df <- read.csv(here::here("Output", "Dive batches dataset - Accelerometer.csv"), header=TRUE) 
df <- df %>% mutate(
  ID = format(ID, nsmall=3),
  start.time = as.POSIXct(start.time, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
  end.time = as.POSIXct(end.time, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
  batch.duration = as.numeric(difftime(end.time, start.time, units = "mins")))

tmp <- df %>% group_by(df$ID) %>% summarise(n=n())
short.trips <- tmp$`df$ID`[which(tmp$n<3)]
df <- df[which(!df$ID %in% short.trips),]

df$seal_ID <- paste0("seal_",df$seal_ID)

##Prepare data for HMM
seal.HMM <- prepData(df, type="UTM", coordNames=c("x", "y"), spatialCovs = list(memory=Brick))

##Flag uncertain locations
gps <- read.table(here::here("Dryad", "pv64-2017_gps_data_with_haulout_&_trip_info.txt"),sep="\t" ,header=TRUE)
gps$time <- as.POSIXct(gps$D_DATE, format="%Y-%m-%d %H:%M:%S", tz="UTC")
gps$trip_code <- format(gps$trip_code, nsmal=3)

#Select only the trips in the analysis
seal.HMM$ID <- format(seal.HMM$ID, nsmall=3)
trips <- unique(seal.HMM$ID)
gps <- gps[gps$trip_code %in% trips,]

#Create a column by which merge them
flag.tmp <- seal.HMM[1,]
flag.tmp[1,] <- NA
flag.tmp$flag <- NA

#Flagging for batch data
for(i in 1:length(trips)){
  HMM_tmp <- subset(seal.HMM, seal.HMM$ID==trips[i])
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
seal.HMM$flag <- flag.tmp$flag
seal.HMM$step <- ifelse(seal.HMM$flag==1, NA, seal.HMM$step)
seal.HMM$angle <- ifelse(seal.HMM$flag==1, NA, seal.HMM$angle)


#Set also to NA the unbelivable step lenghts
quant99.duration <- as.numeric(quantile(seal.HMM$batch.duration, probs=.99))
max.distance <- quant99.duration*60*1.5

seal.HMM$step <- ifelse(seal.HMM$step>max.distance, NA, seal.HMM$step)
seal.HMM$angle <- ifelse(seal.HMM$step>max.distance, NA, seal.HMM$angle)

hist(seal.HMM$step)
hist(seal.HMM$angle)

#Save dataset for model
write.csv(seal.HMM, here::here("Output","Model 3 - dataset.csv"), row.names=TRUE)


## Select model initial parameters ------------------------------------------------------------------------------------
data = seal.HMM
m_list<-list()
n_its<-50
output<-data.frame(iter<-seq(1,n_its), s1_mean = NA, s2_mean = NA,
                   s1_sd = NA, s2_sd = NA, s1_zero= NA, s2_zero =NA,
                   s1_angle= NA, s2_angle = NA, AIC = NA, loglik = NA)
stateNames <- c("state1","state2")
dist <- list(step="gamma", angle="wrpcauchy")
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
plot(output$iter....seq.1..n_its., output$loglik)

#Include the saving output
write.csv(output, here::here("Output","Model 3 - initial parameters selection output.csv"), row.names=TRUE)
output <- read.csv(here::here("Output","Model 3 - initial parameters selection output.csv"))

## Run simple model --------------------------------------------------------------------------------------------------------
data=seal.HMM

nbstates <- 2
stateNames <- c("state1","state2")
dist <- list(step="gamma", angle="wrpcauchy")
b <- as.numeric(which.min(output$AIC))

stepPar0<-c(output$s1_mean[b],output$s2_mean[b],
            output$s1_sd[b],output$s2_sd[b],
            output$s1_zero[b],output$s2_zero[b])
anglePar0<-c(output$s1_angle[b],output$s2_angle[b])

ms <- fitHMM(data=data, nbStates= nbstates, dist=dist, 
             Par0=list(step=stepPar0[1:4], angle=anglePar0),
             stateNames=stateNames)

-ms$mod$minimum
AIC(ms)
-2*(-ms$mod$minimum)+length(ms$mod$wpar)*log(length(data$ID))
plot(ms, plotTracks = FALSE)

## Model selection with covariates -------------------------------------------------------------------------------------
formulas.list <- c(~ batch.foraging.index + memory,
                   ~ batch.foraging.index,
                   ~ memory)

model.selection <- data.frame("covariate formula" = c("null model", "~ batch.foraging.index + memory",
                                                     "~ batch.foraging.index", "~ memory"), 
                              "log-likelihood" = c(-ms$mod$minimum,NA,NA,NA),
                              "AIC" = c(AIC(ms),NA,NA,NA), 
                              "BIC" = c(-2*(-ms$mod$minimum)+length(ms$mod$wpar)*log(length(data$ID)),NA,NA,NA))

for(x in 1:length(formulas.list)){
  formula <- formulas.list[[x]]
  Par0 <- getPar0(ms, nbStates = nbstates, formula=formula)
  m <- fitHMM(data=data, nbStates= nbstates, dist=dist, 
              Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
              stateNames=stateNames, formula= formula) 
  
  model.selection$log.likelihood[x+1] <- -m$mod$minimum
  model.selection$AIC[x+1] <- AIC(m)
  model.selection$BIC[x+1] <- -2*(-m$mod$minimum)+length(m$mod$wpar)*log(length(data$ID))
  print(x)
}

n <- which.min(model.selection$AIC)
model.selection$delta.AIC <- model.selection$AIC - model.selection$AIC[n]
n <- which.min(model.selection$BIC)
model.selection$delta.BIC <- model.selection$BIC - model.selection$BIC[n]

write.csv(model.selection, here::here("Output", "Model 3 - Covariates model selection.csv"), row.names=FALSE)

## Model 2 - drivers of ARS ------------------------------------------------------------------------------------------
formula <- ~ batch.foraging.index+memory
Par0 <- getPar0(ms, nbStates = nbstates, formula=formula)
m3 <- fitHMM(data=data, nbStates= nbstates, dist=dist, 
            Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
            stateNames=stateNames, formula= formula) 

-m3$mod$minimum #Log-likelihood
AIC(m3) #AIC
-2*(-m3$mod$minimum)+length(m3$mod$wpar)*log(length(data$ID)) #BIC

plot(m3, plotCI=TRUE, plotTracks=FALSE)
plotPR(m3)
plotStationary(m3, plotCI=TRUE)

plot(m3, plotCI=TRUE, plotTracks=FALSE)

seal.HMM$HMMstate <- viterbi(m3)
seal.HMM$HMMstate <- ifelse(seal.HMM$flag==1, NA, seal.HMM$HMMstate)

step.state1 <- m3[["mle"]][["step"]][1]
step.state2 <- m3[["mle"]][["step"]][3]

step.ARS <- as.numeric(which.min(c(step.state1, step.state2)))

seal.HMM$state <- ifelse(seal.HMM$HMMstate==step.ARS, "ARS", "Transit")

write.csv(seal.HMM, here::here("Output", "Model 3 - HMM dive batches classified.csv"), row.names = TRUE)
