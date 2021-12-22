## Create memory raster 
## Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
## Purpose: Create seals individual memory raster map to use as covariates in Model 3 for the 5 seals for which we have
##          accelerometer data.
## Output: Individual raster files
#Created on: 04/12/2020
#Updated on: 13/01/2020

#Load necessary packages
library(pacman)
p_load(sp,rgdal, sf, tidyverse, fasterize, raster, lubridate, ggplot2, geosphere)

## Data preparation ------------------------------------------------------------------------------------
mem <- read.table(here::here("Dryad","Outputs", "Model 1 - HMM dive batches classified.txt"), sep="\t", header=TRUE)
seal.IDs <- c(90,158,242,283,285)
mem <- mem[which(mem$seal_ID %in% seal.IDs),]

x <- 1
for(x in 1:length(mem$ID)){
  mid.point <- midPoint(c(mem$batch.start.lon[x], mem$batch.start.lat[x]), c(mem$batch.end.lon[x], mem$batch.end.lat[x])) 
  mem$batch.mid.point.long[x] <- mid.point[1]
  mem$batch.mid.point.lat[x] <- mid.point[2]
}

crs <- "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

llcord <- SpatialPoints(mem[,c("batch.mid.point.long","batch.mid.point.lat")],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
mem$x <- attr(utmcoord, "coords")[,1]
mem$y <- attr(utmcoord, "coords")[,2]
maxx <- max(mem$x)
mem_sf <- st_as_sf(mem, coords=c("x", "y"), crs=crs)

## Create empty raster -----------------------------------------------------------------------------------
grid <- st_read(here::here("Dryad","Moray_Firth_1km_grid_shapefile","MF_grid_1km_UTM30.shp"))
grid <- st_transform(grid, crs=crs)

#Create raster size based on the data to analyse
grid_crop <- grid[which(grid$grid_id>=2080 & grid$grid_id<14000),]

grid.matrix <- matrix(unlist(grid_crop$geometry), ncol=10, byrow=TRUE)
grid.matrix <- as.data.frame(grid.matrix)
grid.matrix$unique.id <- grid_crop$grid_id
cells.to.keep <- grid.matrix$unique.id[which(grid.matrix[,c(1:5)]< maxx+10000)]

grid_crop2 <- grid_crop[which(grid_crop$grid_id %in% cells.to.keep),]
grid.matrix <- matrix(unlist(grid_crop2$geometry), ncol=10, byrow=TRUE)
grid.matrix <- as.data.frame(grid.matrix)

ext <- extent(min(grid.matrix[,c(1:5)]), max(grid.matrix[,c(1:5)]), min(grid.matrix[,c(6:10)]), max(grid.matrix[,c(6:10)]))
gridsize <- 1000 #meters
r <- raster(ext, res=gridsize)

## Select data ----------------------------------------------------------------------------------------
#from a month before the date in which acc start. Manually specify the start of the beginning of the accelerometer data
acc.start.date <- data.frame(ID = c(90,158,242,283,285), acc_start= c("15/05/2017","28/04/2017","23/05/2017","31/05/2017", "16/04/2017")) %>%
  mutate(
    acc_start = as.Date(acc_start, format="%d/%m/%Y"),
    one_month_before = acc_start)
month(acc.start.date$one_month_before) <- month(acc.start.date$acc_start) - 1
acc.start.date$one_month_before[4] <-"2017-04-30"

mem_sf$date <- as.Date(mem_sf$start.time)
dat <- mem_sf[1,]
x <- 1
for(x in 1:length(acc.start.date$ID)){
  trips <- unique(mem_sf$ID[which(mem_sf$seal_ID == acc.start.date$ID[x] & 
                                  mem_sf$date>=acc.start.date$one_month_before[x] & 
                                  mem_sf$date<=acc.start.date$acc_start[x])])
  tmp <- mem_sf[which(mem_sf$ID %in% trips),]
  dat <- rbind(dat, tmp)
}
dat <- dat[-1,]

dat2 <- dat
st_geometry(dat2) <- NULL
write.table(dat2, here::here("Dryad","Outputs","Spatial memory - dataset for rasters - accelerometer seals.txt"), 
            sep="\t", row.names = FALSE)

## Fill empty raster -----------------------------------------------------------------------------------
#Calculate number of transiting and ARS dive batches occurring in each grid cells for each seal
mem_sftmp <- st_join(dat, grid_crop2)
cells <- mem_sftmp %>% 
  group_by(grid_id, seal_ID) %>% 
  summarise(
    count=n(), 
    for.count=length(which(state == "ARS")), #Dive batches as ARS
    trans.count = length(which(state == "Transit")), #Dive batches as Transit
    n.trips= length(unique(ID))
    )

IDs <- c(242,90,285,158,283)

#Create raster brick
rasterIds <- unique(cells$grid_id)
grid_crop2$effort <- NA
y <- 2
for(y in 1:length(IDs)){
  values <- cells[which(cells$seal_ID==IDs[y]),]
  values <- as.data.frame(values)
  
  values$for.prop <- values$for.count/values$count
  
  mu <- mean(values$for.prop)
  
  rasterIds <- unique(values$grid_id)
  for(i in 1:length(grid_crop$grid_id)){
    if(grid_crop$grid_id[i] %in% rasterIds){
      x <- which(values$grid_id == grid_crop$grid_id[i]) 
      grid_crop$effort[i] <- values$for.prop[x]
    } else {
      grid_crop$effort[i] <- mu
    }
  }
  
  rr <- fasterize(grid_crop, r, field="effort")
  names(rr) <- paste0("seal ",IDs[y])
  plot(rr)

  outfile <- writeRaster(rr, filename=here::here("Dryad","Outputs","Accelerometer memory raster maps",paste0( paste0("seal ",IDs[y]),".tiff")),
                         format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"), bylayer=TRUE)
}
