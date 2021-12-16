## Create memory raster 
## Author: Virginia Iorio (v.iorio1.18@abdn.ac.uk)
## Purpose: Create seals individual memory raster map to use as covariates in Model 2. Representing the areas they 
##          searched in during the month of May.
## Output: Individual raster files
#Created on: 04/12/2020
#Updated on: 16/12/2021

#Load necessary packages
library(pacman)
p_load(sp,rgdal, sf, tidyverse, fasterize, raster, lubridate, ggplot2, geosphere)

## Data preparation ------------------------------------------------------------------------------------
mem <- read.csv(here::here("Output", "Model 1 - HMM dive batches classified.csv"), header=TRUE)

mem <- mem[,-1]
x <- 1
#Calculate mid point of each batch 
for(x in 1:length(mem$ID)){
  mid.point <- midPoint(c(mem$batch.start.lon[x], mem$batch.start.lat[x]), c(mem$batch.end.long[x], mem$batch.end.lat[x])) 
  mem$batch.mid.point.long[x] <- mid.point[1]
  mem$batch.mid.point.lat[x] <- mid.point[2]
}

crs <- "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

llcord <- SpatialPoints(mem[,c("batch.start.lon","batch.start.lat")],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
mem$x <- attr(utmcoord, "coords")[,1]
mem$y <- attr(utmcoord, "coords")[,2]
maxx <- max(mem$x)
mem_sf <- st_as_sf(mem, coords=c("x", "y"), crs=crs)

## Create empty raster -----------------------------------------------------------------------------------
grid <- st_read(here::here("Datasets","Moray_Firth_1km_grid_shapefile","MF_grid_1km_UTM30.shp"))
grid <- st_transform(grid, crs=crs)

grid.matrix <- matrix(unlist(grid$geometry), ncol=10, byrow=TRUE)
grid.matrix <- as.data.frame(grid.matrix)
grid.matrix$unique.id <- grid$grid_id

ext <- extent(min(grid.matrix[,c(1:5)]), max(grid.matrix[,c(1:5)]), min(grid.matrix[,c(6:10)]), max(grid.matrix[,c(6:10)]))
gridsize <- 1000 #meters
r <- raster(ext, res=gridsize)

## Select data ----------------------------------------------------------------------------------------
# Select trips that start in April
trip <- read.table(here::here("Dryad","pv64-2017_trip_summaries.txt"),sep="\t", head=TRUE)
trip$month <- lubridate::month(trip$Trip_Start)
trip$Trip_Code <- format(trip$Trip_Code, nsmall=3)

April.trips <- trip$Trip_Code[which(trip$month==4)]

mem_sf$ID <- format(mem_sf$ID, nsmall=3)

mem_apr <- mem_sf[which(mem_sf$ID %in% April.trips),]

write.csv(mem_apr, here::here("Output","Spatial memory - dataset for rasters all seals.csv"), row.names = FALSE)

## Fill empty raster -----------------------------------------------------------------------------------
#Calculate number of transiting and ARS dive batches occurring in each grid cells for each seal

IDs <- unique(mem_apr$seal_ID)
#Create raster brick
rasterIds <- unique(grid$grid_id)
grid$effort <- NA
y <- 1
for(y in 1:length(IDs)){
  seal_mem <- mem_apr[which(mem_apr$seal_ID==IDs[y]),]
  
  mem_tmp <- st_join(seal_mem, grid)
  cells <- mem_tmp %>% 
    group_by(grid_id) %>% 
    summarise(
      count=n(), 
      for.count=length(which(state == "ARS")), #Dive batches as ARS
      trans.count = length(which(state == "Transit")), #Dive batches as Transit
      n.trips= length(unique(ID))
    )
  
  cells$for.prop <- cells$for.count/cells$count
  
  mu <- mean(cells$for.prop)
  
  rasterIds <- unique(cells$grid_id)
  for(i in 1:length(grid$grid_id)){
    if(grid$grid_id[i] %in% rasterIds){
      x <- which(cells$grid_id == grid$grid_id[i]) 
      grid$effort[i] <- cells$for.prop[x]
    } else {
      grid$effort[i] <- mu
    }
  }
  
  rr <- fasterize(grid, r, field="effort")
  names(rr) <- paste0("seal ",IDs[y])
  plot(rr)
  
  outfile <- writeRaster(rr, filename=here::here("Output","Memory_grids","April memory",paste0( paste0("seal ",IDs[y]),".tiff")),
                         format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"), bylayer=TRUE)
}

