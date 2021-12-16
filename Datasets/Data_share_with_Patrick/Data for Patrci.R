depth <- accel_depth[which(accel_depth$ddate>=as.POSIXct("2017-05-26 20:07:00", tz="UTC") & 
                             accel_depth$ddate<=as.POSIXct("2017-05-28 05:35:00", tz="UTC")),]
depth <- depth[,-4]

acc <- pvac[which(pvac$ddate>=as.POSIXct("2017-05-26 20:07:00", tz="UTC") & 
                             pvac$ddate<=as.POSIXct("2017-05-28 05:35:00", tz="UTC")),]

write.table(depth, "C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Datasets/Data_share_with_Patrick/Depth_data_example_14464_trip120.txt",
          sep="\t")

a <- acc$Xs
b <- (acc$Ys^2)+(acc$Zs^2)
acc$pitch <- atan((a/sqrt(b)))
acc$pitch <- acc$pitch*(180/pi)
acc$pitch_deg <- acc$pitch+90

# Calculate roll angle
acc$roll <- atan((-acc$Ys)/acc$Zs)
acc$roll <- acc$roll*(180/pi)


write.table(acc, "C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Datasets/Data_share_with_Patrick/Accelerometer_data_example_14464_trip120.txt",
            sep="\t")

gps <- read.csv("C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Datasets/pv64-2017_gps_data_with_haulout_&_trip_info.csv",
                header=TRUE)
gps1 <- gps[which(gps$PTT==14464 & gps$trip_code==14464.120),]
gps1 <- gps1[,-c(6,7,16:29)]

write.table(gps1, "C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Datasets/Data_share_with_Patrick/GPS_data_example_14464_trip120.txt",
            sep="\t")

dives <- read.csv("C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Datasets/pv64-2017_dive.csv",
                header=TRUE)
dive1 <- dives[which(dives$PTT==14464 & dives$Trip_No==120),]

write.table(dive1, "C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Datasets/Data_share_with_Patrick/Dive_data_example_14464_trip120.txt",
            sep="\t")


llcord <- SpatialPoints(gps1[,c(5,4)],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
gps1$x <- attr(utmcoord, "coords")[,1]
gps1$y <- attr(utmcoord, "coords")[,2]

ggplot()+
  geom_path(data=gps1, aes(x=x, y=y), lwd=0.5, col="red")+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  xlim(min(gps1$x)-20000, max(gps1$x)+120000)+
  ylim(min(gps1$y)-40000, max(gps1$y)+120000)+
  xlab("Longitude")+ylab("Latitude")+
  theme_classic()+
  theme(legend.position="none",
        panel.grid.major = element_line(colour="transparent"),
        text = element_text(size=12))



