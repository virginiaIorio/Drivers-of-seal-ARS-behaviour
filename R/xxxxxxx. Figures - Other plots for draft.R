#Extra things for draft
gpsx <- read.csv("C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/pv64_CSV/pv64-2017_gps_data_with_haulout_&_trip_info.csv", header=TRUE)
gpsx <- gpsx[,c(1:5,11,12)]
gpsx$time <- as.POSIXct(gpsx$D_DATE, tz="UTC")
gpsx$gps_interval <- NA
for(i in 2:length(gpsx$ID)){
  if(gpsx$time[i]>gpsx$time[i-1]){
    gpsx$gps_interval[i] <- difftime(gpsx$time[i], gpsx$time[i-1], unit="mins")
  }
}
mean(gpsx$gps_interval, na.rm=TRUE)
sd(gpsx$gps_interval, na.rm=TRUE)

#Extra track plot
library(ggplot2)
library(ggspatial)
library(ggsn)
library(ggpubr)
library(sf)
library(sp)

coastline <- st_read("C:/Users/r02vi18/PhD_Virginia/Seal behaviour/HMM/R/Bathymetry & Sediment/Coastline_UTM30.shp")


#Start with the plot with all the seals
llcord <- SpatialPoints(gpsx[,c(5,4)],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
gpsx$x <- attr(utmcoord, "coords")[,1]
gpsx$y <- attr(utmcoord, "coords")[,2]

gpsx_all <- gpsx[!gpsx$ID %in% IDs,]
gpsx_seals <- gpsx[which(gpsx$ID %in% IDs),]


map_all <-  ggplot()+ 
  geom_path(data=gpsx_all, aes(x=x, y=y, group=ID), col="grey", lwd=0.5)+
  annotation_spatial(coastline, fill = "lightgrey", lwd = 0)+
  geom_path(data=gpsx_seals, aes(x=x, y=y, group=ID),col="#0072B2", lwd=0.5)+
  xlab("Longitude")+ylab("Latitude")+
  theme_classic()+
  #north(land, location="topright")
  scalebar(y.min =  6375392, y.max= 6581648,
           x.min = 416466.9 , x.max=530568.8,
           dist=20, dist_unit = "km",
           transform = FALSE, model = "WGS84", location = "bottomright",
           st.color="#757575", box.fill=c("#757575","white"), height = 0.009,
           box.color="#757575", st.size=3, border.size=0.1)+
  theme(legend.position="none",
        panel.grid.major = element_line(colour="transparent"),
        text = element_text(size=10))

gpsx$trip_code <- format(gpsx$trip_code, nsmall=3)
trips_analysis <- unique(seal.HMM.flag$ID)
gpsx_seals$trip_in <- ifelse(gpsx_seals$trip_code %in% trips_analysis, 1, 0)

seal_maps <- ggplot()+ 
  geom_path(data=gpsx_seals[which(gpsx_seals$trip_in==0),], aes(x=x, y=y, group=ID), col="gray", lwd=0.5)+
  geom_path(data=gpsx_seals[which(gpsx_seals$trip_in==1),], aes(x=x, y=y, group=ID), col="red", lwd=0.5)+
  annotation_spatial(land, fill = "lightgrey", lwd = 0)+
  facet_wrap(~ID, nrow=2)+
  xlab("Longitude")+
  ylab("Latitude")+
  scalebar(y.min =  6374392, y.max= 6435955,
           x.min = 416466.9 , x.max=470884.8,
           dist=10, dist_unit = "km",
           transform = FALSE, model = "WGS84", location = "bottomright",
           st.color="#757575", box.fill=c("#757575","white"), box.color="#757575",
           st.size=2,height = 0.02, st.dist=0.03, border.size=0.1)+
  theme_bw()+
  theme(legend.position="none",
        panel.grid.major = element_line(colour="transparent"),
        #text = element_text(size=7.5),
        strip.text = element_text(size=10),
        axis.text= element_text(size=6),
        axis.title = element_text(size=10))

arrange <- ggarrange(map_all, seal_maps, ncol=2, widths = c(0.3,0.7), labels=c("A", "B"))
ggsave(plot=arrange, filename="C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/Seal behaviour/HMM/Memory Hp2/Graphs/trips_maps5.tiff",
       device="tiff", width = 10, height=5, units="in", dpi=300)

#####
tripsx <- read.csv("C:/Users/r02vi18/OneDrive - University of Aberdeen/PhD_Virginia/pv64_CSV/pv64-2017_trip_summaries.csv", header=TRUE)
tripsx <- tripsx[which(tripsx$Trip_Duration>=12),]
tripsx_sum <- tripsx %>% group_by(PTT) %>% summarise(n=n(), min_duration = min(Trip_Duration), max_duration = max(Trip_Duration))
tripsx_sum <-tripsx_sum[1:31,]
mean(tripsx_sum$n)

long_tripsx <- tripsx$Trip_Code[which(tripsx$Trip_Duration>12)]
divex <- dive[which(dive$Trip_Code %in% long_tripsx),]
dive_cycle <- divex$DIVE_DUR + divex$SURF_DUR
mean(dive_cycle)

#######
library(ggpubr)


ggscatter(data, x="long.memory", y="batch.foraging.index",
          add="reg.line", conf.int = TRUE,
          cor.coef=TRUE, cor.method = "pearson",
          facet.by = "seal_ID",
          xlab="Memory value", ylab="Foraging index")


tmp <- sealdf[which(!is.na(sealdf$foraging.index)),]
llcord <- SpatialPoints(tmp[,c(9,8)],
                        proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcord, CRS("+proj=utm +zone=30 ellps=WGS84"))
tmp$x <- attr(utmcoord, "coords")[,1]
tmp$y <- attr(utmcoord, "coords")[,2]

seal.id <- 158
tmp1 <- tmp[which(tmp$seal_ID==seal.id),]
tmp1.coords <- tmp1[,c(31,32)]
tmp1$memory <- extract(Longbrick$X158, tmp1.coords)

df <- rbind(df,tmp1)

ggscatter(df, x="memory", y="foraging.index",
          add="reg.line", conf.int = TRUE,
          cor.coef=TRUE, cor.method = "pearson",
          facet.by = "seal_ID",
          xlab="Memory value", ylab="Foraging index")


########
ggplot(seal.batch, aes(x=batch.foraging.index, y=batch.first.foraging))+
  geom_point()+
  xlab("Mean foraging attempts during a dive batch")+
  ylab("Foraging attempts in the first dive of a batch")+
  theme_classic()

ggplot(seal.batch, aes(x=batch.first.foraging, y=batch.last.foraging))+
  geom_point()+
  xlab("Foraging attempts in the first\n dive of a batch")+
  ylab("Foraging attempts in the last\n dive of a batch")+
  theme_classic()
