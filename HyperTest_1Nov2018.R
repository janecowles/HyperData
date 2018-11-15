### jmc hyperspectral processing
### 9 November 2018 BIOCON

setwd("/Volumes/HyperDrive/Google Drive/remote sensing data/20180917/100040_bc_2018_09_17_14_48_50")
list.files()

# install.packages("sf")

### for hyperspectral processing
library(hsdar)
library(caTools)
library(rgdal)
library(raster)
library(spatstat)
library(rgeos)
library(maptools)
library(sf)
### general packages that I use often
library(plyr)
library(nlme)
library(lme4)
library(car)
library(ggplot2)
library(vegan)
library(readxl)
library(data.table)
library(tidyr)
library(gridExtra)
library(RColorBrewer)
# tstr <- raster("raw_1024")

# str(tstr)
# plot(tstr)

imu <- read.delim("imu_gps.txt")
str(imu)
setDT(imu,key="Timestamp")


imuxy <- imu[,c("Lon","Lat")]
imusp <- SpatialPointsDataFrame(coords=imuxy,data=imu,proj4string = CRS("+init=epsg:4326"))
plot(imusp,col="green")


# plot(imu$Lat~imu$Lon,col=imu$Alt*4)
# ggplot(imu,aes(Lon,Lat,color=Alt))+geom_point()+scale_color_continuous(low="red",high="green")


# plot(imu$Timestamp,rep(1,nrow(imu)))
# points(fr.ind$Timestamp,rep(1.01,nrow(fr.ind)),col="red")

fr.ind <- NULL
frames <- list.files(pattern="frameIndex")

for(i in frames){
  tmp <- read.delim(i)
  fr.ind <- rbind(fr.ind,tmp)
}
setDT(fr.ind,key="Timestamp")

fr.indximu <- imu[fr.ind,roll=T]
fr.indximu$Frame.m <- fr.indximu$Frame.+0.5

hist(fr.indximu$Frame.m)

setDT(fr.indximu,key = "Frame.m")

### ok if i can take the frame indices and match to the spatial data frame, a row would be a frame, so we could merge the information from frame 1024 with the y=1 row, 1025 would be y=2, and y=1791 would be frame 2815

# tstr <- raster("raw_1024")

tstreadgdal <-readGDAL("raw_1024")

# rel.fr <- fr.indmatch[fr.indmatch$Frame.m>1024&fr.indmatch$Frame.m<2816,]
tst <- tstreadgdal
str(tst)
plot(tst)

str(tst)

coordinates(tst)
coordinates(tst)[,2]
tst@data$yco <-coordinates(tst)[,2]


###################### 15 November 2018



####### NEED TO FIGURE OUT X COORDINATES...

tstdf <- as.data.frame(tst)
tstdf$Frame.m <- tstdf$y+1024
setDT(tstdf,key = "Frame.m")
names(tstdf)

tstshort <- tstdf[tstdf$Frame.m>1900&tstdf$Frame.m<2100,]

comb <- fr.indximu[tstshort]

summary(comb$Frame.m)

xy <- comb[,c("Lon","Lat")]
combsp <- SpatialPointsDataFrame(coords=xy,data=comb,proj4string = CRS("+init=epsg:4326"))
plot(combsp,add=T,col="blue")


tst.m <- merge(tstshort,fr.indmatch,all.x=T,all.y=F,by.x="y",by.y="Frame.m")
tst.c <- tst.m[!is.na(tst.m$band1),]
###as.data.frame -> merge -> as rasterbrick with new coordinates



#### THIS MAY WORK? Wait, unsure if this will work since data may be a list instead of a matrix. May need to update coordinates by naming coordinates with row=1 as y=1024 and figuring out the x based on the x of the drone (middle) and the extent of the image (the edges)

tstsl <- speclib(tst,wavelength = c(1:272))
plot(tstsl)


tst <- read.ENVI("raw_1024")
head(tst)
tst2 <-speclib(tst)

###### biocon shapefiles
bcshp <- readOGR("~/Dropbox/UMN Postdoc/E141 Shapes/e141_plot_correction(Updating).shp")
bcshpt <- spTransform(bcshp, CRS("+init=epsg:4326"))

str(bcshp)

bcshpl <- readOGR("~/Dropbox/UMN Postdoc/E141 Shapes/e141.shp")
bcshptl <- spTransform(bcshpl, CRS("+init=epsg:4326"))
bcsf<-st_as_sf(bcshptl)
bcsfp <- st_polygonize(bcsf)

plot(bcsfp)

str(bcshp,list.len=20)
bbox(bcshp)
bbox(imusp)
plot(imusp)
plot(bcshpt,add=T,col="magenta")
plot(combsp,add=T,col="green")


plot(combsp)
plot(bcshpt2,add=T,col="magenta")


##### biocon visual shapefile
vissh <- readOGR("/Volumes/HyperDrive/Google Drive/remote sensing data/BioCON visual/BioCON visual 10 August 2018_volume_shapefile/BioCON visual 10 August 2018_volume.shp")
plot(vissh)
contsh <- readOGR("/Volumes/HyperDrive/Google Drive/remote sensing data/BioCON visual/contours/contours.shp")
plot(contsh)
vistif <- raster("/Volumes/HyperDrive/Google Drive/remote sensing data/BioCON visual/ortho.tif")
plot(vistif)

#### big bio shapefile

bbsh <- readOGR("~/Documents/BIG BIO DRONES/big bio shape files from CDR/BigBioCombined.shp")
bbsh$Id
