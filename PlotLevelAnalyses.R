### Plot level ML for species ID, etc
#### Jane M Cowles

library(hsdar)
library(caTools)
library(rgdal)
library(raster)
library(spatstat)
library(rgeos)
library(maptools)
library(sf)
library(ggmap)
library(PBSmapping)
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
library(latticeExtra)


#for parallel processing
library(parallel)

#for making a beep noise
library(beepr)

#for clustering
library(NbClust)
library(cluster)
library(factoextra)

##### I can update this later with Sys.info()[['sysname']] -- the biggest hurdle is the st_write where you need to include the shapefile name as part of the dsn in windows but not in mac.

# #directories for MY MAC
# computer <- 'mac'
# RemoteSenDataLoc <- "/Volumes/HyperDrive/Google Drive/remote sensing data/"
# LocalSource <- "~/Dropbox/UMN Postdoc/Ortho_Proc/" #imu and biocon shapefile (polygons)
# ProcLoc <- "/Volumes/HyperDrive/JaneProc/"
# VisLoc <- "/Volumes/HyperDrive/BioCon....."
#directories for ISBELL PC
computer <- 'pc'
RemoteSenDataLoc <- "F:/RemoteSensing/"
LocalSource <- "~/Ortho_Proc/" #imu and biocon shapefile (polygons)
ProcLoc <- "F:/JaneProc/"
VisLoc <- "F:/BioCON10Aug--VISUAL/"


### biocon vis - read in and crop
vis <- raster(paste0(VisLoc,"ortho_biocon_visual.tif"))


# visproj <- projectRaster(vis,crs="+init=epsg:32615 +proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")






plot.files <- list.files(ProcLoc,pattern=".*BUFF.*Plot.*shp")
# plot.files <- list.files("~/Downloads/plotshapestmp",pattern=".*BUFF.*Plot.*shp")

system.time(POI <- readOGR(paste0(ProcLoc,plot.files[58])))
POI$NDVI <- (POI$nm800_614-POI$nm660_688)/(POI$nm800_614+POI$nm660_688)

POI_DF <- as.data.frame(POI)
colnames(POI_DF)[1:272]<-gsub("_",".",colnames(POI_DF)[1:272])
colnames(POI_DF)[1:272]<-substr(colnames(POI_DF)[1:272],3,9)

setDT(POI_DF)
POI_DFsub <- POI_DF[,c(1:272,276,277)]
POI_DFsub$Unique <- 1:nrow(POI_DFsub)
POI_DFsub$NDVI <- (POI_DFsub$`800.614`-POI_DFsub$`660.688`)/(POI_DFsub$`800.614`+POI_DFsub$`660.688`)
hist(POI_DFsub$NDVI)


rast <- raster()
extent(rast) <- extent(POI) 
ncol(rast) <-25
nrow(rast) <- 25

POI_rast <- rasterize(POI, rast, POI$NDVI, fun=mean) 
plot(POI_rast)

system.time(RELPLOT <- crop(vis,extent(bbox(subset(plotshp,PLOTID==104)))+1))

plot(PLOT)
breakpoints <- c(minValue(POI_rast),minValue(POI_rast)+50,minValue(POI_rast)+100,maxValue(POI_rast))
mycol <- rgb(0, 0, 255, max = 255, alpha = 5, names = "blue50")

plot(RELPLOT)
plot(POI_rast,breaks=breakpoints,col=c("blue","yellow","red"),add=T)


POI_DFL<- melt(POI_DFsub,id=c("Lat2","Lon2","Unique","NDVI"))
POI_DFL$nm <- as.numeric(as.character(POI_DFL$variable))

plot(value~nm,POI_DFL)

ggplot(POI_DFL,aes(nm,value,group=Unique,color=Unique))+geom_line()

# POI_DFL<- melt(POI_DF,patterns=c("nm"),id=c("Lat2","Lon2"))
str(POI_DFL)

ggplot(POI_DFsub,aes(Lon2,Lat2,color=NDVI))+geom_point()

testCluster <- kmeans(POI_DFsub[, 1:272], 2, nstart = 20)
fvis_cluster()

testCluster2 <- kmeans(POI_DFsub[, 1:272], 9, nstart = 20)
testCluster2$betweenss

testCluster3 <- kmeans(POI_DFsub[, 1:272], 8, nstart = 20)
testCluster3$betweenss


nb <- NbClust(POI_DFsub[, 1:272], diss=NULL, distance = "euclidean", 
                          min.nc=2, max.nc=10, method = "kmeans", 
                          index = "all", alphaBeale = 0.1)
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

?cluster

clust_plot <- cbind(POI_DFsub,nb$Best.partition)
colnames(clust_plot)[ncol(clust_plot)]<-"clust"
str(clust_plot)
names(clust_plot)
ggplot(clust_plot,aes(Lon2,Lat2,color=factor(clust)))+geom_point()


