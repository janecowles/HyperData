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

plot.filesdf <- data.frame(plot.files)

substrRight <- function(x, n){
  substr(x, nchar(x)-4-n+1, nchar(x)-4)
}
plot.filesdf$Plots <- substrRight(as.character(plot.filesdf$plot.files),3)
plot.filesdf$Plots <- as.numeric(gsub("[^0-9]", "", plot.filesdf$Plots)) 
plot.filesdf[plot.filesdf$Plots%in%c(46,77,142,146,229,237,271,327),]
plot.filesdf[plot.filesdf$Plots%in%c(25,61,73,182,217,284,289,367,368,369,370,371),]
plot.filesdf[plot.filesdf$Plots%in%c(256),]

plot.filesdf$Frames1 <- substrRight(as.character(plot.filesdf$plot.files),13)
plot.filesdf$Frames <- substr(as.character(plot.filesdf$Frames1),1,8)
plot.filesdf$Frames <- as.numeric(gsub("[^0-9]", "", plot.filesdf$Frames))
unique(plot.filesdf$Frames)

system.time(POI <- readOGR(paste0(ProcLoc,plot.filesdf[197,1])))
POI$NDVI <- (POI$nm800_614-POI$nm660_688)/(POI$nm800_614+POI$nm660_688)

POI_DF <- as.data.frame(POI)
colnames(POI_DF)[1:272]<-gsub("_",".",colnames(POI_DF)[1:272])
colnames(POI_DF)[1:272]<-substr(colnames(POI_DF)[1:272],3,9)

setDT(POI_DF)
POI_DFsub <- POI_DF[,c(1:272,276,277)]
POI_DFsub$Unique <- 1:nrow(POI_DFsub)
POI_DFsub$NDVI <- (POI_DFsub$`800.614`-POI_DFsub$`660.688`)/(POI_DFsub$`800.614`+POI_DFsub$`660.688`)
hist(POI_DFsub$NDVI)
POI_long <- melt(POI_DFsub, id.vars=c("Lat2","Lon2","Unique","NDVI"),measure.vars=c(1:272))
POI_long$nm <- as.numeric(as.character(POI_long$variable))
ggplot(POI_long,aes(nm,value,group=Unique,color=Unique))+geom_line()+labs(title='POI')


rast <- raster()
extent(rast) <- extent(POI) 
ncol(rast) <-25
nrow(rast) <- 25

POI_rast <- rasterize(POI, rast, POI$NDVI, fun=mean) 
plot(POI_rast)


system.time(POI2 <- readOGR(paste0(ProcLoc,plot.filesdf[293,1])))
POI2$NDVI <- (POI2$nm800_614-POI2$nm660_688)/(POI2$nm800_614+POI2$nm660_688)

POI2_DF <- as.data.frame(POI2)
colnames(POI2_DF)[1:272]<-gsub("_",".",colnames(POI2_DF)[1:272])
colnames(POI2_DF)[1:272]<-substr(colnames(POI2_DF)[1:272],3,9)

setDT(POI2_DF)
POI2_DFsub <- POI2_DF[,c(1:272,276,277)]
POI2_DFsub$Unique <- 1:nrow(POI2_DFsub)
POI2_DFsub$NDVI <- (POI2_DFsub$`800.614`-POI2_DFsub$`660.688`)/(POI2_DFsub$`800.614`+POI2_DFsub$`660.688`)
hist(POI2_DFsub$NDVI)
POI2_long <- melt(POI2_DFsub, id.vars=c("Lat2","Lon2","Unique","NDVI"),measure.vars=c(1:272))
POI2_long$nm <- as.numeric(as.character(POI2_long$variable))

ggplot(POI2_long,aes(nm,value,group=Unique,color=Unique))+geom_line()+labs(title='POI2')


rast <- raster()
extent(rast) <- extent(POI2) 
ncol(rast) <-25
nrow(rast) <- 25

POI2_rast <- rasterize(POI2, rast, POI2$NDVI, fun=mean) 
par(mfrow=c(1,2))
plot(POI_rast,main="4-14CorrProcBUFF39773Plot256")
plot(POI2_rast,main="4-15CorrProcBUFF37901Plot256")




system.time(RELPLOT <- crop(vis,extent(bbox(subset(plotshp,PLOTID==256)))+1))

plot(RELPLOT)
breakpoints <- c(minValue(POI_rast),minValue(POI_rast)+0.130,minValue(POI_rast)+0.150,maxValue(POI_rast))
mycol <- rgb(0, 0, 255, max = 255, alpha = 5, names = "blue50")

plot(RELPLOT)
plot(POI_rast,breaks=breakpoints,col=c(mycol,"blue","yellow","red"),add=T)
plot(RELPLOT)
plot(POI2_rast,breaks=breakpoints,col=c(mycol,"blue","yellow","red"),add=T)


POI_DFL<- melt(POI_DFsub,id=c("Lat2","Lon2","Unique","NDVI"))
POI_DFL$nm <- as.numeric(as.character(POI_DFL$variable))

plot(value~nm,POI_DFL)

ggplot(POI_DFL,aes(nm,value,group=Unique,color=Unique))+geom_line()

# POI_DFL<- melt(POI_DF,patterns=c("nm"),id=c("Lat2","Lon2"))
str(POI_DFL)

ggplot(POI2_DFsub,aes(Lon2,Lat2,color=NDVI))+geom_point()

testCluster <- kmeans(POI2_DFsub[, 1:272], 2, nstart = 20)
fviz_cluster(testCluster,POI2_DFsub[,1:272])

testCluster <- kmeans(POI_DFsub[, 1:272], 2, nstart = 20)
fviz_cluster(testCluster,POI_DFsub[,1:272])


testCluster2 <- kmeans(POI_DFsub[, 1:272], 9, nstart = 20)
testCluster2$betweenss

testCluster3 <- kmeans(POI_DFsub[, 1:272], 2, nstart = 20)
testCluster3$betweenss


nb <- NbClust(POI_DFsub[, 20:119], diss=NULL, distance = "euclidean", 
                          min.nc=2, max.nc=10, method = "kmeans", 
                          index = "all", alphaBeale = 0.1)
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

?cluster

clust_plot <- cbind(POI_DFsub,nb$Best.partition)
colnames(clust_plot)[ncol(clust_plot)]<-"clust"
ggplot(clust_plot,aes(Lon2,Lat2,color=factor(clust)))+geom_point(size=4)

nb2 <- NbClust(POI2_DFsub[, 20:119], diss=NULL, distance = "euclidean", 
              min.nc=2, max.nc=4, method = "kmeans", 
              index = "all", alphaBeale = 0.1)
hist(nb2$Best.nc[1,], breaks = max(na.omit(nb2$Best.nc[1,])))

?cluster

clust_plot2 <- cbind(POI2_DFsub,nb2$Best.partition)
colnames(clust_plot2)[ncol(clust_plot2)]<-"clust"
ggplot(clust_plot2,aes(Lon2,Lat2,color=factor(clust)))+geom_point(size=4)


par(mfrow=c(1,2))
plot(POI_rast,main="4-14CorrProcBUFF39773Plot256")
plot(POI2_rast,main="4-15CorrProcBUFF37901Plot256")
ggplot(clust_plot,aes(Lon2,Lat2,color=factor(clust)))+geom_point(size=4)+labs(title="4-14CorrProcBUFF39773Plot256")
ggplot(clust_plot2,aes(Lon2,Lat2,color=factor(clust)))+geom_point(size=4)+labs(title="4-15CorrProcBUFF37901Plot256")

