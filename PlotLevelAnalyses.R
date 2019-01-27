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

plot.files <- list.files(ProcLoc,pattern=".*BUFF.*Plot.*shp")
# plot.files <- list.files("~/Downloads/plotshapestmp",pattern=".*BUFF.*Plot.*shp")

system.time(POI <- readOGR(paste0(ProcLoc,plot.files[2])))

# system.time(POI <- readOGR(paste0("~/Downloads/plotshapestmp/",plot.files[2])))
plot(POI)

POI_DF <- as.data.frame(POI)
colnames(POI_DF)[1:272]<-gsub("_",".",colnames(POI_DF)[1:272])
colnames(POI_DF)[1:272]<-substr(colnames(POI_DF)[1:272],3,9)

setDT(POI_DF)
names(POI_DF)
POI_DFsub <- POI_DF[,c(1:272,276,277)]
POI_DFsub$Unique <- 1:nrow(POI_DFsub)
POI_DFsub$NDVI <- (POI_DFsub$`800.614`-POI_DFsub$`660.688`)/(POI_DFsub$`800.614`+POI_DFsub$`660.688`)
POI_DFL<- melt(POI_DFsub,id=c("Lat2","Lon2","Unique"))
POI_DFL$nm <- as.numeric(substr(POI_DFL$variable,3,9))

plot(value~nm,POI_DFL)

ggplot(POI_DFL,aes(nm,value,group=Unique,color=Unique))+geom_line()

# POI_DFL<- melt(POI_DF,patterns=c("nm"),id=c("Lat2","Lon2"))
str(POI_DFL)

ggplot(POI_DFsub,aes(Lon2,Lat2,color=NDVI))+geom_point()

testCluster <- kmeans(POI_DFsub[, 1:272], 10, nstart = 20)
testCluster$betweenss

testCluster2 <- kmeans(POI_DFsub[, 1:272], 9, nstart = 20)
testCluster2$betweenss

testCluster3 <- kmeans(POI_DFsub[, 1:272], 8, nstart = 20)
testCluster3$betweenss


nb <- NbClust(POI_DFsub[, 1:272], diss=NULL, distance = "euclidean", 
                          min.nc=2, max.nc=10, method = "kmeans", 
                          index = "all", alphaBeale = 0.1)
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))


clust_plot <- cbind(POI_DFsub,testCluster$cluster)
colnames(clust_plot)[ncol(clust_plot)]<-"clust"
str(clust_plot)
names(clust_plot)
ggplot(clust_plot,aes(Lon2,Lat2,color=factor(clust)))+geom_point()


