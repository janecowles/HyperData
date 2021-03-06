#### Cleaning of orthorectified data.

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

# plotsdatalocation <- "/Users/cowl0037/Documents/Hyper/Final PLOTS/" #mac
plotsdatalocation <- "/Volumes/HyperDrive/Final PLOTS/" #mac
# plotsdatalocation <- "F:/Final PLOTS/" #isbell pc
substrRight <- function(x, n){
  substr(x, nchar(x)-4-n+1, nchar(x)-4)
}
plot.files <- list.files(plotsdatalocation,pattern=".*shp")
plot.filesdf<-as.data.frame(plot.files)
plot.filesdf$Plot <- substrRight(as.character(plot.filesdf$plot.files),3)
plot.filesdf$Plot <- as.numeric(gsub("[^0-9]", "", plot.filesdf$Plot)) 
bc <- read.csv("/Users/cowl0037/Dropbox/UMN Postdoc/CDR DATA/BioCON Master Harvest_190109_USE.csv")
bc18 <- bc[bc$year==2018&bc$Season=="August",]
bc_df <- merge(bc18,plot.filesdf,by="Plot")

bc_df0 <- bc_df[bc_df$CountOfSpecies==0,]
bc_df16 <- bc_df[bc_df$CountOfSpecies==16,]
bc_df16sub <- bc_df16[bc_df16$plot.files%in%sample(bc_df16$plot.files,8),]
bc_df016<-bc_df[bc_df$CountOfSpecies%in%c(0,16),]

#read in plot files and compare 0 to 16 
readinfunc <- function(x){
  POI <- st_read(paste0(plotsdatalocation,x))
  POIdf <- as.data.frame(POI)
  POIdf$Plot <- as.numeric(gsub("[^0-9]", "", substrRight(as.character(x),3))) 
  return(POIdf)
}

# POI <- st_read(paste0(plotsdatalocation,plot.files[100]))

# system.time(zeroplots <- rbindlist(lapply(unique(bc_df0$plot.files),readinfunc)))

# system.time(sixteenplots <- rbindlist(lapply(unique(bc_df16sub$plot.files),readinfunc)))
# system.time(zero16plots <- rbindlist(lapply(unique(bc_df016$plot.files),readinfunc)))

# system.time(plots.files250to447 <- rbindlist(lapply(unique(bc_df$plot.files)[250:447],readinfunc)))
# system.time(plots.files500to554 <- rbindlist(lapply(unique(bc_df$plot.files)[500:554],readinfunc)))
plots.files250to447 <- fread("/Users/cowl0037/Downloads/plots_files250to447.csv")
plots.files447to555 <- fread("/Users/cowl0037/Downloads/plots_files447to555.csv")
#accidentally made an overlap!
plots.filesno448 <-plots.files447to555[!(plots.files447to555$Plot==297&plots.files447to555$StartFrame==39773),]
# unique(plots.files447to555$StartFrame[plots.files447to555$Plot==297])
# unique(plots.files250to447$StartFrame[plots.files250to447$Plot==297])

plots.files1to249 <- fread("/Users/cowl0037/Downloads/plots_files1to249.csv")


allplots <- rbind(plots.files1to249,plots.files250to447,plots.filesno448)



rm(plots.files1to249)
rm(plots.files250to447)
rm(plots.files447to555)
rm(plots.filesno448)

df <- merge(allplots,bc18,by="Plot")

rm(allplots)

# df$mean430660 <- rowMeans(df[,17:120])
# hist(df$mean430660)
df.1 <- df[df$Ring==1,] #rm(df)
df.2 <- df[df$Ring==2,]
df.3 <- df[df$Ring==3,]
df.4 <- df[df$Ring==4,]
df.5 <- df[df$Ring==5,]
df.6 <- df[df$Ring==6,]



ggplot(df.1[df$mean430660<800&df$Ring==1,],aes(Lon2,Lat2,color=mean430660))+geom_point(size=2.5)+scale_color_continuous(low="red",high="blue")+theme_classic()
ggplot(df[df$mean430660<800&df$Ring==1,],aes(Lon2,Lat2,color=factor(CountOfSpecies)))+geom_point(size=2.5)+theme_classic()

table(df$Plot)


zero16 <- rbind(zeroplots,sixteenplots)
zero16$mean430660 <- rowMeans(zero16[,16:119])
zero16$mean600660 <- rowMeans(zero16[,92:119])
zeroplots$mean430660 <- rowMeans(zeroplots[,16:119])
ggplot(zeroplots,aes(Lon2,Lat2,color=nm540_751))+geom_point(size=2.5)+scale_color_continuous(low="red",high="blue")+theme_classic()+facet_wrap("Plot",scales="free")
ggplot(zeroplots,aes(Lon2,Lat2,color=mean430660))+geom_point(size=2.5)+scale_color_continuous(low="red",high="blue")+theme_classic()+facet_wrap("Plot",scales="free")



ggplot(sixteenplots,aes(Lon2,Lat2,color=nm540_751))+geom_point(size=2.5)+scale_color_continuous(low="red",high="blue")+theme_classic()+facet_wrap("Plot",scales="free")
ggplot(zero16,aes(Lon2,Lat2,color=mean430660))+geom_point(size=2.5)+scale_color_continuous(low="red",high="blue")+theme_classic()+facet_wrap("Plot",scales="free")
ggplot(zero16,aes(Lon2,Lat2,color=mean600660))+geom_point(size=2.5)+scale_color_continuous(low="red",high="blue")+theme_classic()+facet_wrap("Plot",scales="free")


system.time(POI <- readOGR(paste0(plotsdatalocation,plot.files[197])))
plot(POI,col=POI$nm398_604)
POIdf <- as.data.frame(POI)
ggplot(POIdf,aes(Lon2,Lat2,color=nm540_751))+geom_point(size=2.5)+scale_color_continuous(low="forestgreen",high="tan")

POIdf$mean440660 <- rowMeans(POIdf[,20:119])

# Read in all plots, determine classifications


