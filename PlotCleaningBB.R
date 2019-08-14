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



# plotsdatalocation <- "/Users/cowl0037/Documents/Hyper/Final PLOTS/" #mac
# plotsdatalocation <- "/Volumes/HyperDrive/Final PLOTS/" #mac
plotsdatalocation <- "F:/BB_30July/" #isbell pc
substrRight <- function(x, n){
  substr(x, nchar(x)-4-n+1, nchar(x)-4)
}
plot.files <- list.files(plotsdatalocation,pattern=".*shp")
plot.files2 <- plot.files[grepl("Plot",plot.files)]
plot.filesdf<-as.data.frame(plot.files2)
plot.filesdf$Plot <- substrRight(as.character(plot.filesdf$plot.files),3)
plot.filesdf$Plot <- as.numeric(gsub("[^0-9]", "", plot.filesdf$Plot)) 
bb <- read.csv("~/Ortho_proc/e120mega96-18.csv")
bb18 <- bb[bb$Year==2018&bb$Month==8,]
bb_df <- merge(bb18,plot.filesdf,by="Plot")



#read in plot files and compare 0 to 16 
readinfunc <- function(x){
  POI <- st_read(paste0(plotsdatalocation,x))
  POIdf <- as.data.frame(POI)
  POIdf$Plot <- as.numeric(gsub("[^0-9]", "", substrRight(as.character(x),3))) 
  return(POIdf)
}

    Sys.time()
    cl<-makeCluster(no_cores)
    clusterExport(cl,c("rbindlist","st_read","substrRight"))
    plots.filesC <- rbindlist(lapply(unique(plot.filesdf$plot.files)[201:300],readinfunc))
    stopCluster(cl)
    fwrite(plots.filesC,"~/Ortho_Proc/PlotFiles/plots_files201to300.csv")
    Sys.time()
    
    Sys.time()
    cl<-makeCluster(no_cores)
    clusterExport(cl,c("rbindlist","st_read","substrRight"))
    plots.filesD <- rbindlist(lapply(unique(plot.filesdf$plot.files)[301:400],readinfunc))
    stopCluster(cl)
    fwrite(plots.filesD,"~/Ortho_Proc/PlotFiles/plots_files301to400.csv")
    Sys.time()

    Sys.time()
    cl<-makeCluster(no_cores)
    clusterExport(cl,c("rbindlist","st_read","substrRight"))
    plots.filesE <- rbindlist(lapply(unique(plot.filesdf$plot.files)[401:556],readinfunc))
    stopCluster(cl)
    fwrite(plots.filesE,"~/Ortho_Proc/PlotFiles/plots_files401to556.csv")
    Sys.time()
    
    
# system.time(plots.filesA <- rbindlist(lapply(unique(bc_df$plot.files)[1:100],readinfunc)))
plots.filesA <- fread("~/Ortho_Proc/PlotFiles/plots_files1to100.csv")
plots.filesB <- fread("~/Ortho_Proc/PlotFiles/plots_files101to200.csv")
plots.filesC <- fread("~/Ortho_Proc/PlotFiles/plots_files201to300.csv")
plots.filesD <- fread("~/Ortho_Proc/PlotFiles/plots_files301to400.csv")
plots.filesE <- fread("~/Ortho_Proc/PlotFiles/plots_files401to556.csv")


allplots <- rbind(plots.filesA,plots.filesB,plots.filesC,plots.filesD,plots.filesE)


df <- merge(allplots,bb18,by="Plot")

rm(allplots)




