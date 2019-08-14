#white ref
# load FinalFULL1024.shp
# load BrightestWhiteLoc.shp
# cut brightest white from final full - take averages.
# end product, 1 row with values for each wavelength
library(sf)
library(data.table)
library(rgdal)
# finalfull1024 <- st_read("/Volumes/HyperDrive/Final FULL/Final1024full.shp")
sppoi <- as(finalfull1024,"Spatial")

# cut out white ref by hand in R?
library(raster)
# poi <- st_read("/Volumes/HyperDrive/Final PLOTS/Final1024Plot62.shp")
# sppoi <- as(poi,"Spatial")
# 
# finalfull1024 <- st_read("/Volumes/HyperDrive/Final FULL/Final1024full.shp")
# sppoi <- as(finalfull1024,"Spatial")
# 
# gr <- sppoi$nm540_751
# clye <- rgb(0, 0, 255, max = 255, alpha = 50, names = "yellow")
# clbl <- rgb(0, 0, 255, max = 255, alpha = 5, names = "blue")
# clgr <- rgb(0, 0, 255, max = 255, alpha = 5, names = "green")
# color <- ifelse(gr > 0 & gr <= 300, "blue",ifelse(gr > 300 & gr <= 400, "green",ifelse(gr > 400, "yellow", NA)))
# plot(sppoi,col=color)
# refareagenerally <- click(sppoi)
# sprefxy <- refareagenerally[,c("Lon2","Lat2")]
# spref <- SpatialPointsDataFrame(coords=sprefxy,data=refareagenerally,proj4string = CRS("+init=epsg:32615")) 
# color <- ifelse(gr > 0 & gr <= 400, clbl,ifelse(gr > 400 & gr <= 800, "green",ifelse(gr > 800 & gr <= 1000, "yellow", ifelse(gr > 1000, "red", NA))))
# plot(spref,col=color)


bright <- st_read("~/Ortho_Proc/BrightWhiteRef_BBFab_30July.shp")
brightdf <- as.data.frame(bright)
setDT(brightdf)
bright_means <- brightdf[,lapply(.SD, mean),  .SDcols = 1:272]
bright_means$Type <- "BrightRef"
fwrite(bright_means,"~/Ortho_Proc/BrightRef_meanvalBB.csv",row.names=F)
dark <- read.ENVI("F:/RemoteSensing/BigBio-eDNA-2019/BigBio-FAB1-2019-07-30-Hyper/100061_dark_2019_07_30_16_04_17/raw_0")
init.dim <- dim(dark)
dim(dark) <- c(dim(dark)[1]*dim(dark)[2],dim(dark)[3])

darkdf <- as.data.frame(dark)
rm(dark)
bandtowave <- read.csv("~/Ortho_Proc/BandNumWavelength.csv")
colnames(darkdf)[1:272]<-gsub("\\.","_",paste0("nm",bandtowave$Wavelength))
setDT(darkdf)
dark_means <- darkdf[,lapply(.SD, mean),  .SDcols = 1:272]
dark_means$Type <- "DarkRef"
fwrite(dark_means,"~/Ortho_Proc/DarkRef_meanvalBB.csv")
dark_means <- read.csv("~/Ortho_Proc/DarkRef_meanval.csv")
### for each wavelength -- (target-dark)/(light-dark)

light_dark_df1 <- rbind.fill(allplots,bright_means,dark_means)
names(light_dark_df1)
i=3
system.time(for(i in 1:272){
  # print((light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="BrightRef",i]-light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="DarkRef",i]))
  tmp <- (light_dark_df1[,i]-light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="DarkRef",i])/(light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="BrightRef",i]-light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="DarkRef",i])
  light_dark_df1<- cbind(light_dark_df1,tmp)
  colnames(light_dark_df1)[ncol(light_dark_df1)]<- gsub("_","\\.",colnames(light_dark_df1)[i])
  
})

# corrDF <- light_dark_df1[,c(1,275:294,310:315,317:334,338:362,383:404,470:742)]
corrDF <- light_dark_df1[,c(273:554)]


fwrite(corrDF,"~/REFCORRECTED_DATABIGBIO.csv")

