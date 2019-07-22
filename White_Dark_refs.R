#white ref
# load FinalFULL1024.shp
# load BrightestWhiteLoc.shp
# cut brightest white from final full - take averages.
# end product, 1 row with values for each wavelength
library(sf)
library(data.table)
library(rgdal)
# finalfull1024 <- st_read("/Volumes/HyperDrive/Final FULL/Final1024full.shp")
bright <- st_read("~/Dropbox/UMN Postdoc/Ortho_Proc/BrightestWhiteRef_JMCCut_31May2019_fromFinalFULL1024.shp")
brightdf <- as.data.frame(bright)
setDT(brightdf)
bright_means <- brightdf[,lapply(.SD, mean),  .SDcols = 1:272]
bright_means$Type <- "BrightRef"
fwrite(bright_means,"~/Dropbox/UMN Postdoc/Ortho_Proc/BrightRef_meanval.csv",row.names=F)
dark <- readGDAL("/Volumes/HyperDrive/Google Drive/remote sensing data/20180917/100039_darkBioCON_2018_09_17_14_46_37/raw_0")
darkdf <- as.data.frame(dark)
rm(dark)
bandtowave <- read.csv("~/Dropbox/UMN Postdoc/Ortho_Proc/BandNumWavelength.csv")
colnames(darkdf)[1:272]<-gsub("\\.","_",paste0("nm",bandtowave$Wavelength))
setDT(darkdf)
dark_means <- darkdf[,lapply(.SD, mean),  .SDcols = 1:272]
dark_means$Type <- "DarkRef"
fwrite(dark_means,"~/Dropbox/UMN Postdoc/Ortho_Proc/DarkRef_meanval.csv")
dark_means <- read.csv("~/Dropbox/UMN Postdoc/Ortho_Proc/DarkRef_meanval.csv")
### for each wavelength -- (target-dark)/(light-dark)

light_dark_df1 <- rbind.fill(df,bright_means,dark_means)

i=3
system.time(for(i in 3:274){
  # print((light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="BrightRef",i]-light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="DarkRef",i]))
  tmp <- (light_dark_df1[,i]-light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="DarkRef",i])/(light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="BrightRef",i]-light_dark_df1[!is.na(light_dark_df1$Type)&light_dark_df1$Type=="DarkRef",i])
  light_dark_df1<- cbind(light_dark_df1,tmp)
  colnames(light_dark_df1)[ncol(light_dark_df1)]<- gsub("_","\\.",colnames(light_dark_df1)[i])
  
})

corrDF <- light_dark_df1[,c(1,275:294,310:315,317:334,338:362,383:404,470:742)]


fwrite(corrDF,"~/Documents/REFCORRECTED_DATA.csv")

