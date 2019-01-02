#### Trig Calcs for Drone
#### Jane M Cowles
#### Dec 22 2018
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

#for parallel processing
library(parallel)

#for making a beep noise
library(beepr)


##### I can update this later with Sys.info()[['sysname']] -- the biggest hurdle is the st_write where you need to include the shapefile name as part of the dsn in windows but not in mac.

#directories for MY MAC
computer <- 'mac'
RemoteSenDataLoc <- "/Volumes/HyperDrive/Google Drive/remote sensing data/"
LocalSource <- "~/Dropbox/UMN Postdoc/Ortho_Proc/" #imu and biocon shapefile (polygons)
ProcLoc <- "/Volumes/HyperDrive/JaneProc/"

#directories for ISBELL PC
computer <- 'pc'
RemoteSenDataLoc <- "F:/RemoteSensing/"
LocalSource <- "~/Ortho_Proc/" #imu and biocon shapefile (polygons)
ProcLoc <- "F:/JaneProc/"

# fakefunction <- function(x){
#   print(x)
#   print(computer)
#   }
# fakefunction("will the computer know what I set as computer outside of the loop???")


no_cores <- detectCores()-1

degtorad <- function(x){x*pi/180} #convert degrees to radians for trig functions in R


#coords.epsg = "32615" is utm 15n which is what I want here. "4326" is lat lon

imu <- read.delim(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/imu_gps.txt"))
overallIMUmin <- min(imu$Alt)
rm(imu)#keep things tidy
imu.framematch <- read.csv(paste0(LocalSource,"FramexIMU.csv"))

imu_proc <- function(imu.datafile,FOVAngle,GroundLevel=0,degree=TRUE,coords.epsg){
  
 tmp <- imu.datafile
 
 tmp$HeightAboveGround <- tmp$Alt-GroundLevel
 
    if(degree==TRUE){
    tmp$FOVAngle <- degtorad(FOVAngle)
    tmp$Roll <- degtorad(tmp$Roll)
    tmp$Pitch <- degtorad(tmp$Pitch)
    tmp$Yaw <- degtorad(tmp$Yaw)
    }
 
 #at each point, this is the FOV
  tmp$UncorrectedFOVmeters <- 2*tmp$HeightAboveGround*tan(tmp$FOVAngle/2)
  
  #ROLL
   tmp$Rolloffset <- tmp$HeightAboveGround*tan(tmp$Roll) #parallel to the frame's long direction (as specified by yaw), therefore:
   tmp$RollYComponent <- tmp$Rolloffset*sin(tmp$Yaw) ###(FOR NOW REMOVING THIS NEGATIVE SIGN TO SEE HOW THINGS LOOK!!!!!)
   tmp$RollXComponent <- tmp$Rolloffset*cos(tmp$Yaw)
   
  #PITCH
  tmp$Pitchoffset <- tmp$HeightAboveGround*tan(tmp$Pitch) # orthogonal to the frame's long direction (as specified by the yaw), therefore:
  tmp$PitchYComponent <-tmp$Pitchoffset*cos(tmp$Yaw)
  tmp$PitchXComponent <-tmp$Pitchoffset*sin(tmp$Yaw)
  
  #TOTAL
  tmp$TotMidPointCorrectionY <- tmp$RollYComponent+tmp$PitchYComponent
  tmp$TotMidPointCorrectionX <- tmp$RollXComponent+tmp$PitchXComponent

  
  #Change in FOV length due to ROLL. Potential to add PITCH to this to make it even more accurate. If pitch is not 0, the distance away will be even farther.
  #The new distance will be the base of two right triangles (apex of each is the FOV angle + or - the Roll angle)
  
  #triangle 1
  Tri1ApexAngle <- 0.5*tmp$FOVAngle + tmp$Roll
  Tri1Base <- tmp$HeightAboveGround*tan(Tri1ApexAngle)
  #triangle 2
  Tri2ApexAngle <- 0.5*tmp$FOVAngle - tmp$Roll
  Tri2Base <- tmp$HeightAboveGround*tan(Tri2ApexAngle)

#if roll=0, these will be identical

#totalbase = FOV on the ground
tmp$RollCorrectedFOVmeters <- Tri1Base + Tri2Base

#Yaw offset -- used in finding the endpoints -- left side will use -x and +y, right side will use +x and -y (with +yaw) -- Or vice versa -- check on this later -- I believe I swapped this around in the later code.
tmp$YawYOffset<-0.5*tmp$RollCorrectedFOVmeters*sin(tmp$Yaw)
tmp$YawXOffset<-0.5*tmp$RollCorrectedFOVmeters*cos(tmp$Yaw)
  
sp.tmp <- SpatialPointsDataFrame(coords=tmp[,c("Lon","Lat")],data=tmp,proj4string = CRS(paste0("+init=epsg:",coords.epsg)))
sp.tmpTR <- spTransform(sp.tmp,"+init=epsg:32615")

sp.out <- sp.tmpTR[,!(names(sp.tmpTR)%in%c("Lon","Lat"))]

out <- data.frame(sp.out)

return(out)


}

#note I have min(imu$Alt) -- only works because I have the full imu loaded (instead of the imu file of just the positions with images taken)
Proc_IMU <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326)
names(Proc_IMU)
rm(imu.framematch)



spacing_fun <- function(pix_i,matchimu,specdfframe){
  specdfpix <- specdfframe[specdfframe$x==pix_i,]
  
  ########### THINK ABOUT THIS
  ########### Do I need to flip the + and - here to get these in the right direction. Should it be +Y and -X?????
  ########### I ask because the images are just perfectly flipped -- gasp!
  
  
  #Original!!!!
  # specdfpix$Lat2<-matchimu$Lat+matchimu$TotMidPointCorrectionY+matchimu$YawYOffset-(pix_i*2*matchimu$YawYOffset/640)
  # specdfpix$Lon2<-matchimu$Lon+matchimu$TotMidPointCorrectionX-matchimu$YawXOffset+(pix_i*2*matchimu$YawXOffset/640)

  #ALTERED -- should be complete reversal of the yaw stuff. Maybe also the midpoint, but I did NOT do that here.
  specdfpix$Lat2<-matchimu$Lat+matchimu$TotMidPointCorrectionY-matchimu$YawYOffset+(pix_i*2*matchimu$YawYOffset/640)
  specdfpix$Lon2<-matchimu$Lon+matchimu$TotMidPointCorrectionX+matchimu$YawXOffset-(pix_i*2*matchimu$YawXOffset/640)
  
  
  ####new attempt -- minus the midpoint correction, too!
  
  return(specdfpix)
}

centerandends_corr <- function(framex,spectral.data.frame,ProcessedIMU){
   specdfframe <- spectral.data.frame[spectral.data.frame$frame==framex,]
    matchimu <- ProcessedIMU[ProcessedIMU$Frame.==framex,]
    specdfframebypix <- rbindlist(lapply(sort(unique(specdfframe$x)),spacing_fun,matchimu,specdfframe))
    print(framex)
    return(specdfframebypix)
}


shpfile_plotloop <-function(plotnum,specdfOUT_sf,PlotShapeFile,filenumber){
  plottmp <- subset(PlotShapeFile,PLOTID==plotnum)
  plot_sf <- st_as_sf(plottmp)
  #for now, I want to test if I can see the plot edges as a square within the buffered region
  # plot_buff <- gBuffer(subset(PlotShapeFile,PLOTID==plotnum),width = 0.5)
  # plot_sf <- st_as_sf(plot_buff)
  st_crs(plot_sf)<-st_crs(specdfOUT_sf)
  plot_clip <- st_intersection(specdfOUT_sf,plot_sf)

  if(dim(plot_clip)[1]>0){
    print(plotnum)
    if(computer=="mac"){
  st_write(plot_clip,dsn=paste0(ProcLoc),layer=paste0("BUFFProc",filenumber,"Plot",plotnum),driver="ESRI Shapefile",delete_layer=TRUE)}

  if(computer=="pc"){
    st_write(plot_clip,dsn=paste0(ProcLoc,"BUFFProc",filenumber,"Plot",plotnum,".shp"),layer=paste0("BUFFProc",filenumber,"Plot",plotnum),driver="ESRI Shapefile",delete_layer=TRUE)}
  }
}

# shpfile_ringloop <- function(ringnum,specdfOUT_sf,PlotShapeFile,filenumber){
#   
#     ring_poly <- as(extent(bbox(subset(PlotShapeFile,RingID==ringnum))),"SpatialPolygons")
#   ring_buff <- gBuffer(ring_poly,width = 2)
#   ring_sf <- st_as_sf(ring_buff)
#   st_crs(ring_sf)<-st_crs(specdfOUT_sf)
#   ring_clip <- st_intersection(specdfOUT_sf,ring_sf)
# 
#   if(dim(ring_clip)[1]>0){
#     print(ringnum)
#     if(computer=="mac"){
#       st_write(ring_clip,dsn=paste0(ProcLoc),layer=paste0("BUFFProc",filenumber,"ring",ringnum),driver="ESRI Shapefile",update=TRUE)}
#     
#     if(computer=="pc"){
#       st_write(ring_clip,dsn=paste0(ProcLoc,"BUFFProc",filenumber,"ring",ringnum,".shp"),layer=paste0("BUFFProc",filenumber,"ring",ringnum),driver="ESRI Shapefile",update=TRUE)}
# }
# }


ortho_fun <- function(filenumber,ProcessedIMU,PlotShapeFile){
  orig_sp <-readGDAL(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/raw_",filenumber))
  spectral.data.frame <- as.data.frame(orig_sp)
  orig_sp <- NULL
  
  #### here is where I want to FLIP IT!!!!
  ###flipped version: (NEXT: Do i want this - 0.5 anymore??)
    spectral.data.frame$frame <- filenumber+max(spectral.data.frame$y)-spectral.data.frame$y
### original version
  # spectral.data.frame$frame <- filenumber+spectral.data.frame$y-0.5
  spectral.data.frame$Lat2 <- NA
  spectral.data.frame$Lon2 <- NA
  cl<-makeCluster(no_cores)
  clusterExport(cl,c("rbindlist","spacing_fun"))
  specdfOUT<- rbindlist(parLapply(cl,sort(unique(spectral.data.frame$frame)),centerandends_corr,spectral.data.frame,ProcessedIMU))
  stopCluster(cl)
  specdfOUT$StartFrame <- filenumber
  # specdfOUT$Lat <- format(specdfOUT$Lat2,digits=20)
  # specdfOUT$Lon <- format(specdfOUT$Lon2,digits=20)

  # write.csv(specdfOUT,paste0("/Volumes/HyperDrive/JaneProc/Proc",filenumber,".csv"))
specdfOUT_xy <- specdfOUT[,c("Lon2","Lat2")]
specdfOUT_sp <- SpatialPointsDataFrame(coords=specdfOUT_xy,data=specdfOUT,proj4string = CRS("+init=epsg:32615")) 
specdfOUT_sf <- st_as_sf(specdfOUT_sp)

plotshp <- spTransform(plotshp,proj4string(specdfOUT_sp))



lapply(sort(unique(plotshp$PLOTID)),shpfile_plotloop,specdfOUT_sf,PlotShapeFile,filenumber)

# lapply(sort(unique(plotshp$RingID)),shpfile_ringloop,specdfOUT_sf,PlotShapeFile,filenumber)


  if(computer=="mac"){
    st_write(specdfOUT_sf,dsn=paste0(ProcLoc),layer=paste0("Proc",filenumber,"full"),driver="ESRI Shapefile",delete_layer=TRUE)
    print(Sys.time())
    print(filenumber) }else{print("notmac!")}
  
  if(computer=="pc"){
    st_write(specdfOUT_sf,dsn=paste0(ProcLoc,"Proc",filenumber,"full.shp"),layer=paste0("Proc",filenumber,"full"),driver="ESRI Shapefile",delete_layer=TRUE)
    print(Sys.time())
    print(filenumber)}else{print("notpc!")}

means_out <- over(plotshp,specdfOUT_sp,fn=mean)
means_out$Plot <- 1:nrow(means_out)
means_out$File <- filenumber
return(means_out)
  

}


listoffilenums <- sort(unique(as.numeric(gsub("\\D", "",list.files(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/"))))))


plotshp <- readOGR(paste0(LocalSource,"e141_poly.shp"))


system.time(out_df4 <- rbindlist(lapply(listoffilenums[c(4)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp)));beep(2)

system.time(out_df2 <- rbindlist(lapply(listoffilenums[c(2)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp)));beep(2)
system.time(out_df3 <- rbindlist(lapply(listoffilenums[c(3)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp)));beep(2)


system.time(out_df1_2_3_4 <- rbindlist(lapply(listoffilenums[c(1,2,3,4)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp)));beep(2)
system.time(out_df5_6_7_8 <- rbindlist(lapply(listoffilenums[c(5,6,7,8)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp)));beep(2)
system.time(out_df <- rbindlist(lapply(listoffilenums[c(17)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp)));beep(2)



out_df_xy <- out_df[,c("Lon2","Lat2")]
out_df_sp <- SpatialPointsDataFrame(coords=out_df_xy,data=out_df,proj4string = CRS("+init=epsg:32615")) 
plot(out_df_sp)
plot(plotshp,col="blue",add=T)
proc <- readOGR("/Volumes/HyperDrive/JaneProc/Proc4816.shp")
plot(proc,col="red",add=T)

system.time(tst <- readOGR("/Volumes/HyperDrive/JaneProc/Proc1024Plot69.shp"))
plot(tst)
spplot(tst,"band10")





