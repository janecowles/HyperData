#### Code to orthorectify UAV captured hyperspectral data (Headwall Nano)
#### Code is specific toreturning cropped experimental plots for analyses in addition to returning full orthorectified data cube
#### but can be altered to just orthorectify and return data
#### Jane M Cowles
#### jcowles@umn.edu
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
library(ggplot2)
library(data.table)

#for parallel processing
library(parallel)

#for making a beep noise
library(beepr)


##### I can update this later with Sys.info()[['sysname']] -- the biggest hurdle is the st_write where you need to include the shapefile name as part of the dsn in windows but not in mac.

# #directories for MY MAC
computer <- 'mac'
RemoteSenDataLoc <- "/Volumes/HyperDrive/Google Drive/remote sensing data/"
LocalSource <- "~/Dropbox/UMN Postdoc/Ortho_Proc/" #imu and biocon shapefile (polygons)
ProcLoc <- "/Volumes/HyperDrive/JaneProc/"
VisLoc <- "/Volumes/HyperDrive/BioCON10Aug2018/"
#directories for ISBELL PC
computer <- 'pc'
# RemoteSenDataLoc <- "F:/RemoteSensing/"
LocalSource <- "~/Ortho_Proc/" #imu and biocon shapefile (polygons)
ProcLoc <- "F:/BB_30July/"
VisLoc <- "F:/BioCON10Aug--VISUAL/"
RemoteSenDataLoc <- "F:/RemoteSensing/BigBio-eDNA-2019/" #2019July30BigBioFab1/100062_bf_2019_07_30_16_05_48
FolderLoc <- "BigBio-FAB1-2019-07-30-Hyper/100062_bf_2019_07_30_16_05_48/"
CoordSystem <- 32615

no_cores <- detectCores()-1

degtorad <- function(x){x*pi/180} #convert degrees to radians for trig functions in R


framecombine <- function(x){fi <- fread(paste0(RemoteSenDataLoc,FolderLoc,x));return(fi)}

framelist <- list.files(paste0(RemoteSenDataLoc,FolderLoc),pattern = "frameIndex")
framematch <- rbindlist(lapply(framelist,framecombine))


#coords.epsg = "32615" is utm 15n which is what I want here. "4326" is lat lon

imu <- read.delim(paste0(RemoteSenDataLoc,FolderLoc,"imu_gps.txt"))

imu <- imu[!is.na(imu$Alt),]
overallIMUmin <- min(imu$Alt)
imuxy <- imu[,c("Lon","Lat")]
imusp <- SpatialPointsDataFrame(coords=imuxy,data=imu,proj4string = CRS("+init=epsg:4326")) 
imusp <- spTransform(imusp,"+init=epsg:32615")
imu <- as.data.frame(imusp[,!names(imusp) %in% c("Lon","Lat")])
latMinIMU <- imu$Lat[imu$Alt==min(imu$Alt)]
lonMinIMU <- imu$Lon[imu$Alt==min(imu$Alt)]

setDT(imu,key="Timestamp")
setDT(framematch,key="Timestamp")
imu.framematch <- imu[framematch,roll=T]

rm(imu)
rm(imuxy)
rm(imusp)

bandtowave <- read.csv(paste0(LocalSource,"BandNumWavelength.csv"))

plotshp <- readOGR(paste0(LocalSource,"CORRECTE120/BigBio Shape/BB_Shape_JC.shp"))

#don't need to do this again
imu.framexy <- imu.framematch[,c("Lon","Lat")]
imu.framesp <- SpatialPointsDataFrame(coords=imu.framexy,data=imu.framematch,proj4string = CRS("+init=epsg:32615"))
# imu.framesp <- spTransform(imu.framesp,"+init=epsg:32615")

dem1m <- readGDAL(paste0(LocalSource,"BBFAB 1m DEM/dem_1m_m.bil"))
dem1m <- spTransform(dem1m,"+init=epsg:32615")
dem_rel <- crop(dem1m,extent(imu.framesp)+40)
class(dem_rel)
dem_sf <- st_as_sf(dem_rel)
# st_write(dem_sf,dsn=paste0(LocalSource,"dem_BBFAB.shp"),layer="dem_BBFAB.shp",driver="ESRI Shapefile",delete_layer=TRUE)
rm(dem1m)
# dem_rel <- readOGR(paste0(LocalSource,"dem_BBFAB.shp"))
# dem_rel <- spTransform(dem_rel,"+init=epsg:32615")

rast <- raster()
extent(rast) <- extent(dem_rel) 
ncol(rast) <- round((extent(dem_rel)@xmax-extent(dem_rel)@xmin)/2)
nrow(rast) <- round((extent(dem_rel)@ymax-extent(dem_rel)@ymin)/2)

# And then ... rasterize it! This creates a grid version 
# of your points using the cells of rast, values from the IP field:
dem_rast <- rasterize(dem_rel, rast, dem_rel$band1, fun=mean) 
crs(dem_rast)<-crs(dem_rel)
minAlt_dem_atminIMU <- raster::extract(dem_rast,SpatialPoints(cbind(lonMinIMU,latMinIMU)),buffer=2,fun=mean,na.rm=T)
plot(dem_rast)





imu_proc <- function(imu.datafile,FOVAngle,GroundLevel=0,minAlt_dem_atminIMU,degree=TRUE,coords.epsg,dem_rast,YawCorrFactor=0,PitchCorrFactor=0,RollCorrFactor=0){
  
 tmp <- imu.datafile
 # tmp<-imu.framematch
 sp.tmpfordem <- SpatialPointsDataFrame(coords=tmp[,c("Lon","Lat")],data=tmp,proj4string = CRS(paste0("+init=epsg:",coords.epsg)))
sp.tmpTRfordem <- spTransform(sp.tmpfordem,"+init=epsg:32615")
imuTR <- as.data.frame(sp.tmpTRfordem[,!(names(sp.tmpTRfordem)%in%c("Lon","Lat"))])

 dem_alt <- raster::extract(dem_rast,imuTR[,c("Lon","Lat")],df=T)
 dem_alt$Frame. <- imuTR$Frame.
 dem_alt$DEMAlt <- dem_alt$layer
 tmp <- merge(imuTR,dem_alt,by="Frame.")
 # tmp <- merge(imu.framematch,dem_alt,by="Frame.")
 tmp$AdjGroundLevel <- GroundLevel+(tmp$DEMAlt-minAlt_dem_atminIMU)
 tmp$HeightAboveGround <- tmp$Alt-tmp$AdjGroundLevel
 tmp$HeightAboveGroundOLD <- tmp$Alt-GroundLevel
 
    if(degree==TRUE){
    tmp$FOVAngle <- degtorad(FOVAngle)
    tmp$Roll <- degtorad(tmp$Roll)
    tmp$Pitch <- degtorad(tmp$Pitch)
    tmp$Yaw <- degtorad(tmp$Yaw)
    }
 
 #at each point, this is the FOV
  tmp$UncorrectedFOVmeters <- 2*tmp$HeightAboveGround*tan(tmp$FOVAngle/2)
  
  ############## CORRECTIONS/OFFSETS

  tmp$Roll <- tmp$Roll + RollCorrFactor
  tmp$Pitch <- tmp$Pitch + PitchCorrFactor
  tmp$Yaw <- tmp$Yaw + YawCorrFactor
  
  tmp$RollCorrFactor <-  RollCorrFactor
  tmp$PitchCorrFactor <-  PitchCorrFactor
  tmp$YawCorrFactor <-  YawCorrFactor
  
  # TESTING (ignore)
  # tmp<-data.frame(nrow=1)
  # tmp$Roll <- degtorad(5)
  # tmp$Pitch <- degtorad(5)
  # tmp$Yaw <- degtorad(45)
  # tmp$HeightAboveGround<-50
  # tmp$FOVAngle <- degtorad(15.9619)


  
  #ROLL -- Negative sign added 1/3/2019 in order to account for the reversal of what happens to the drone and which way the sensor points!
   tmp$Rolloffset <- -tmp$HeightAboveGround*tan(tmp$Roll) #parallel to the frame's long direction (as specified by yaw), therefore:
   tmp$RollYComponent <- -tmp$Rolloffset*sin(tmp$Yaw) ###(FOR NOW REMOVING THIS NEGATIVE SIGN TO SEE HOW THINGS LOOK!!!!!)#### THAT WAS WRONG PUTTING IT BACK IN
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
tmp$YawXOffset<- -0.5*tmp$RollCorrectedFOVmeters*cos(tmp$Yaw) #this is negative of y
  
sp.tmp <- SpatialPointsDataFrame(coords=tmp[,c("Lon","Lat")],data=tmp,proj4string = CRS("+init=epsg:32615"))

sp.out <- sp.tmp[,!(names(sp.tmp)%in%c("Lon","Lat"))]

out <- data.frame(sp.out)

return(out)


}

#note I have min(imu$Alt) -- only works because I have the full imu loaded (instead of the imu file of just the positions with images taken)
Proc_IMU <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast)





spacing_fun <- function(pix_i,matchimu,specdfframe){
  specdfpix <- specdfframe[specdfframe$x==pix_i,]

  
# if things are numbered the reverse of what they are (everything in these files is completely turned around because frame = 2000 is Y=0 so x goes 640->0), then needs to be switching +/- of the yaw stuff.
  specdfpix$Lat2<-matchimu$Lat+matchimu$TotMidPointCorrectionY-matchimu$YawYOffset+(pix_i*2*matchimu$YawYOffset/640)
  specdfpix$Lon2<-matchimu$Lon+matchimu$TotMidPointCorrectionX-matchimu$YawXOffset+(pix_i*2*matchimu$YawXOffset/640)
  specdfpix$Heading <- matchimu$Heading

  return(specdfpix)
}


byframe_corr <- function(framex,sdf,ProcessedIMU){
   specdfframe <- sdf[sdf$frame==framex,]
    matchimu <- ProcessedIMU[ProcessedIMU$Frame.==framex,]
    # specdfframebypix <- rbindlist(lapply(sort(unique(specdfframe$x)),spacing_fun,matchimu,specdfframe))
    
    specdfframe$Lat2 <- base::seq(from=matchimu$Lat+matchimu$TotMidPointCorrectionY-matchimu$YawYOffset,to=matchimu$Lat+matchimu$TotMidPointCorrectionY+matchimu$YawYOffset,length.out=nrow(specdfframe))
    specdfframe$Lon2 <- base::seq(from=matchimu$Lon+matchimu$TotMidPointCorrectionX-matchimu$YawXOffset,to=matchimu$Lon+matchimu$TotMidPointCorrectionX+matchimu$YawXOffset,length.out=nrow(specdfframe))    
    
    return(specdfframe)
}


shpfile_plotloop <-function(plotnum,specdfOUT_sf,PlotShapeFile,filenumber,computer=computer,ProcLoc=ProcLoc){
  plottmp <- subset(PlotShapeFile,Plot==plotnum)
  # plot_sf <- st_as_sf(plottmp)
  #for now, I want to test if I can see the plot edges as a square within the buffered region
  plot_buff <- gBuffer(plottmp,width = -2)
  plot_sf <- st_as_sf(plot_buff)
  st_crs(plot_sf)<-st_crs(specdfOUT_sf)
  suppressWarnings(plot_clip <- st_intersection(specdfOUT_sf,plot_sf))

  if(dim(plot_clip)[1]>0){
    print(plotnum)
  #   if(computer=="mac"){
  # st_write(plot_clip,dsn=paste0(ProcLoc),layer=paste0("Final",filenumber,"Plot",plotnum),driver="ESRI Shapefile",delete_layer=TRUE)}
  # 
  # if(computer=="pc"){
    st_write(plot_clip,dsn=paste0(ProcLoc,"Final",filenumber,"Plot",plotnum,".shp"),layer=paste0("Final",filenumber,"Plot",plotnum),driver="ESRI Shapefile",update = TRUE)}
  # }
}


ortho_fun <- function(filenumber,ProcessedIMU,PlotShapeFile,bandtowave){
  system.time(orig_sp <-caTools::read.ENVI(paste0(RemoteSenDataLoc,FolderLoc,"raw_",filenumber)))
  init.dim <- dim(orig_sp)
  dim(orig_sp) <- c(dim(orig_sp)[1]*dim(orig_sp)[2],dim(orig_sp)[3])
  sdf <- as.data.frame(orig_sp)
  #rep x 00000, 11111, 2222 -- columns (aka the width of the pushbroom)
  sdf$x <- rep((1:init.dim[2])-1, each = init.dim[1])

  # rep y 0 1 2 3 4 #individual scans
  sdf$y <- rep((1:init.dim[1])-1,init.dim[2])
  orig_sp <- NULL

  ### change colnames using the bandname to wavelength csv I made
  colnames(sdf)[1:272]<-paste0("nm",bandtowave$Wavelength)
  
  #### frame is filenumber+y
  
    sdf$frame <- filenumber+sdf$y
    
    #to fill in with the next function.
    sdf$Lat2 <- NA
    sdf$Lon2 <- NA
    sdf$Heading <- NA
    
    Sys.time()
    cl<-makeCluster(no_cores)
    clusterExport(cl,c("rbindlist","spacing_fun"))
    specdfOUT<- rbindlist(parLapply(cl,sort(unique(sdf$frame)),byframe_corr,sdf,ProcessedIMU))
    stopCluster(cl)
    Sys.time()
    
    
    specdfOUT$StartFrame <- filenumber
    specdfOUT_xy <- specdfOUT[,c("Lon2","Lat2")]
    specdfOUT_sp <- SpatialPointsDataFrame(coords=specdfOUT_xy,data=specdfOUT,proj4string = CRS("+init=epsg:32615")) 
    specdfOUT_sf <- st_as_sf(specdfOUT_sp)



plotshp <- spTransform(plotshp,proj4string(specdfOUT_sp))

# #comment this out for faster test runs
cl2<-makeCluster(no_cores)
clusterExport(cl2,c("st_crs","st_crs<-","subset","st_as_sf","st_intersection","st_write","gBuffer"),envir=environment())
parLapply(cl2,sort(unique(plotshp$Plot)),shpfile_plotloop,specdfOUT_sf,PlotShapeFile,filenumber,computer,ProcLoc)
stopCluster(cl2)


  # if(computer=="mac"){
  #   st_write(specdfOUT_sf,dsn=paste0(ProcLoc),layer=paste0("NEWProc",filenumber,"full"),driver="ESRI Shapefile",delete_layer=TRUE)
  #   print(Sys.time())
  #   print(filenumber) }else{print("notmac!")}
  # 
  # if(computer=="pc"){

### 11 april - temporary comment-out to see how much time is saved. for 24510 it took 3464.97 (TOTAL) WITH this in. Without? 400. 22510 was 1000 without, with??
     st_write(specdfOUT_sf,dsn=paste0(ProcLoc,"Final",filenumber,"full.shp"),layer=paste0("Final",filenumber,"full"),driver="ESRI Shapefile",update = TRUE)
     

    print(Sys.time())
    print(filenumber)
    # }else{print("notpc!")}

#commenting out for faster test runs
means_out <- over(plotshp,specdfOUT_sp,fn=mean)
means_out$Plot <- 1:nrow(means_out)
means_out$File <- filenumber

cv_out <- over(plotshp,specdfOUT_sp,fn=cv,na.rm=T)
colnames(cv_out)<- paste(colnames(cv_out),"cv",sep="_")
cv_out$Plot <- 1:nrow(cv_out)
cv_out$File <- filenumber

tot_out <-cbind(means_out,cv_out)

# tot_out <- data.frame(lol=c(1,2,3,4))
return(tot_out)
  

}

#froze after 62 so starting 63 now
listoffilenums <- sort(unique(as.numeric(gsub("\\D", "",list.files(paste0(RemoteSenDataLoc,FolderLoc))))))
system.time(out_dfTEST <- rbindlist(lapply(listoffilenums[c(1:9)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave)))

str(out_dfTEST)


Rprof("file.out")
system.time(out_dfTEST <- rbindlist(lapply(listoffilenums[c(4)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave)))
Rprof(NULL)

summaryRprof("file.out")

Proc_IMU2816 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = .45, RollCorrFactor = -0.06,PitchCorrFactor = 0.005);system.time(out_df2816 <- rbindlist(lapply(listoffilenums[c(3)],ortho_fun,ProcessedIMU=Proc_IMU2816,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(3)
write.csv(out_df2816,"~/out_df2816.csv")
rm(list=c("Proc_IMU2816","out_df2816"))

Rprof("file.out")

Proc_IMU9200 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.16, RollCorrFactor = 0.01,PitchCorrFactor = -0.04);system.time(out_df9200 <- rbindlist(lapply(listoffilenums[c(9)],ortho_fun,ProcessedIMU=Proc_IMU9200,PlotShapeFile=plotshp,bandtowave=bandtowave)))

Rprof(NULL)
summaryRprof("file.out")

write.csv(out_df9200,"~/out_df9200.csv")
rm(list=c("Proc_IMU9200","out_df9200"))


Proc_IMU1024 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.4, RollCorrFactor = 0.01,PitchCorrFactor = -0.04);system.time(out_df1024 <- rbindlist(lapply(listoffilenums[c(2)],ortho_fun,ProcessedIMU=Proc_IMU1024,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df1024,"~/out_df1024.csv")
rm(list=c("Proc_IMU1024","out_df1024"))

Proc_IMU22510 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = -0.028,PitchCorrFactor = 0.015);system.time(out_df22510 <- rbindlist(lapply(listoffilenums[c(18)],ortho_fun,ProcessedIMU=Proc_IMU22510,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df22510,"~/out_df22510.csv")
rm(list=c("Proc_IMU22510","out_df22510"))



Proc_IMU24510 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = 0.025,PitchCorrFactor = -0.025);system.time(out_df24510 <- rbindlist(lapply(listoffilenums[c(19)],ortho_fun,ProcessedIMU=Proc_IMU24510,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df24510,"~/out_df24510.csv")
rm(list=c("Proc_IMU24510","out_df24510"))

Proc_IMU45229 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = -0.045,PitchCorrFactor = 0.02);system.time(out_df45229 <- rbindlist(lapply(listoffilenums[c(31)],ortho_fun,ProcessedIMU=Proc_IMU45229,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df45229,"~/out_df45229.csv")
rm(list=c("Proc_IMU45229","out_df45229"))

Proc_IMU39773 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = 0.04,PitchCorrFactor = -0.03);system.time(out_df39773 <- rbindlist(lapply(listoffilenums[c(28)],ortho_fun,ProcessedIMU=Proc_IMU39773,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df39773,"~/out_df39773.csv")
rm(list=c("Proc_IMU39773","out_df39773"))



Proc_IMU37901 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = -0.04,PitchCorrFactor = 0.02);system.time(out_df37901 <- rbindlist(lapply(listoffilenums[c(27)],ortho_fun,ProcessedIMU=Proc_IMU37901,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df37901,"~/out_df37901.csv")
rm(list=c("Proc_IMU37901","out_df37901"))

Proc_IMU41565 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = 0.025,PitchCorrFactor = 0.01);system.time(out_df41565 <- rbindlist(lapply(listoffilenums[c(29)],ortho_fun,ProcessedIMU=Proc_IMU41565,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df41565,"~/out_df41565.csv")
rm(list=c("Proc_IMU41565","out_df41565"))


Proc_IMU36109 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = -0.045,PitchCorrFactor = -0.02);system.time(out_df36109 <- rbindlist(lapply(listoffilenums[c(26)],ortho_fun,ProcessedIMU=Proc_IMU36109,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df36109,"~/out_df36109.csv")
rm(list=c("Proc_IMU36109","out_df36109"))

Proc_IMU34317 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = 0.025,PitchCorrFactor = 0.005);system.time(out_df34317 <- rbindlist(lapply(listoffilenums[c(25)],ortho_fun,ProcessedIMU=Proc_IMU34317,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df34317,"~/out_df34317.csv")
rm(list=c("Proc_IMU34317","out_df34317"))

### THESE NEXT (TOMORROW MORNING!)
Proc_IMU4832 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.48, RollCorrFactor = -0.032,PitchCorrFactor = 0.055);system.time(out_df4832 <- rbindlist(lapply(listoffilenums[c(5)],ortho_fun,ProcessedIMU=Proc_IMU4832,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df4832,"~/out_df4832.csv")
rm(list=c("Proc_IMU4832","out_df4832"))

Proc_IMU6912 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.3, RollCorrFactor = -0.02,PitchCorrFactor = -0.085);system.time(out_df6912 <- rbindlist(lapply(listoffilenums[c(7)],ortho_fun,ProcessedIMU=Proc_IMU6912,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df6912,"~/out_df6912.csv")
rm(list=c("Proc_IMU6912","out_df6912"))


Proc_IMU12239 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = -0.01,PitchCorrFactor = 0.05);system.time(out_df12239 <- rbindlist(lapply(listoffilenums[c(11)],ortho_fun,ProcessedIMU=Proc_IMU12239,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df12239,"~/out_df12239.csv")
rm(list=c("Proc_IMU12239","out_df12239"))


Proc_IMU16847 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = 0.015,PitchCorrFactor = -0.02);system.time(out_df16847 <- rbindlist(lapply(listoffilenums[c(15)],ortho_fun,ProcessedIMU=Proc_IMU16847,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df16847,"~/out_df16847.csv")
rm(list=c("Proc_IMU16847","out_df16847"))


Proc_IMU20718 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = -0.05, RollCorrFactor = -0.015,PitchCorrFactor = -0.02);system.time(out_df20718 <- rbindlist(lapply(listoffilenums[c(17)],ortho_fun,ProcessedIMU=Proc_IMU20718,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df20718,"~/out_df20718.csv")
rm(list=c("Proc_IMU20718","out_df20718"))

Proc_IMU26510 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = 0.015,PitchCorrFactor = 0.003);system.time(out_df26510 <- rbindlist(lapply(listoffilenums[c(20)],ortho_fun,ProcessedIMU=Proc_IMU26510,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df26510,"~/out_df26510.csv")
rm(list=c("Proc_IMU26510","out_df26510"))

Proc_IMU28237 <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coord.epsg=CoordSystem,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0, RollCorrFactor = -0.03,PitchCorrFactor = -0.02);system.time(out_df28237 <- rbindlist(lapply(listoffilenums[c(21)],ortho_fun,ProcessedIMU=Proc_IMU28237,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
write.csv(out_df28237,"~/out_df28237.csv")
rm(list=c("Proc_IMU28237","out_df28237"))

plot(c(1,1))
