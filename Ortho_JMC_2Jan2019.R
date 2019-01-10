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
library(latticeExtra)

#for parallel processing
library(parallel)

#for making a beep noise
library(beepr)


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
imuxy <- imu[,c("Lon","Lat")]
imusp <- SpatialPointsDataFrame(coords=imuxy,data=imu,proj4string = CRS("+init=epsg:4326")) 
imusp <- spTransform(imusp,"+init=epsg:32615")
imu <- as.data.frame(imusp[,!names(imusp) %in% c("Lon","Lat")])
latMinIMU <- imu$Lat[imu$Alt==min(imu$Alt)]
lonMinIMU <- imu$Lon[imu$Alt==min(imu$Alt)]
rm(imu)#keep things tidy
rm(imuxy)
rm(imusp)
imu.framematch <- read.csv(paste0(LocalSource,"FramexIMU.csv"))

bandtowave <- read.csv(paste0(LocalSource,"BandNumWavelength.csv"))

plotshp <- readOGR(paste0(LocalSource,"e141_poly.shp"))

#don't need to do this again
# imu.framexy <- imu.framematch[,c("Lon","Lat")]
# imu.framesp <- SpatialPointsDataFrame(coords=imu.framexy,data=imu.framematch,proj4string = CRS("+init=epsg:4326")) 
# imu.framesp <- spTransform(imu.framesp,"+init=epsg:32615")

# dem1m <- readGDAL(paste0(LocalSource,"DEM 1m/dem_1m_m.bil"))
# dem1m <- spTransform(dem1m,"+init=epsg:32615")
# dem_rel <- crop(dem1m,extent(imu.framesp)+10)
# class(dem_rel)
# dem_sf <- st_as_sf(dem_rel)
# st_write(dem_sf,dsn=paste0(ProcLoc,"dem_BIOCON.shp"),layer="dem_BIOCON.shp",driver="ESRI Shapefile",delete_layer=TRUE)

dem_rel <- readOGR(paste0(LocalSource,"dem_BIOCON.shp"))
dem_rel <- spTransform(dem_rel,"+init=epsg:32615")

rast <- raster()
extent(rast) <- extent(dem_rel) 
ncol(rast) <- round((extent(dem_rel)@xmax-extent(dem_rel)@xmin)/2)
nrow(rast) <- round((extent(dem_rel)@ymax-extent(dem_rel)@ymin)/2)

# And then ... rasterize it! This creates a grid version 
# of your points using the cells of rast, values from the IP field:
dem_rast <- rasterize(dem_rel, rast, dem_rel$band1, fun=mean) 
crs(dem_rast)<-crs(dem_rel)
minAlt_dem_atminIMU <- raster::extract(dem_rast,SpatialPoints(cbind(lonMinIMU,latMinIMU)),buffer=2,fun=mean)


### biocon vis - read in and crop
vis <- raster(paste0(VisLoc,"ortho_biocon_visual.tif"))


# visproj <- projectRaster(vis,crs="+init=epsg:32615 +proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

system.time(ring1vis <- crop(vis,extent(bbox(subset(plotshp,RingID==1)))+30))
system.time(ring2vis <- crop(vis,extent(bbox(subset(plotshp,RingID==2)))+30))
system.time(ring3vis <- crop(vis,extent(bbox(subset(plotshp,RingID==3)))+30))
system.time(ring4vis <- crop(vis,extent(bbox(subset(plotshp,RingID==4)))+30))
system.time(ring5vis <- crop(vis,extent(bbox(subset(plotshp,RingID==5)))+30))
system.time(ring6vis <- crop(vis,extent(bbox(subset(plotshp,RingID==6)))+30))
rm(vis)




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
  
  # TESTERS
  # tmp<-data.frame(nrow=1)
  # tmp$Roll <- 1
  # tmp$Pitch <- 1
  # tmp$Yaw <- degtorad(45)
  # tmp$HeightAboveGround<-50
  
  
  #ROLL -- Negative sign added 1/3/2019 in order to account for the reversal of what happens to the drone and which way the sensor points!
   tmp$Rolloffset <- -tmp$HeightAboveGround*tan(tmp$Roll) #parallel to the frame's long direction (as specified by yaw), therefore:
   tmp$RollYComponent <- -tmp$Rolloffset*sin(tmp$Yaw) ###(FOR NOW REMOVING THIS NEGATIVE SIGN TO SEE HOW THINGS LOOK!!!!!)#### PUTTING IT BACK IN
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
  
sp.tmp <- SpatialPointsDataFrame(coords=tmp[,c("Lon","Lat")],data=tmp,proj4string = CRS("+init=epsg:32615"))

sp.out <- sp.tmp[,!(names(sp.tmp)%in%c("Lon","Lat"))]

out <- data.frame(sp.out)

return(out)


}

#note I have min(imu$Alt) -- only works because I have the full imu loaded (instead of the imu file of just the positions with images taken)
Proc_IMU <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast)

Proc_IMUYawCorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.2)

Proc_IMURollCorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,RollCorrFactor = 0.01)

Proc_IMUPitchCorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,PitchCorrFactor = 0.01)
rm(imu.framematch)


# write.csv(Proc_IMU,paste0(LocalSource,"Proc_8Jan2019.csv"),row.names = F)


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
  specdfpix$Heading <- matchimu$Heading
  
  
  ####new attempt -- minus the midpoint correction, too!
  
  return(specdfpix)
}

centerandends_corr <- function(framex,spectral.data.frame,ProcessedIMU){
   specdfframe <- spectral.data.frame[spectral.data.frame$frame==framex,]
    matchimu <- ProcessedIMU[ProcessedIMU$Frame.==framex,]
    specdfframebypix <- rbindlist(lapply(sort(unique(specdfframe$x)),spacing_fun,matchimu,specdfframe))
    return(specdfframebypix)
}


shpfile_plotloop <-function(plotnum,specdfOUT_sf,PlotShapeFile,filenumber,computer=computer,ProcLoc=ProcLoc){
  plottmp <- subset(PlotShapeFile,PLOTID==plotnum)
  plot_sf <- st_as_sf(plottmp)
  #for now, I want to test if I can see the plot edges as a square within the buffered region
  # plot_buff <- gBuffer(subset(PlotShapeFile,PLOTID==plotnum),width = 0.5)
  # plot_sf <- st_as_sf(plot_buff)
  st_crs(plot_sf)<-st_crs(specdfOUT_sf)
  suppressWarnings(plot_clip <- st_intersection(specdfOUT_sf,plot_sf))

  if(dim(plot_clip)[1]>0){
    print(plotnum)
    if(computer=="mac"){
  st_write(plot_clip,dsn=paste0(ProcLoc),layer=paste0("BUFFProc",filenumber,"Plot",plotnum),driver="ESRI Shapefile",delete_layer=TRUE)}

  if(computer=="pc"){
    st_write(plot_clip,dsn=paste0(ProcLoc,"BUFFProc",filenumber,"Plot",plotnum,".shp"),layer=paste0("BUFFProc",filenumber,"Plot",plotnum),driver="ESRI Shapefile",delete_layer=TRUE)}
  }
}

#######################################################################
#############to run smaller test runs
ProcessedIMU=Proc_IMU
PlotShapeFile=plotshp

filenumber <- 34317
orig_sp <-readGDAL(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/raw_",filenumber))
spectral.data.frame <- as.data.frame(orig_sp)
orig_sp <- NULL
colnames(spectral.data.frame)[1:272]<-paste0("nm",bandtowave$Wavelength)
spectral.data.frame$frame <- filenumber+max(spectral.data.frame$y)-spectral.data.frame$y
sub34317 <- spectral.data.frame[spectral.data.frame$frame%in%c(35000:35500),]

filenumber <- 2816
orig_sp <-readGDAL(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/raw_",filenumber))
spectral.data.frame <- as.data.frame(orig_sp)
orig_sp <- NULL
colnames(spectral.data.frame)[1:272]<-paste0("nm",bandtowave$Wavelength)
spectral.data.frame$frame <- filenumber+max(spectral.data.frame$y)-spectral.data.frame$y
sub2816 <- spectral.data.frame[spectral.data.frame$frame%in%c(3100:4600),]
rm(spectral.data.frame)


ortho_funSUB <- function(subdataframe,ProcessedIMU,PlotShapeFile,bandtowave,framesofinterest){
  spectral.data.frame<-subdataframe
  filenumber<-min(framesofinterest)
  spectral.data.frame$Lat2 <- NA
  spectral.data.frame$Lon2 <- NA
  spectral.data.frame$Heading <- NA
  
  Sys.time()
  cl<-makeCluster(no_cores)
  clusterExport(cl,c("rbindlist","spacing_fun"))
  specdfOUT<- rbindlist(parLapply(cl,sort(unique(spectral.data.frame$frame)),centerandends_corr,spectral.data.frame,ProcessedIMU))
  stopCluster(cl)
  Sys.time()
  
  
  specdfOUT$StartFrame <- filenumber
  specdfOUT_xy <- specdfOUT[,c("Lon2","Lat2")]
  specdfOUT_sp <- SpatialPointsDataFrame(coords=specdfOUT_xy,data=specdfOUT,proj4string = CRS("+init=epsg:32615")) 
  # specdfOUT_sf <- st_as_sf(specdfOUT_sp)
  
  
  rast <- raster()
  extent(rast) <- extent(specdfOUT_sp) 
  ncol(rast) <- round(640/3)
  nrow(rast) <- round(500/3)
  
  # And then ... rasterize it! This creates a grid version 
  # of your points using the cells of rast, values from the IP field:
  rast_specdfOUT_sp <- rasterize(specdfOUT_sp, rast, specdfOUT_sp$nm540.751, fun=mean) 
  crs(rast_specdfOUT_sp)<-crs(specdfOUT_sp)
  
  
  return(rast_specdfOUT_sp)
  
}


Proc_IMU <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast)

# Proc_IMUYawCorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.5)
# 
# Proc_IMURollCorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,RollCorrFactor = -0.02)
# 
# Proc_IMUPitchCorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,PitchCorrFactor = -0.02)

system.time(rast_3100<-ortho_funSUB(sub2816,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave,framesofinterest = c(3100:4600)));beep(2)
# system.time(rast_3100YAW<-ortho_funSUB(sub2816,ProcessedIMU=Proc_IMUYawCorr,PlotShapeFile=plotshp,bandtowave=bandtowave,framesofinterest = c(3100:4600)));beep(2)
# system.time(rast_3100ROLL<-ortho_funSUB(sub2816,ProcessedIMU=Proc_IMURollCorr,PlotShapeFile=plotshp,bandtowave=bandtowave,framesofinterest = c(3100:4600)));beep(2)
# system.time(rast_3100PITCH<-ortho_funSUB(sub2816,ProcessedIMU=Proc_IMUPitchCorr,PlotShapeFile=plotshp,bandtowave=bandtowave,framesofinterest = c(3100:4600)));beep(2)



Proc_IMUMULTICorr <- imu_proc(imu.datafile = imu.framematch,GroundLevel=overallIMUmin,FOVAngle = 15.9619, degree=T,coords.epsg=4326,minAlt_dem_atminIMU=minAlt_dem_atminIMU,dem_rast=dem_rast,YawCorrFactor = 0.5, RollCorrFactor = -.07, PitchCorrFactor = 0)

system.time(rast_3100MULTI<-ortho_funSUB(sub2816,ProcessedIMU=Proc_IMUMULTICorr,PlotShapeFile=plotshp,bandtowave=bandtowave,framesofinterest = c(3100:4600)));beep(2)


# system.time(proc_35000<-ortho_funSUB(sub34317,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave,framesofinterest = c(35000:35500)));beep(2)



# plot(ring2vis)
# ring2rel <- crop(ring2vis,extent(rast_3100)+5)

plot(ring2rel)

breakpoints <- c(minValue(rast_3100),minValue(rast_3100)+150,minValue(rast_3100)+250,maxValue(rast_3100))
mycol <- rgb(0, 0, 255, max = 255, alpha = 5, names = "blue50")
plot(rast_3100,breaks=breakpoints,col=c(mycol,"yellow","red"),add=T)
plot(rast_3100MULTI,breaks=breakpoints,col=c(mycol,"blue","darkblue"),add=T)



plot(ring2rel)
plot(rast_3100,breaks=breakpoints,col=c(mycol,"yellow","red"),add=T)
plot(rast_3100YAW,breaks=breakpoints,col=c(mycol,"white","darkblue"),add=T)
plot(ring2rel)
plot(rast_3100,breaks=breakpoints,col=c(mycol,"yellow","red"),add=T)
plot(rast_3100PITCH,breaks=breakpoints,col=c(mycol,"black","grey"),add=T)
plot(ring2rel)
plot(rast_3100,breaks=breakpoints,col=c(mycol,"yellow","red"),add=T)
plot(rast_3100ROLL,breaks=breakpoints,col=c(mycol,"white","grey"),add=T)

plot(rast_3100,breaks=breakpoints,col=colors)
# plot(rast_3100,add=T,col=c("red","orange"))
plot(plotshp,add=T)
plot(proc_3100,add=T)
spplot(ring2vis)+as.layer(spplot(rast_3100))
spplot(rast_3100)+as.layer(spplot(ring2vis),under=TRUE)


plot(subset(plotshp,RingID==2))
spplot(rast_3100)+as.layer(subset(plotshp,RingID==2))

##############################################end of test area
##############################################
##############################################



ortho_fun <- function(filenumber,ProcessedIMU,PlotShapeFile,bandtowave){
  orig_sp <-readGDAL(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/raw_",filenumber))
  spectral.data.frame <- as.data.frame(orig_sp)
  orig_sp <- NULL

  ### change colnames using the bandname to wavelength csv I made
  colnames(spectral.data.frame)[1:272]<-paste0("nm",bandtowave$Wavelength)
  
  #### here is where I want to FLIP IT!!!!
  ###flipped version: (NEXT: Do i want this - 0.5 anymore??)
    spectral.data.frame$frame <- filenumber+max(spectral.data.frame$y)-spectral.data.frame$y
### original version
  # spectral.data.frame$frame <- filenumber+spectral.data.frame$y-0.5
    
    
    ###FOR SUBSET TESTS
    # spectral.data.frame <- sub34317
    # filenumber <- 25000
    # spectral.data.frame <- sub2816
    # filenumber <- 3100
    #run each of these with different corrections on roll,    
    spectral.data.frame$Lat2 <- NA
    spectral.data.frame$Lon2 <- NA
    spectral.data.frame$Heading <- NA
    
    Sys.time()
    cl<-makeCluster(no_cores)
    clusterExport(cl,c("rbindlist","spacing_fun"))
    specdfOUT<- rbindlist(parLapply(cl,sort(unique(spectral.data.frame$frame)),centerandends_corr,spectral.data.frame,ProcessedIMU))
    stopCluster(cl)
    Sys.time()
    
    
    specdfOUT$StartFrame <- filenumber
    specdfOUT_xy <- specdfOUT[,c("Lon2","Lat2")]
    specdfOUT_sp <- SpatialPointsDataFrame(coords=specdfOUT_xy,data=specdfOUT,proj4string = CRS("+init=epsg:32615")) 
    specdfOUT_sf <- st_as_sf(specdfOUT_sp)



plotshp <- spTransform(plotshp,proj4string(specdfOUT_sp))

#commenting out for faster test runs
# cl2<-makeCluster(no_cores)
# clusterExport(cl2,c("st_crs","st_crs<-","subset","st_as_sf","st_intersection","st_write"),envir=environment())
# parLapply(cl2,sort(unique(plotshp$PLOTID)),shpfile_plotloop,specdfOUT_sf,PlotShapeFile,filenumber,computer,ProcLoc)
# stopCluster(cl2)

### lapply(sort(unique(plotshp$PLOTID)),shpfile_plotloop,specdfOUT_sf,PlotShapeFile,filenumber)

  if(computer=="mac"){
    st_write(specdfOUT_sf,dsn=paste0(ProcLoc),layer=paste0("NEWProc",filenumber,"full"),driver="ESRI Shapefile",delete_layer=TRUE)
    print(Sys.time())
    print(filenumber) }else{print("notmac!")}
  
  if(computer=="pc"){
    st_write(specdfOUT_sf,dsn=paste0(ProcLoc,"NEWProc",filenumber,"full.shp"),layer=paste0("NEWProc",filenumber,"full"),driver="ESRI Shapefile",delete_layer=TRUE)
    print(Sys.time())
    print(filenumber)}else{print("notpc!")}

#commenting out for faster test runs
# means_out <- over(plotshp,specdfOUT_sp,fn=mean)
# means_out$Plot <- 1:nrow(means_out)
# means_out$File <- filenumber
# return(means_out)
  

}


listoffilenums <- sort(unique(as.numeric(gsub("\\D", "",list.files(paste0(RemoteSenDataLoc,"20180917/100040_bc_2018_09_17_14_48_50/"))))))



system.time(out_df2 <- rbindlist(lapply(listoffilenums[c(2)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
system.time(out_df5 <- rbindlist(lapply(listoffilenums[c(5)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)
system.time(out_df25 <- rbindlist(lapply(listoffilenums[c(25)],ortho_fun,ProcessedIMU=Proc_IMU,PlotShapeFile=plotshp,bandtowave=bandtowave)));beep(2)

