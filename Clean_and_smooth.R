#### pre processed df
library(data.table)
library(tidyverse)
library(stringr)
library(sf)
#for clustering
library(NbClust)
library(cluster)
library(factoextra)

#multipanel ggplot
library(gridExtra)

#classification 
library(randomForest)

#spectral stuff
library(spectrolab)

library(nlme)
library(car)

df <- fread("~/REFCORRECTED_DATA.csv")
# df <- merge(bb18,corrDF,by="Plot")
df <- corrDF
#remove incomplete plots (do this before removing bareground because that'll artificially 
df$PlotStartFrame <- paste(df$StartFrame,df$Plot,sep="_")
# dt <- df
setDT(df)
dt.stats<- df[,.(preremovalnumberofpixels=length(x)),by=c("Plot","StartFrame","PlotStartFrame")]

dt.stats <- dt.stats[order(-dt.stats$preremovalnumberofpixels)]
dt.stats1 <- dt.stats[!duplicated(dt.stats$Plot)]
hist(dt.stats1$preremovalnumberofpixels)
dt.stats1 <- dt.stats1[dt.stats1$preremovalnumberofpixels>1000]

includeplots<-dt.stats1$PlotStartFrame


df_keep <- df[df$PlotStartFrame%in%includeplots,]
df_keep <- merge(bb18,df_keep,all=T)

df_keep$NDVI <- (df_keep$nm800.614-df_keep$nm680.677)/(df_keep$nm800.614+df_keep$nm680.677)
# NDVI wavelengths checked from Wang et al 2016 (Alberta prairie), and they should range from -1 to 1 but I'm leaving all in at this point.
# fwrite(df_keep,"~/BB_NoDups_WITHBG10Aug2019.csv")
# df_keep <- fread("~/Documents/BB_NoDups_WITHBG10Aug2019.csv")
df_keep$AGB <- df_keep$AbvBioAnnProd
# df_keep$SuperBarePlots <- ifelse(df_keep$Plot%in%unique(df_keep$Plot[df_keep$AGB<5]),"SuperBare","NormalPlot")
# df_keep$SuperBarePlots[df_keep$Plot %in% unique(df_keep$Plot[df_keep$AGB>1000])]<-"SuperHigh"

df_keep$SuperBarePlots <- ifelse(df_keep$Plot%in%c(126,166,244,248),"SuperBare","NormalPlot")
df_keep$SuperBarePlots[df_keep$Plot %in% unique(df_keep$Plot[df_keep$AGB>350])]<-"SuperHigh"

# ggplot(df_keep[df_keep$Ring==2&df_keep$Use.for.Main=="Y",],aes(Lon2,Lat2,color=factor(CountOfSpecies)))+geom_point()+coord_fixed()+ scale_color_manual(values=c("darkorange","deepskyblue","forestgreen","darkslategrey"))+theme_minimal()
# ggsave("~/Documents/Ring2_dfkeepAKAallpointsMAIN.pdf")
unique(df_keep$Plot[df_keep$X..cover.Total.Planted.Species>75])






df_train <- df_keep[df_keep$SuperBarePlots%in%c("SuperBare","SuperHigh"),]
df_train$SuperBarePlots<-factor(df_train$SuperBarePlots)

wavelengthcols <- paste(colnames(df_train[,c(280:551)]),collapse="+")
formula1 <- as.formula(paste("SuperBarePlots~",wavelengthcols,sep=""))
system.time(rf1 <- randomForest(formula1,data = df_train))
df_e120 <- df_keep[!is.na(df_keep$AGB),]

system.time(df_e120$preds <- predict(rf1,df_e120))

table(df_e120$preds)


ggplot(df_e120,aes(preds,AGB))+geom_violin()

ggplot(df_e120,aes(SuperBarePlots,NDVI))+geom_violin()


ExpBGremove <- df_e120[df_e120$preds=="SuperHigh",]
# fwrite(ExpBGremove,"~RF_BGRemoval_BB_12AUG2019.csv")


ggplot(ExpBGremove,aes(Lon2,Lat2,color=factor(NumSp)))+geom_point()+coord_fixed()+scale_color_manual(values=c("darkorange","deepskyblue","forestgreen","darkslategrey","yellow","red"))+theme_minimal()


# ExpBGremove_r2 <- ExpBGremove[ExpBGremove$Ring==2,]
# ggplot(ExpBGremove_r2,aes(Lon2,Lat2,color=NDVI))+geom_point()
# table(ExpBGremove_r2$Plot)
# ggplot(ExpBGremove_r2,aes(Lon2,Lat2,color=NDVI))+geom_point()
setDT(df_e120)
setDT(ExpBGremove)

cols <- colnames(df_e120)[c(1:270)]
preremoval_nrows <- df_e120[,.(preremovalnumberofpixels=length(NDVI)),by=cols]
ggplot(preremoval_nrows,aes(preremovalnumberofpixels,AbvBioAnnProd))+geom_point()


cols <- colnames(ExpBGremove)[c(1:270)]
mean_nrows <- ExpBGremove[,.(numberofpixels=length(NDVI)),by=cols]
ggplot(mean_nrows,aes(numberofpixels,AbvBioAnnProd))+geom_point()

preandpost <- merge(preremoval_nrows,mean_nrows,by=cols)
ggplot(preandpost,aes(preremovalnumberofpixels,numberofpixels,color=NumSp))+geom_point()
ggplot(preandpost[preandpost$preremovalnumberofpixels>50000],aes(preremovalnumberofpixels,numberofpixels,color=factor(Plot)))+geom_point()


ggplot(preandpost,aes(preremovalnumberofpixels,numberofpixels,color=X..cover.Total.Planted.Species,size=X..cover.Total.Planted.Species))+geom_point()

fwrite(preandpost,"~/PlotLevelPixelsRemoved_BB_13Aug2019.csv")


ggsave("~/Documents/Pixelsxremoval.pdf")

dffin <- ExpBGremove[ExpBGremove$Use.for.Main=="Y",]
dffin$PlotXY <- paste(dffin$Plot,dffin$x,dffin$y,sep="_")

colnames(dffin)[94:365] <- as.numeric(str_replace(colnames(dffin)[94:365],"nm",""))

#next, I'll want to try to read in the files I made below and match them up with the true lat lon... this will be the final smoothed file.
dffin_loc <- dffin[,c("PlotXY","Lon2","Lat2")]
colnames(dffin_loc)<-c("PlotXY","Longitude","Latitude")
# ggplot(dffin_loc[1:75000,],aes(Longitude,Latitude))+geom_point()
i=1
for(i in 1:6){
dffin_tmp <- dffin[dffin$Ring==i,]
df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
df_smooth <- smooth(df_spec)
# df_norm<-normalize(df_smooth)
# smdf <- as.data.frame(df_norm)
smdf <- as.data.frame(df_smooth)
fwrite(smdf,paste0("~/Documents/SmoothRing",i,"_23June2019.csv"))
print(dim(df_smooth))
rm(dffin_tmp,df_spec,df_smooth,df_norm,smdf)
}



# 
# 
# smdf <- fread("~/Documents/SmoothRing6_18June2019.csv")
# names(df_smdf)
# mean(smdf$`400.825`)
# sd(smdf$`400.825`)/mean(smdf$`400.825`)

#### here's the start for tomorrow! Read in all the ring files as above, and go forth!
smdf1 <- fread("~/Documents/SmoothRing1_23June2019.csv")
smdf2 <- fread("~/Documents/SmoothRing2_23June2019.csv")
smdf3 <- fread("~/Documents/SmoothRing3_23June2019.csv")
smdf4 <- fread("~/Documents/SmoothRing4_23June2019.csv")
smdf5 <- fread("~/Documents/SmoothRing5_23June2019.csv")
smdf6 <- fread("~/Documents/SmoothRing6_23June2019.csv")

smdf <- rbind(smdf1,smdf2,smdf3,smdf4,smdf5,smdf6)

ggplot(smdf1,aes(Lon2,Lat2,color=NDVI))+geom_point()

names(smdf)

sm_fin <- merge(dffin_loc,smdf,by.x="PlotXY",by.y="sample_name")
sm_fin <-sm_fin[,!c("Lon2","Lat2")]
# fwrite(sm_fin,file="~/Documents/HYPER-FinalSmoothedCorrectLocs_23June2019.csv")

sm_fin <- fread("~/Documents/HYPER-FinalSmoothedCorrectLocs_23June2019.csv")
ggplot(sm_fin[sm_fin$Ring==2,],aes(Longitude,Latitude,color=NDVI))+geom_point()

