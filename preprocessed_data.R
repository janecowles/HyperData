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

df <- fread("~/Documents/REFCORRECTED_DATA.csv")
#remove incomplete plots (do this before removing bareground because that'll artificially 
df$PlotStartFrame <- paste(df$StartFrame,df$Plot,sep="_")
completeplots <- NULL
for(i in unique(df$PlotStartFrame)){
  if(nrow(df[df$PlotStartFrame==i,])>1800){completeplots <- append(completeplots,i)
  }
}


df_comp <- df[df$PlotStartFrame%in%completeplots,]


# tapply(df_comp$StartFrame,df_comp$Plot,unique)

plotcols <- df_comp[,c("Plot","StartFrame","PlotStartFrame")]
plotcols <- plotcols[!duplicated(plotcols$PlotStartFrame)]
plotcols <- plotcols[order(plotcols$StartFrame,plotcols$Plot)]
plotcols_first <- plotcols[!duplicated(plotcols$Plot),]
keepers <- plotcols_first$PlotStartFrame
df_keep <- df_comp[df_comp$PlotStartFrame%in%keepers,]
df_keep$NDVI <- (df_keep$nm800.614-df_keep$nm680.677)/df_keep$nm800.614+df_keep$nm680.677
# NDVI #s checked from Wang et al 2016 (Alberta prairie)
df_keep <- df_keep[df_keep$NDVI<1&df_keep$NDVI>-1,]
# fwrite(df_keep,"~/Documents/FINAL_NoDups_WITHBG5June2019.csv")
df_keep <- fread("~/Documents/FINAL_NoDups_WITHBG5June2019.csv")
df_keep$AGB <- df_keep$AbovegroundTotal.Biomass..g.m.2
df_keep$SuperBarePlots <- ifelse(df_keep$Plot%in%c(260,276,338,73,77,31,46),"SuperBare","NormalPlot")
df_keep$SuperBarePlots[df_keep$Plot %in% unique(df_keep$Plot[df_keep$AGB>1000])]<-"SuperHigh"
df_keep$`Plot Diversity` <- as.factor(df_keep$CountOfSpecies)
ggplot(df_keep[df_keep$Ring==2,],aes(Lon2,Lat2,color=factor(CountOfSpecies)))+geom_point()+coord_fixed()+ scale_color_manual(values=c("grey87","darkorange","deepskyblue","forestgreen","darkslategrey"))+theme_minimal()
ggsave("~/Documents/Ring2_dfkeepAKAallpoints.pdf")
unique(df_keep$Plot[df_keep$X..cover.Total.Planted.Species>75])


df_keep <- fread("~/Documents/FINAL_NoDups_WITHBG5June2019.csv")
set.seed(10)
train_plots <- sample(unique(df_keep$Plot),50)
# df_train2 <- df_keep[df_keep$Plot%in%train_plots]
df_train2 <- df_keep[df_keep$Ring==2,]

ggplot(df_train2,aes(Lon2,Lat2,color=AbovegroundTotal.Biomass..g.m.2
))+geom_point()

clusts <- kmeans(df_train2[, 94:365], 2,nstart = 20)


df_train2c <- cbind(df_train2,clusts$cluster)
colnames(df_train2c)[ncol(df_train2c)]<- "cluster"
ggplot(df_train2c[df_train2$Ring==2,],aes(Lon2,Lat2,color=cluster))+geom_point()



# full df
clusts <- kmeans(df_keep[, 94:365], 2,nstart = 20)


df_keepc <- cbind(df_keep,clusts$cluster)
colnames(df_keepc)[ncol(df_keepc)]<- "cluster"
ggplot(df_keepc[df_keepc$Ring==2,],aes(Lon2,Lat2,color=cluster))+geom_point()

ggplot(df_keep[df_keep$Ring==2,],aes(Lon2,Lat2,color=factor(CountOfSpecies)))+geom_point()

tapply(df_keep$CountOfSpecies[df_keep$Ring==2],df_keep$Plot[df_keep$Ring==2],mean)


df_train <- df_keep[df_keep$SuperBarePlots%in%c("SuperBare","SuperHigh"),]
df_train$SuperBarePlots<-factor(df_train$SuperBarePlots)

wavelengthcols <- paste(colnames(df_train[,c(94:365)]),collapse="+")
formula1 <- as.formula(paste("SuperBarePlots~",wavelengthcols,sep=""))
system.time(rf1 <- randomForest(formula1,data = df_train))
df_keep$preds <- predict(rf1,df_keep)

table(df_keep$preds)


ggplot(df_keep,aes(preds,AGB))+geom_violin()

ggplot(df_keep,aes(SuperBarePlots,NDVI))+geom_violin()


ExpBGremove <- df_keep[df_keep$preds=="SuperHigh"]
# fwrite(ExpBGremove,"~/Documents/RF_BGRemoval6June2019.csv")

ExpBGremove <- fread("~/Documents/RF_BGRemoval6June2019.csv")
ggplot(ExpBGremove[ExpBGremove$Ring==1,],aes(Lon2,Lat2,color=NDVI))+geom_point()
ggsave("~/Documents/Ring1.pdf")

ggplot(ExpBGremove[ExpBGremove$Ring==2,],aes(Lon2,Lat2,color=NDVI))+geom_point()
ggplot(ExpBGremove[ExpBGremove$Ring==2,],aes(Lon2,Lat2,color=factor(CountOfSpecies)))+geom_point()+coord_fixed()+scale_color_manual(values=c("darkorange","deepskyblue","forestgreen","darkslategrey"))+theme_minimal()
ggsave("~/Documents/Ring2_div.pdf")

ggplot(ExpBGremove[ExpBGremove$Ring==3,],aes(Lon2,Lat2,color=NDVI))+geom_point()
ggsave("~/Documents/Ring3.pdf")

ggplot(ExpBGremove[ExpBGremove$Ring==4,],aes(Lon2,Lat2,color=NDVI))+geom_point()
ggsave("~/Documents/Ring4.pdf")

ggplot(ExpBGremove[ExpBGremove$Ring==5,],aes(Lon2,Lat2,color=NDVI))+geom_point()
ggsave("~/Documents/Ring5.pdf")

ggplot(ExpBGremove[ExpBGremove$Ring==6,],aes(Lon2,Lat2,color=NDVI))+geom_point()
ggsave("~/Documents/Ring6.pdf")

# ExpBGremove_r2 <- ExpBGremove[ExpBGremove$Ring==2,]
# ggplot(ExpBGremove_r2,aes(Lon2,Lat2,color=NDVI))+geom_point()
# table(ExpBGremove_r2$Plot)
# ggplot(ExpBGremove_r2,aes(Lon2,Lat2,color=NDVI))+geom_point()


cols <- colnames(df_keep)[c(11:92,366)]
preremoval_nrows <- df_keep[,.(preremovalnumberofpixels=length(NDVI)),by=cols]
ggplot(preremoval_nrows,aes(preremovalnumberofpixels,X..cover.Total.Planted.Species))+geom_point()


cols <- colnames(ExpBGremove)[c(11:92,366)]
mean_nrows <- ExpBGremove[,.(numberofpixels=length(NDVI)),by=cols]
ggplot(mean_nrows,aes(numberofpixels,X..cover.Total.Planted.Species))+geom_point()

preandpost <- merge(preremoval_nrows,mean_nrows,by=cols)
ggplot(preandpost,aes(preremovalnumberofpixels,numberofpixels,color=X..cover.Total.Planted.Species))+geom_point()

ggplot(preandpost,aes(preremovalnumberofpixels,numberofpixels,color=X..cover.Total.Planted.Species))+geom_point()


ggsave("~/Documents/Pixelsxremoval.pdf")

dffin <- ExpBGremove[ExpBGremove$Use.for.Main=="Y",]
dffin$PlotXY <- paste(dffin$Plot,dffin$x,dffin$y,sep="_")

colnames(dffin)[94:365] <- as.numeric(str_replace(colnames(dffin)[94:365],"nm",""))

#next, I'll want to try to read in the files I made below and match them up with the true lat lon... this will be the final smoothed file.
dffin_loc <- dffin[,c("PlotXY","Lon2","Lat2")]
colnames(dffin_loc)<-c("PlotXY","Longitude","Latitude")
ggplot(dffin_loc[1:75000,],aes(Longitude,Latitude))+geom_point()

# dffin_tmp <- dffin[dffin$Ring==1,]
# df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
# df_smooth <- smooth(df_spec)
# smdf <- as.data.frame(df_smooth)
# fwrite(smdf,"~/Documents/SmoothRing1_11June2019.csv")
# rm(dffin_tmp)
# rm(df_spec)
# rm(df_smooth)
# rm(smdf)
# dffin_tmp <- dffin[dffin$Ring==2,]
# df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
# df_smooth <- smooth(df_spec)
# smdf <- as.data.frame(df_smooth)
# fwrite(smdf,"~/Documents/SmoothRing2_11June2019.csv")
# rm(dffin_tmp)
# rm(df_spec)
# rm(df_smooth)
# rm(smdf)
# 
# dffin_tmp <- dffin[dffin$Ring==3,]
# df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
# df_smooth <- smooth(df_spec)
# smdf <- as.data.frame(df_smooth)
# fwrite(smdf,"~/Documents/SmoothRing3_11June2019.csv")
# rm(dffin_tmp)
# rm(df_spec)
# rm(df_smooth)
# rm(smdf)
# 
# dffin_tmp <- dffin[dffin$Ring==4,]
# df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
# df_smooth <- smooth(df_spec)
# smdf <- as.data.frame(df_smooth)
# fwrite(smdf,"~/Documents/SmoothRing4_11June2019.csv")
# rm(dffin_tmp)
# rm(df_spec)
# rm(df_smooth)
# rm(smdf)
# 
# dffin_tmp <- dffin[dffin$Ring==5,]
# df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
# df_smooth <- smooth(df_spec)
# smdf <- as.data.frame(df_smooth)
# fwrite(smdf,"~/Documents/SmoothRing5_11June2019.csv")
# rm(dffin_tmp)
# rm(df_spec)
# rm(df_smooth)
# rm(smdf)
# 
# dffin_tmp <- dffin[dffin$Ring==6,]
# df_spec <- as.spectra(dffin_tmp,name_idx = 371,meta_idxs = c(1:93,366:370))
# df_smooth <- smooth(df_spec)
# smdf <- as.data.frame(df_smooth)
# fwrite(smdf,"~/Documents/SmoothRing6_11June2019.csv")
# rm(dffin_tmp)
# rm(df_spec)
# rm(df_smooth)
# rm(smdf)
# 
# 
# smdf <- fread("~/Documents/SmoothRing6_11June2019.csv")
# names(df_smdf)
# mean(smdf$`400.825`)
# sd(smdf$`400.825`)/mean(smdf$`400.825`)

#### here's the start for tomorrow! Read in all the ring files as above, and go forth!
smdf1 <- fread("~/Documents/SmoothRing1_11June2019.csv")
smdf2 <- fread("~/Documents/SmoothRing2_11June2019.csv")
smdf3 <- fread("~/Documents/SmoothRing3_11June2019.csv")
smdf4 <- fread("~/Documents/SmoothRing4_11June2019.csv")
smdf5 <- fread("~/Documents/SmoothRing5_11June2019.csv")
smdf6 <- fread("~/Documents/SmoothRing6_11June2019.csv")

smdf <- rbind(smdf1,smdf2,smdf3,smdf4,smdf5,smdf6)

ggplot(smdf1,aes(Lon2,Lat2,color=NDVI))+geom_point()

names(smdf)

sm_fin <- merge(dffin_loc,smdf,by.x="PlotXY",by.y="sample_name")
sm_fin <-sm_fin[,!c("Lon2","Lat2")]
# fwrite(sm_fin,file="~/Documents/HYPER-FinalSmoothedCorrectLocs_18June2019.csv")

sm_fin <- fread("~/Documents/HYPER-FinalSmoothedCorrectLocs_18June2019.csv")

# smdf5_long <- melt(smdf5,id.vars=c(1:99),measure.vars = c(100:371),variable.name="wavelength",value.name="reflectance")
# smdf5_long$wavelength<-as.numeric(as.character(smdf5_long$wavelength))
# 
# ggplot(smdf5_long[smdf5_long$wavelength==400.825,],aes(Lon2,reflectance))+geom_point()
# ggplot(smdf5_long,aes(wavelength,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")


cols <- colnames(smdf)[c(2,12:93,97)]
smdf_means <- smdf[,lapply(.SD, mean),  .SDcols = 100:371,by=cols]
names(smdf_means)

smdf_mean_long <- melt(smdf_means, id.vars = c(1:84),measure.vars = c(85:356),variable.name = "wavelength",value.name = "reflectance")

smdf_mean_long$wavelength<-as.numeric(as.character(smdf_mean_long$wavelength))

ggplot(smdf_mean_long,aes(wavelength,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")+facet_grid(CO2.Treatment~Nitrogen.Treatment)+scale_color_manual(values=c("red","purple","orange","blue"))
ggsave("~/Documents/smdf_mean_long1.pdf")
ggplot(smdf_mean_long,aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary")+facet_wrap("CountOfSpecies")+scale_color_manual(values=c("red","purple","orange","blue"))
ggsave("~/Documents/smdf_mean_long2.pdf")


cv <- function(x){sd(x)/mean(x)}
cols <- colnames(smdf)[c(2,12:93,97)]
smdf_cv <- smdf[,lapply(.SD, cv),  .SDcols = 100:371,by=cols]
names(smdf_cv)

smdf_cv$meanCV <- rowMeans(smdf_cv[,85:356])
smdf_cv$meanCVcons <- rowMeans(smdf_cv[,98:317])
ggplot(smdf_cv,aes(CountOfSpecies,meanCV))+geom_point()+geom_point(stat="summary",col="red",size=3)
ggsave("~/Documents/smdf_CV1.pdf")

ggplot(smdf_cv,aes(CountOfSpecies,meanCVcons))+geom_point()+geom_point(stat="summary",col="red",size=3)
ggsave("~/Documents/smdf_cv2.pdf")

ggplot(smdf_cv,aes(meanCVcons,AbovegroundTotal.Biomass..g.m.2))+geom_point()
ggplot(smdf_cv,aes(meanCVcons,AbovegroundTotal.Biomass..g.m.2,color=C.and.N.treatment,group=C.and.N.treatment))+geom_point()+geom_smooth(method="lm")
smdf_cv[smdf_cv$AbovegroundTotal.Biomass..g.m.2<10&smdf_cv$meanCVcons>0.5,]

smdf_cv$RingFac <- as.factor(smdf_cv$Ring)
smdf_cv$PlotFac <- as.factor(smdf_cv$Plot)
mod <- lme(AbovegroundTotal.Biomass..g.m.2~meanCVcons,random=~1|Ring/Plot,data=smdf_cv)
summary(mod)
Anova(mod)

mod <- lme(meanCVcons~CountOfSpecies,random=~1|Ring/Plot,data=smdf_cv)
summary(mod)
Anova(mod)
cor(smdf_cv$meanCVcons,smdf_cv$CountOfSpecies)
cor(smdf_cv$meanCV,smdf_cv$CountOfSpecies)

ggplot(smdf_mean_long[smdf_mean_long$CountOfSpecies==1&smdf_mean_long$wavelength<680,],aes(wavelength,reflectance,group=monospecies,color=monospecies))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/smdf_mean_longsp1.pdf")
ggplot(smdf_mean_long[smdf_mean_long$CountOfSpecies==1&smdf_mean_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("monospecies")
ggsave("~/Documents/smdf_mean_longsp2.pdf")


ggplot(smdf_mean_long[smdf_mean_long$CountOfSpecies==1&smdf_mean_long$wavelength<680,],aes(wavelength,reflectance,group=Monogroup,color=Monogroup))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/smdf_mean_longgr1.pdf")

ggplot(smdf_mean_long[smdf_mean_long$CountOfSpecies==1&smdf_mean_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("Monogroup")
ggsave("~/Documents/smdf_mean_longgr2.pdf")


####






ggplot(dffin_long,aes(wavelengthnum,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")
ggplot(dffin_long,aes(wavelengthnum,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")+facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggplot(dffin_long,aes(wavelengthnum,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary")+facet_grid(1~CountOfSpecies)

ggplot(dffin_long[dffin_long$wavelengthnum>800,],aes(wavelengthnum,reflectance,group=Plot,color=(CountOfSpecies)))+geom_point(stat="summary")
ggplot(dffin_long,aes(wavelengthnum,reflectance,group=Plot,color=factor(CountOfSpecies)))+geom_point(stat="summary")



# cv <- function(x){sd(x)/mean(x)}

plot_cv <- dffin[,lapply(.SD, cv),  .SDcols = 94:365,by=cols]
plot_cv$meanCV <- rowMeans(plot_cv[,84:355])

ggplot(plot_cv[plot_cv$Ring==5&plot_cv$Plot%in%c(252:254),],aes(CountOfSpecies,meanCV,color=factor(Plot)))+geom_point()

plot_cv[plot_cv$Plot==253,]

sd(dffin$nm400.825[dffin$Plot==253])/mean(dffin$nm400.825[dffin$Plot==253])
length(df_keep$nm400.825[df_keep$Plot==253])
#




































#############################################
df_keepMAIN <- df_keep[df_keep$Use.for.Main=="Y",]
  
ggplot(df_keep,aes(StartFrame,NDVI,color=factor(CountOfSpecies),group=factor(CountOfSpecies)))+geom_point(stat="summary")

hist(df_keepMAIN$NDVI)
hist(df_keepMAIN$nm800.614)
hist(df_keepMAIN$nm680.677)
clean_ndvi <- df_keepMAIN$NDVI
clean_ndvi <- clean_ndvi[clean_ndvi>0]
clean_ndvi <- clean_ndvi[clean_ndvi<1]
hist(clean_ndvi)
#remove NDVI <0.4
# df_biggerthan1 <- df_keepMAIN[df_keepMAIN$NDVI>1,]
df_nobg <- df_keepMAIN[df_keepMAIN$NDVI>0.4&df_keepMAIN$NDVI<1,]
# df_nobg <- df_nobg[df_nobg$NDVI<2,]
# fwrite(df_nobg,"~/Documents/FINAL_NoDups_NoBG2June2019.csv")



#### clean out by AGB = 0




df_keep0 <- df_keep[df_keep$Plot %in% unique(df_keep$Plot[df_keep$AbovegroundTotal.Biomass..g.m.2<1&df_keep$CountOfSpecies==0])]

ggplot(df_keep0,aes(Lon2,Lat2,color=NDVI))+geom_point()+facet_wrap("Plot",scales = "free")

ggplot(df_keep0,aes(Lon2,Lat2,color=AGB))+geom_point()+facet_wrap("Plot",scales = "free")

df_keep016 <- df_keep[df_keep$Plot %in% unique(df_keep$Plot[(df_keep$AGB<1&df_keep$CountOfSpecies==0)|(df_keep$AGB>800&df_keep$CountOfSpecies==16)])]

ggplot(df_keep016,aes(Lon2,Lat2,color=NDVI))+geom_point()+facet_wrap("Plot",scales = "free")

ggplot(df_keep016,aes(Lon2,Lat2,color=AGB))+geom_point()+facet_grid(CountOfSpecies~Plot,scales = "free")
keeplong <- melt(df_keep016, id.vars = c(1:92,368),measure.vars = c(94:365),variable.name = "wavelength",value.name = "reflectance")
keeplong$wavelengthnum <- as.numeric(str_replace(keeplong$wavelength,"nm",""))
ggplot(keeplong,aes(wavelengthnum,reflectance,group=Plot,color=CountOfSpecies))+geom_point(stat="summary")
ggplot(keeplong[!keeplong$Plot%in%c(61,289),],aes(wavelengthnum,reflectance,group=Plot,color=CountOfSpecies))+geom_point(stat="summary")

ggplot(df_keep016[df_keep016$CountOfSpecies==0],aes(Lon2,Lat2,color=nm540.751))+geom_point()+facet_wrap("Plot",scales="free")

testCluster <- kmeans(df_keep016[, c(113:212)], 3, nstart = 20)
fviz_cluster(testCluster,r1[,c(113:212)],geom=c("point"),pointsize = 0.5,labelsize=1)
df_keep016$CLUSTERONWAVE<-testCluster$cluster

ggplot(df_keep016[df_keep016$CLUSTERONWAVE==1,],aes(Lon2,Lat2,color=factor(CountOfSpecies)))+geom_point()+facet_wrap("Plot",scales="free")

table(factor(df_keep016$CLUSTERONWAVE),df_keep016$CountOfSpecies)
tapply(df_keep016$NDVI,factor(df_keep016$CLUSTERONWAVE),max)


tapply(df_keep$NDVI,factor(df_keep$SuperBarePlots),quantile,.75)
cols <- colnames(df_keep)[c(11:92,368:369)]
keep_means <- df_keep[,lapply(.SD, mean),  .SDcols = 94:365,by=cols]
plotlong <- melt(keep_means, id.vars = c(1:84),measure.vars = c(85:356),variable.name = "wavelength",value.name = "reflectance")
plotlong$wavelengthnum<-as.numeric(str_replace(plotlong$wavelength,"nm",""))
ggplot(plotlong,aes(wavelengthnum,reflectance,group=SuperBarePlots,color=SuperBarePlots))+geom_point(stat="summary")
ggplot(plotlong[plotlong$X..cover.Total.Planted.Species<10,],aes(wavelengthnum,reflectance,group=X..cover.Total.Planted.Species,color=X..cover.Total.Planted.Species))+geom_point(stat="summary")





#################

# df_nobg <- fread("~/Documents/FINAL_NoDups_NoBG2June2019.csv")


# plotx <- df_nobg[df_nobg$CountOfSpecies%in%c(0,16)&df_nobg$Ring%in%c(1,2),]
# 
# plotlong <- melt(plotx, id.vars = c(1:92,367),measure.vars = c(94:365),variable.name = "wavelength",value.name = "reflectance")
# plotlong$wavelengthnum <- as.numeric(str_replace(plotlong$wavelength,"nm",""))
# ggplot(plotlong,aes(wavelengthnum,reflectance,group=CountOfSpecies,color=CountOfSpecies))+geom_point(stat="summary")

plotx <- df_nobg[df_nobg$Plot%in%c(100:102),]
ggplot(plotx,aes(Lon2,Lat2,color=NDVI))+geom_point()+coord_fixed()

head(plotx)
plotx_line <- plotx[plotx$Lat2>5027867.1&plotx$Lat2<5027867.12,]
ggplot(plotx_line,aes(Lon2,Lat2,color=NDVI))+geom_point()+coord_fixed()
plotx_linelong <- melt(plotx_line, id.vars = c(1:92,367),measure.vars = c(94:365),variable.name = "wavelength",value.name = "reflectance")
plotx_linelong <- plotx_linelong[plotx_linelong$reflectance>0&plotx_linelong$reflectance<1,]
plotx_linelong$wavelengthnum <- as.numeric(str_replace(plotx_linelong$wavelength,"nm",""))
plotx_linelong$Lon2Fac <- as.factor(plotx_linelong$Lon2)
ggplot(plotx_linelong[plotx_linelong$wavelengthnum>440&plotx_linelong$wavelengthnum<660,], aes(Lon2Fac, wavelengthnum)) + geom_tile(aes(fill = reflectance)) + scale_fill_gradient(low = "white",high = "steelblue")


cv <- function(x){sd(x)/mean(x)}

cols <- colnames(df_nobg)[11:92]
plot_means <- df_nobg[,lapply(.SD, mean),  .SDcols = 94:366,by=cols]
plot_cv <- df_nobg[,lapply(.SD, cv),  .SDcols = 94:366,by=cols]
plot_cv$meanCV <- rowMeans(plot_cv[,83:354])
ggplot(plot_means,aes(NDVI,AbovegroundTotal.Biomass..g.m.2))+geom_point()

ggplot(plot_cv,aes(CountOfSpecies,meanCV))+geom_point()


names(plot_cv)

plotlong <- melt(plot_means, id.vars = c(1:82),measure.vars = c(83:354),variable.name = "wavelength",value.name = "reflectance")
plotlong$wavelengthnum <- as.numeric(str_replace(plotlong$wavelength,"nm",""))
plotlong$Lon2Fac <- as.factor(plotlong$Lon2)

ggplot(plotlong,aes(wavelengthnum,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")
ggplot(plotlong[plotlong$wavelengthnum>440&plotlong$wavelengthnum<660,],aes(wavelengthnum,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")


trt_means <- df_nobg[,lapply(.SD, mean),  .SDcols = 94:366,by=.(CountOfSpecies,C.and.N.treatment)]
trtlong <- melt(trt_means, id.vars = c(1:2),measure.vars = c(3:275),variable.name = "wavelength",value.name = "reflectance")
trtlong$wavelengthnum <- as.numeric(str_replace(trtlong$wavelength,"nm",""))
trtlong$Lon2Fac <- as.factor(trtlong$Lon2)

ggplot(trtlong[trtlong$wavelengthnum>440&trtlong$wavelengthnum<660,], aes(factor(CountOfSpecies), wavelengthnum)) + geom_tile(aes(fill = reflectance)) + scale_fill_gradient(low = "white",high = "steelblue")

ggplot(trtlong[trtlong$wavelengthnum>440&trtlong$wavelengthnum<660,], aes(factor(C.and.N.treatment), wavelengthnum)) + geom_tile(aes(fill = reflectance)) + scale_fill_gradient(low = "white",high = "steelblue")+facet_wrap("CountOfSpecies")

plot0 <- df_nobg[df_nobg$CountOfSpecies==0,]

ggplot(plot0,aes(Lon2,Lat2,color=NDVI))+geom_point()+facet_wrap("Plot",scales = "free")

r1 <- df_nobg[df_nobg$Ring==1,] 
ggplot(r1[r1$NDVI<1,],aes(Lon2,Lat2,color=NDVI))+geom_point()

testCluster <- kmeans(r1[, c("Lat2","Lon2")], 61, nstart = 20)
fviz_cluster(testCluster,r1[,c("Lat2","Lon2")],geom=c("point"),pointsize = 0.5,labelsize=1)

testCluster <- kmeans(r1[, c(94:365)], 17, nstart = 20)
fviz_cluster(testCluster,r1[,c("Lat2","Lon2")],geom=c("point"),pointsize = 0.5,labelsize=1)

r2 <- df_nobg[df_nobg$Ring==2,]
r2$AGB <-r2$AbovegroundTotal.Biomass..g.m.2

g1 <- ggplot(r2,aes(Lon2,Lat2,color=factor(NDVI)))+geom_point()+coord_fixed()
g2 <- ggplot(r2,aes(Lon2,Lat2,color=factor(AGB)))+geom_point()+coord_fixed()

grid.arrange(g1,g2)

plots1sp <- df_nobg[df_nobg$CountOfSpecies==1,]

