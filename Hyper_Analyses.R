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

library(FD)


###data sets I may need
#with BG
df_keep <- fread("~/Documents/FINAL_NoDups_WITHBG24June2019.csv")

#before smoothing
ExpBGremove <- fread("~/Documents/RF_BGRemoval24June2019.csv")

#pixels removed for bareground metric
pixels_removed <- fread("~/Documents/PlotLevelPixelsRemoved_24June2019.csv")
#smoothed!
sm_fin <- fread("~/Documents/HYPER-FinalSmoothedCorrectLocs_23June2019.csv")
ggplot(sm_fin[sm_fin$Ring==2,],aes(Longitude,Latitude,color=factor(monospecies)))+geom_point()
ggplot(ExpBGremove[ExpBGremove$Ring==2,],aes(Lon2,Lat2,color=factor(Use.for.Main)))+geom_point()

ggplot(sm_fin[sm_fin$Ring==2,],aes(Longitude,Latitude,color=factor(CountOfSpecies)))+geom_point()+coord_fixed()+scale_color_manual(values=c("darkorange","deepskyblue","forestgreen","darkslategrey"))+theme_minimal()
ggsave("~/Documents/Ring2_div18June.pdf")


############## plot level stuff!!
#figure out convex hull volume
i=1
for(i in unique(sm_fin$Plot)){
  tmp <- sm_fin[sm_fin$Plot==i,]
    testing <- dbFD(tmp[,c(200:235)])
}








#### for CV, need to figure out w, which is the correction factor for spectral diversity. This is in pixels_removed. Calculate 
#soil adjusted spectral diversity = SD X w
# w = (total pixels - total soil pixels) / total pixels
# w = non bareground pixels/ total pixels
pixels_removed$w <- pixels_removed$numberofpixels/pixels_removed$preremovalnumberofpixels

cv <- function(x){sd(x)/mean(x)}
cols <- colnames(sm_fin)[c(4,12:93,97)]
sm_fin_cv <- sm_fin[,lapply(.SD, cv),  .SDcols = 100:371,by=cols]
names(sm_fin_cv)

sm_fin_cv$meanCV <- rowMeans(sm_fin_cv[,85:356])
sm_fin_cv$meanCVcons <- rowMeans(sm_fin_cv[,100:322])
sm_fin_cvP <- merge(sm_fin_cv,pixels_removed[,c("Plot","w","preremovalnumberofpixels", "numberofpixels")],by="Plot")
sm_fin_cvP$SoilCorCV <- sm_fin_cvP$meanCVcons*sm_fin_cvP$w

# sm_fin_cvP<-sm_fin_cvP[!is.na(sm_fin_cvP$meanCV)]



names(sm_fin_cvP)
ggplot(sm_fin_cvP,aes(CountOfSpecies,SoilCorCV))+geom_point()+geom_point(stat="summary",col="red",size=3)
ggsave("~/Documents/sm_fin_CV121July.pdf")

ggplot(sm_fin_cvP,aes(SoilCorCV, AbovegroundTotal.Biomass..g.m.2, color=factor(CountOfSpecies)))+geom_point()
cor(sm_fin_cvP$SoilCorCV,sm_fin_cvP$AbovegroundTotal.Biomass..g.m.2)

ggplot(sm_fin_cvP,aes(SoilCorCV,X..cover.Total.Planted.Species, color=factor(CountOfSpecies)))+geom_point()
cor(sm_fin_cvP$SoilCorCV,sm_fin_cvP$X..cover.Total.Planted.Species)

ggplot(sm_fin_cvP,aes(meanCV, AbovegroundTotal.Biomass..g.m.2, color=factor(CountOfSpecies)))+geom_point()
cor(sm_fin_cvP$meanCV,sm_fin_cvP$AbovegroundTotal.Biomass..g.m.2)
ggplot(sm_fin_cvP,aes(meanCV,X..cover.Total.Planted.Species, color=factor(CountOfSpecies)))+geom_point()
cor(sm_fin_cvP$meanCV,sm_fin_cvP$X..cover.Total.Planted.Species)

ggplot(sm_fin_cvP,aes(meanCVcons, AbovegroundTotal.Biomass..g.m.2, color=factor(CountOfSpecies)))+geom_point()+geom_smooth(method="lm")
cor(sm_fin_cvP$meanCVcons,sm_fin_cvP$AbovegroundTotal.Biomass..g.m.2)
ggplot(sm_fin_cvP,aes(meanCVcons,X..cover.Total.Planted.Species, color=factor(CountOfSpecies)))+geom_point()+geom_smooth(method="lm")
cor(sm_fin_cvP$meanCVcons,sm_fin_cvP$X..cover.Total.Planted.Species)

ggplot(sm_fin_cvP,aes(meanCVcons, AbovegroundTotal.Biomass..g.m.2, color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm")
ggplot(sm_fin_cvP,aes(meanCVcons, X..cover.Total.Planted.Species, color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm")



# # ggplot(sm_fin_cv,aes(CountOfSpecies,meanCVcons))+geom_point()+geom_point(stat="summary",col="red",size=3)
# # ggsave("~/Documents/sm_fin_cv2.pdf")
# 
# ggplot(sm_fin_cv,aes(meanCVcons,AbovegroundTotal.Biomass..g.m.2))+geom_point()
# ggplot(sm_fin_cv,aes(meanCVcons,AbovegroundTotal.Biomass..g.m.2,color=C.and.N.treatment,group=C.and.N.treatment))+geom_point()+geom_smooth(method="lm")
# sm_fin_cv[sm_fin_cv$AbovegroundTotal.Biomass..g.m.2<10&sm_fin_cv$meanCVcons>0.5,]
# 
sm_fin_cvP$RingFac <- as.factor(sm_fin_cvP$Ring)
sm_fin_cvP$PlotFac <- as.factor(sm_fin_cvP$Plot)

mod <- lme(AbovegroundTotal.Biomass..g.m.2~meanCVcons,random=~1|Ring/Plot,data=sm_fin_cvP,na.action=na.exclude)
summary(mod)
Anova(mod)


mod <- lme(AbovegroundTotal.Biomass..g.m.2~SoilCorCV,random=~1|Ring/Plot,data=sm_fin_cvP)
summary(mod)
Anova(mod)

mod <- lme(meanCVcons~CountOfSpecies,random=~1|Ring/Plot,data=sm_fin_cvP)
summary(mod)
Anova(mod)

mod <- lme(SoilCorCV~CountOfSpecies,random=~1|Ring/Plot,data=sm_fin_cvP)
summary(mod)
Anova(mod)

cor(sm_fin_cvP$SoilCorCV,sm_fin_cvP$CountOfSpecies)
cor(sm_fin_cvP$meanCVcons,sm_fin_cvP$CountOfSpecies)





#############means/spectra
cols <- colnames(sm_fin)[c(4,12:93,97)]
sm_fin_means <- sm_fin[,lapply(.SD, mean),  .SDcols = 100:371,by=cols]
names(sm_fin_means)
sm_fin_means$NDVI <- (sm_fin_means$`800.614`-sm_fin_means$`660.688`)/(sm_fin_means$`800.614`+sm_fin_means$`660.688`)


#merge with cvP
plotdf <- merge(sm_fin_means,sm_fin_cvP[,c("Plot","meanCVcons","meanCV","SoilCorCV","w","preremovalnumberofpixels", "numberofpixels")],by="Plot")
plotdf<-plotdf[plotdf$numberofpixels>10,]

ggplot(plotdf,aes(SoilCorCV,NDVI,color=factor(CountOfSpecies)))+geom_point()+geom_smooth(method="lm")
ggplot(plotdf,aes(meanCV,NDVI,color=factor(CountOfSpecies)))+geom_point()+geom_smooth(method="lm")
ggplot(plotdf,aes(meanCVcons,NDVI,color=factor(CountOfSpecies)))+geom_point()+geom_smooth(method="lm")

plotdf$Ring <- as.factor(plotdf$Ring)

cor(plotdf$meanCVcons,plotdf$X..cover.Total.Planted.Species)
cor(plotdf$meanCVcons,plotdf$AGB)
cor(plotdf$meanCVcons,plotdf$NDVI)

###inppt
ggplot(plotdf,aes(factor(CountOfSpecies),meanCVcons))+geom_violin()+geom_point(aes(color=factor(CountOfSpecies)))+labs(x="Planted Species Richness",y="Coefficient of Variance (430-925 nm)")+theme_minimal()+theme(legend.position = "none")
###inppt
ggplot(plotdf,aes(Species.count.from..Cover,meanCVcons))+geom_point()+geom_smooth(method="lm",size=1)+theme_minimal()+labs(x="Observed Species Richness",y="Coefficient of Variance (430-925 nm)")


ggplot(plotdf,aes(factor(CountOfSpecies),meanCVcons))+geom_violin()+geom_point(aes(color=numberofpixels))+labs(x="Planted Species Richness",y="Coefficient of Variance (430-925 nm)")

ggplot(plotdf,aes(CountOfSpecies,meanCVcons,color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm",size=2)

ggplot(plotdf,aes(Species.count.from..Cover,meanCVcons,color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm",size=2)

# ggplot(plotdf,aes(Species.count.from.clip.strip,meanCVcons,color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm",size=2)

ggplot(plotdf,aes(meanCVcons,AGB,color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm",size=2)
ggplot(plotdf,aes(meanCVcons,X..cover.Total.Planted.Species,color=factor(C.and.N.treatment)))+geom_point()+geom_smooth(method="lm",size=2)

ggplot(plotdf,aes(meanCVcons,AGB))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$meanCVcons,plotdf$AGB)

ggplot(plotdf,aes(meanCVcons,X..cover.Total.Planted.Species,size=numberofpixels))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$meanCVcons,plotdf$X..cover.Total.Planted.Species)

ggplot(plotdf,aes(SoilCorCV,X..cover.Total.Planted.Species,size=numberofpixels))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$SoilCorCV,plotdf$X..cover.Total.Planted.Species)

ggplot(plotdf,aes(meanCVcons,AGB,size=numberofpixels))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$meanCVcons,plotdf$AGB)

ggplot(plotdf,aes(SoilCorCV,AGB,size=numberofpixels))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$SoilCorCV,plotdf$AGB)

ggplot(plotdf,aes(meanCVcons,AGB,color=C.and.N.treatment))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$meanCVcons,plotdf$AGB)

ggplot(plotdf,aes(SoilCorCV,AGB,color=C.and.N.treatment))+geom_point()+geom_smooth(method="lm",size=2)
cor(plotdf$SoilCorCV,plotdf$AGB)


mod <- lme(AGB~meanCVcons*Nitrogen.Treatment*CO2.Treatment,  random=~1|Ring,data=plotdf)
summary(mod)
Anova(mod)
mod <- lme(X..cover.Total.Planted.Species~meanCVcons*Nitrogen.Treatment*CO2.Treatment,  random=~1|Ring,data=plotdf)
summary(mod)
Anova(mod)


mod <- lme(Species.count.from..Cover~meanCVcons,  random=~1|Ring,data=plotdf)
summary(mod)
Anova(mod)




sm_fin_cv_long <- melt(sm_fin_cvP, id.vars = c(1:84,359),measure.vars = c(85:356),variable.name = "wavelength",value.name = "cv")



sm_fin_mean_long <- melt(sm_fin_means, id.vars = c(1:84),measure.vars = c(85:356),variable.name = "wavelength",value.name = "reflectance")

# sm_fin_mean_long$wavelength<-as.numeric(as.character(sm_fin_mean_long$wavelength))

meancv_long <- merge(sm_fin_cv_long,sm_fin_mean_long[,c("Plot","wavelength","reflectance")],by=c("Plot","wavelength"))
meancv_long<-meancv_long[meancv_long$w>0.02,]
meancv_long$wavelength<-as.numeric(as.character(meancv_long$wavelength))
ggplot(meancv_long[meancv_long$wavelength>430&meancv_long$wavelength<925,],aes(wavelength,reflectance,color=monospecies))+geom_point()
ggplot(meancv_long[meancv_long$wavelength>430&meancv_long$wavelength<925,],aes(wavelength,reflectance,color=w))+geom_point()

meancv_long$Plot[meancv_long$wavelength>595&meancv_long$wavelength<602&meancv_long$reflectance>0.19]

ggplot(meancv_long[meancv_long$wavelength>430&meancv_long$wavelength<925,],aes(wavelength,cv,color=CountOfSpecies))+geom_point(stat="summary")
ggplot(meancv_long[meancv_long$wavelength>400&meancv_long$wavelength<925,],aes(wavelength,reflectance,color=CountOfSpecies))+geom_point(stat="summary")


ggplot(meancv_long[meancv_long$wavelength>430&meancv_long$wavelength<925,],aes(wavelength,cv,color=reflectance))+geom_point()
ggplot(meancv_long[meancv_long$wavelength>430&meancv_long$wavelength<925,],aes(wavelength,cv,color=w))+geom_point()

meancv_long$Plot[meancv_long$w<0.002&meancv_long$wavelength>430&meancv_long$wavelength<925]
####

ggplot(meancv_long[meancv_long$CountOfSpecies==1&meancv_long$wavelength<680,],aes(wavelength,reflectance,group=monospecies,color=monospecies))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/meancv_longsp1.pdf")
ggplot(meancv_long[meancv_long$CountOfSpecies==1&meancv_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_line(stat="summary") +facet_wrap("monospecies",nrow=5)+theme_minimal()+theme(legend.position = "bottom")
ggsave("~/Documents/meancv_longsp2.pdf")

ggplot(meancv_long[meancv_long$CountOfSpecies==1,],aes(wavelength,reflectance,group=Monogroup,color=Monogroup))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/meancv_longsp1.pdf")
ggplot(meancv_long[meancv_long$CountOfSpecies==1&meancv_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("monospecies")
ggsave("~/Documents/meancv_longsp2.pdf")









ggplot(sm_fin_mean_long,aes(wavelength,reflectance,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_point(stat="summary")+facet_grid(CO2.Treatment~Nitrogen.Treatment)+scale_color_manual(values=c("red","purple","orange","blue"))
ggsave("~/Documents/sm_fin_mean_long1.pdf")
ggplot(sm_fin_mean_long,aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary")+facet_wrap("CountOfSpecies")+scale_color_manual(values=c("red","purple","orange","blue"))+theme_minimal()
ggsave("~/Documents/sm_fin_mean_long2.pdf")


ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=monospecies,color=monospecies))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/sm_fin_mean_longsp1.pdf")
ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("monospecies")
ggsave("~/Documents/sm_fin_mean_longsp2.pdf")


ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=Monogroup,color=Monogroup))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/sm_fin_mean_longgr1.pdf")

ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("Monogroup")
ggsave("~/Documents/sm_fin_mean_longgr2.pdf")





####################### ok -- take sm_fin and subset to monocultures for the training sets
sm_fin2 <- sm_fin
sm_fin2$GrassorNot <- "Not"
sm_fin2$GrassorNot[sm_fin2$Monogroup%in%c("C-3","C-4")]<- "Grass"
sm_fin2$GrassorNot<-as.factor(sm_fin2$GrassorNot)
colnames(sm_fin2)[c(100:371)]<-paste0("nm",colnames(sm_fin2)[c(100:371)])
sm_train<- sm_fin2[sm_fin2$CountOfSpecies==1&sm_fin2$C.and.N.treatment=="AC, AN",]
sm_train$Monogroup<-factor(sm_train$Monogroup)
sm_train$monospecies<-factor(sm_train$monospecies)
sm_train$Nitrogen.Treatment <- factor(sm_train$Nitrogen.Treatment)
sm_trainmini <- sm_train[5000:12000,]
sm_trainmini$Monogroup<-factor(sm_trainmini$Monogroup)
sm_trainmini$monospecies<-factor(sm_trainmini$monospecies)

wavelengthcols <- paste(colnames(sm_train[,c(100:371)]),collapse="+")
formula1 <- as.formula(paste("GrassorNot~",wavelengthcols,sep=""))
system.time(rf1 <- randomForest(formula1,data = sm_train,ntree=1000, mtry=200))
rf1
beep()
system.time(sm_fin2$preds <- predict(rf1,sm_fin2))
table(sm_fin2$GrassorNot,sm_fin2$preds)
ggplot(sm_fin2[sm_fin2$Ring==1,],aes(Longitude,Latitude,color=preds))+facet_grid(1~Nitrogen.Treatment)+geom_point()

########
dffin <- ExpBGremove[ExpBGremove$Use.for.Main=="Y",]
dffin2 <- dffin
fin_train<- dffin2[dffin2$CountOfSpecies==1,]
fin_train$Monogroup<-factor(fin_train$Monogroup)
fin_train$monospecies<-factor(fin_train$monospecies)
fin_trainmini <- fin_train[1:15000,]
fin_trainmini$Monogroup<-factor(fin_trainmini$Monogroup)
fin_trainmini$monospecies<-factor(fin_trainmini$monospecies)

wavelengthcols <- paste(colnames(fin_train[,c(94:365)]),collapse="+")
formula1 <- as.formula(paste("Monogroup~",wavelengthcols,sep=""))
system.time(rf1 <- randomForest(formula1,data = fin_trainmini))
system.time(dffin2$preds <- predict(rf1,df_keep))

ggplot(sm_fin_cvP[sm_fin_cvP$CountOfSpecies==1,],aes(Monogroup,w))+geom_point()
ggplot(sm_fin_cvP[sm_fin_cvP$CountOfSpecies==1,],aes(w,X..cover.Total.Planted.Species))+geom_point()






