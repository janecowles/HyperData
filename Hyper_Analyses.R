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

sm_fin <- fread("~/Documents/HYPER-FinalSmoothedCorrectLocs_18June2019.csv")
ggplot(sm_fin[sm_fin$Ring==6,],aes(Longitude,Latitude,color=factor(monospecies)))+geom_point()
cv <- function(x){sd(x)/mean(x)}
cols <- colnames(sm_fin)[c(4,12:93,97)]
sm_fin_cv <- sm_fin[,lapply(.SD, cv),  .SDcols = 100:371,by=cols]
names(sm_fin_cv)

sm_fin_cv$meanCV <- rowMeans(sm_fin_cv[,85:356])
sm_fin_cv$meanCVcons <- rowMeans(sm_fin_cv[,98:317])
ggplot(sm_fin_cv,aes(CountOfSpecies,meanCV))+geom_point()+geom_point(stat="summary",col="red",size=3)
ggsave("~/Documents/sm_fin_CV1.pdf")

ggplot(sm_fin_cv,aes(CountOfSpecies,meanCVcons))+geom_point()+geom_point(stat="summary",col="red",size=3)
ggsave("~/Documents/sm_fin_cv2.pdf")

ggplot(sm_fin_cv,aes(meanCVcons,AbovegroundTotal.Biomass..g.m.2))+geom_point()
ggplot(sm_fin_cv,aes(meanCVcons,AbovegroundTotal.Biomass..g.m.2,color=C.and.N.treatment,group=C.and.N.treatment))+geom_point()+geom_smooth(method="lm")
sm_fin_cv[sm_fin_cv$AbovegroundTotal.Biomass..g.m.2<10&sm_fin_cv$meanCVcons>0.5,]

sm_fin_cv$RingFac <- as.factor(sm_fin_cv$Ring)
sm_fin_cv$PlotFac <- as.factor(sm_fin_cv$Plot)
mod <- lme(AbovegroundTotal.Biomass..g.m.2~meanCVcons,random=~1|Ring/Plot,data=sm_fin_cv)
summary(mod)
Anova(mod)

mod <- lme(meanCVcons~CountOfSpecies,random=~1|Ring/Plot,data=sm_fin_cv)
summary(mod)
Anova(mod)
cor(sm_fin_cv$meanCVcons,sm_fin_cv$CountOfSpecies)
cor(sm_fin_cv$meanCV,sm_fin_cv$CountOfSpecies)

ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=monospecies,color=monospecies))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/sm_fin_mean_longsp1.pdf")
ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("monospecies")
ggsave("~/Documents/sm_fin_mean_longsp2.pdf")


ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=Monogroup,color=Monogroup))+geom_point(stat="summary") +facet_grid(CO2.Treatment~Nitrogen.Treatment)
ggsave("~/Documents/sm_fin_mean_longgr1.pdf")

ggplot(sm_fin_mean_long[sm_fin_mean_long$CountOfSpecies==1&sm_fin_mean_long$wavelength<680,],aes(wavelength,reflectance,group=C.and.N.treatment,color=C.and.N.treatment))+geom_point(stat="summary") +facet_wrap("Monogroup")
ggsave("~/Documents/sm_fin_mean_longgr2.pdf")

