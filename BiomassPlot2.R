str(out_df37901)
names(out_df37901)

dfall <- NULL
readin_fun <- function(filename){

dftmp <- read.csv(paste0("~/",filename))
dfall <- rbind(dfall,dftmp)
return(dfall)
}

outfi <- list.files("~/",pattern="out_df")

df<-rbindlist(lapply(outfi,readin_fun))
df<-df[!is.na(df$nm398.604)]
setDT(df)
df[,NumReps:=length(unique(StartFrame)),by=.(Plot)]

dflong <- melt(df,id.vars=c("Lat2","Lon2","Plot","StartFrame"),measure.vars = c(2:273))
ggplot(dflong,aes(variable,value,group=Plot,color=Plot))+geom_line()

# ndvi (800-670) / (800+670)
out_df37901$NDVI_37901 <- (out_df37901$nm800.614 - out_df37901$nm669.572)/(out_df37901$nm800.614 + out_df37901$nm669.572)
out_df41565$NDVI_41565 <- (out_df41565$nm800.614 - out_df41565$nm669.572)/(out_df41565$nm800.614 + out_df41565$nm669.572)


out_df37901$CV_ALL_37901 <- rowMeans(out_df37901[,282:553])
out_df41565$CV_ALL_41565 <- rowMeans(out_df41565[,282:553])

bc <- read.csv(paste0(LocalSource,"BioCON Master Harvest_190104 for DB.csv"))
bc18 <- bc[bc$year==2018&bc$Season=="August",]

bc18a <- merge(bc18,out_df37901[,c("Plot","NDVI_37901","CV_ALL_37901")],by="Plot")
df <- merge(bc18a,out_df41565[,c("Plot","NDVI_41565","CV_ALL_41565")],by="Plot")
df$NDVI_DIFF <- df$NDVI_37901-df$NDVI_41565
plot(NDVI_37901 ~ AbovegroundTotal.Biomass..g.m.2,df,ylab="NDVI",xlab="Aboveground Biomass (g/m2)",pch=16)
points(NDVI_41565 ~ AbovegroundTotal.Biomass..g.m.2,df,col="red",pch=16)


plot(NDVI_37901 ~ CountOfSpecies,df,ylab="NDVI",xlab="Planted Species Richness",pch=16)
points(NDVI_41565 ~ CountOfSpecies,df,col="red",pch=16)




plot(CV_ALL_37901 ~ CountOfSpecies,df,ylab="NDVI VAR",xlab="Planted Species Richness",pch=16)
points(CV_ALL_41565 ~ CountOfSpecies,df,col="red",pch=16)


plot(NDVI_37901 ~ CV_ALL_37901,df,ylab="NDVI",xlab="CV (mean across wavelengths)",pch=16)
points(NDVI_41565 ~ CV_ALL_41565,df,col="red",pch=16)


plot(NDVI_37901 ~ Plot,df,ylab="NDVI",xlab="PLOT",pch=16,xlim=c(185,244))
points(NDVI_41565 ~ Plot,df,col="red",pch=16)

plot(NDVI_DIFF~Plot,df)
