

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
dflong$nm <- as.numeric(substr(dflong$variable,3,11))
ggplot(dflong,aes(nm,value,group=Plot,color=Plot))+geom_line()


bc <- read.csv(paste0(LocalSource,"BioCON Master Harvest_190109_USE.csv"))
bc18 <- bc[bc$year==2018&bc$Season=="August",]
bc_df <- merge(bc18,df,by="Plot")

bc_dflong <- merge(bc18,dflong,by="Plot",all=T)

ggplot(bc_dflong,aes(nm,value,group=factor(CountOfSpecies),color=factor(CountOfSpecies)))+geom_line(stat="summary")+theme_classic()

ggplot(bc_dflong,aes(nm,value,group=Nitrogen.Treatment,color=Nitrogen.Treatment))+geom_line(stat="summary")+theme_classic()
ggplot(bc_dflong,aes(nm,value,group=Plot,color=factor(CountOfSpecies)))+geom_line(stat="summary")+theme_classic()+facet_grid(1~Nitrogen.Treatment)
ggplot(bc_dflong,aes(nm,value,group=factor(StartFrame),color=factor(StartFrame)))+geom_line(stat="summary")+theme_classic()+facet_grid(1~factor(Ring))

testCluster <- kmeans(bc_df[, 209:307], 2, nstart = 20)
fviz_cluster(testCluster,bc_df[, c(101,51)])
