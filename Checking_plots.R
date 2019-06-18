library(sf)
library(tidyverse)
check2816 <- st_read("/Volumes/HyperDrive/Final FULL/Final2816full.shp")
plot(check2816$Lat2~check2816$Lon2,col=check2816$nm540_751)
df <- as.data.frame(check2816)
ggplot(df,aes(Lon2,Lat2,color=nm540_751))+geom_point()

dfs <- df[df$Lat2>mean(df$Lat2),]
dfs <- dfs[dfs$Lon2<mean(df$Lon2),]
hist(dfs$Lat2)
ggplot(dfs,aes(Lon2,Lat2,color=nm540_751))+geom_point()
