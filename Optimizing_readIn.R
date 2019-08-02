vec <- readtext(file="~/Desktop/raw_4816")

library(sf)
library(caTools)
library(rgdal)
library(hyperSpec)
try <- st_read("~/Desktop/raw_4816.hdr")

atry <- read.csv(file="~/Desktop/raw_4816 copy.txt",header=F)
atry <- readBin(con="~/Desktop/raw_4816",what="character",n=100000)

ev <- caTools::read.ENVI("~/Desktop/raw_4816")
hs <- hyperSpec::read.ENVI("~/Desktop/raw_4816",headerfile = "~/Desktop/raw_4816.hdr")
init.dim <- dim(ev)
dim(ev) <- c(dim(ev)[1]*dim(ev)[2],dim(ev)[3])
#rep x 00000, 11111, 2222
ev$x <- rep((1:init.dim[1])-1, each = init.dim[2])

# rep y 0 1 2 3 4
ev$y <- rep((1:init.dim[2])-1,init.dim[1])
readLines("~/Desktop/raw_4816.hdr")
ev_df <- as.data.frame(ev)

gd <- readGDAL("~/Desktop/raw_4816")
gd_df <-as.data.frame(gd)
gd_df$y

