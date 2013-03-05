library(plyr)
library(ggplot2)
library(reshape)
library(MASS)
library(car)

hs.setwd <- function(DIR, session){
  DATADIR = sprintf("%s%s/csv/", DIR, session)
  setwd(DATADIR)
  getwd()
}

hs.source <- function(sourcefile) {
 FULLSOURCE = sprintf("%s/%s", scriptdir, sourcefile)
 source(FULLSOURCE)
}

setwd("/opt/MATLAB_WORKSPACE/hs/dump/few_big_groups-DIM-vs-ALPHA/")

params <- read.table('params.csv', head=TRUE, sep=",")

clusters <- read.table('clusters.csv', head=TRUE, sep=",")

clu <- merge(params, clusters, by=c("sim","run"))

p <- ggplot(clu, aes(t, count))

p + geom_smooth()


p.facets <- p + geom_smooth() + facet_grid(init.vscaling ~ ., margins = T);
p.facets <- p.facets + ggtitle("Number of clusters as a function of alpha")
p.facets

