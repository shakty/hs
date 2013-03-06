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
setwd("/opt/MATLAB_WORKSPACE/hs/dump/sigma_tau_sigma_tau-2013-3-6-10-36")


params <- read.table('params2.csv', head=TRUE, sep=",")

clusters <- read.table('clusters.csv', head=TRUE, sep=",")

params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)



clusters$simname <- as.factor(clusters$simname)
clusters$simcount <- as.factor(clusters$simcount)
clusters$run <- as.factor(clusters$run)


clu <- merge(params, clusters, by=c("simname","simcount","run"))

p <- ggplot(clu, aes(t, count))

p + geom_smooth()


p.facets <- p + geom_smooth() + facet_grid(init.vscaling ~ ., margins = T);
p.facets <- p.facets + ggtitle("Number of clusters as a function of alpha")
p.facets

