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

setwd("/opt/MATLAB_WORKSPACE/hs/dump/sigma_tau-2013-3-6-10-36")

params <- read.table('params.csv', head=TRUE, sep=",")
clusters <- read.table('clusters.csv', head=TRUE, sep=",")

params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)



clusters$simname <- as.factor(clusters$simname)
clusters$simcount <- as.factor(clusters$simcount)
clusters$run <- as.factor(clusters$run)


clu <- merge(params, clusters, by=c("simname","simcount","run"))

p <- ggplot(clu, aes(t, count, group=sigma, colour=sigma))
p + geom_smooth()

p <- ggplot(clu, aes(sigma, count, group=sigma, colour=sigma))
p + geom_boxplot()

p <- ggplot(clu, aes(sigma, fromtruth.avg, group=sigma, colour=sigma))
p + geom_boxplot()

p <- ggplot(clu, aes(sigma, fromtruth.sd, group=sigma, colour=sigma))
p + geom_boxplot()

p <- ggplot(clu, aes(sigma, size.avg, group=sigma, colour=sigma))
p + geom_boxplot()

p <- ggplot(clu, aes(sigma, size.sd, group=sigma, colour=sigma))
p + geom_boxplot()

clu$sigmaf <- as.factor(clu$sigma)

p <- ggplot(clu, aes(t, count, group=sigmaf, colour=sigmaf))
p + geom_smooth()

clu$tauf <- as.factor(clu$tau)
p <- ggplot(clu, aes(t, fromtruth.avg, group=tauf, colour=tauf))
p + geom_smooth()

clu2 <- clu

clu2.sigma <- as.factor(clu2$sigma)

p.conv <- ggplot(clu2, aes(t, fromtruth.avg))
p.conv + geom_point(aes(colour = sigma)) + geom_jitter(aes(colour = sigma))

p <- ggplot(clu, aes(t, count, group=sigmaf, colour=sigmaf))
p.facets <- p + geom_smooth() + facet_grid(sigma ~ tau, margins = T);
p.facets <- p.facets + ggtitle("Number of clusters as a function of sigma and tau")
p.facets

