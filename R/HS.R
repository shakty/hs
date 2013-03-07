library(plyr)
library(ggplot2)
library(gridExtra)
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
for (n in names(clu[-1:-23])) {
  clu[, n] <- as.factor(clu[, n])      
}
clu$t <- as.factor(clu$t)



printParams <- function(varlist) {
  out <- '';
  for (n in names(params[c(-1:-9,-19:-length(params))])) {
    if (!(n %in% varlist)) {
      value <- params[1, n]
      if (class(params[,n]) == "integer") {
        out <- sprintf('%s %s: %d, ', out, n, value)
      }
      else {
         out <- sprintf('%s %s: %2.3f, ', out, n, value)
       }
    }
  }
  
  return(substr(out, 2, nchar(out)-2))
}

# START


paramAnnotate <- annotation_custom(grob = paramString,  xmin = 8, xmax = 8, ymin = -3, ymax = -3)
plotMargins <- theme(plot.margin = unit(c(1,3,8,1), "lines"))


## AXIS SCALES
plotScaleDis <- scale_y_discrete(breaks=seq(0, 1, 0.05))
plotScaleSize <- scale_y_discrete(breaks=seq(0, 100, 5))
plotScaleCount <- scale_y_discrete(breaks=seq(0, 100, 5))

reducedXScale <- scale_x_discrete(breaks=seq(5, 25, 10))

reducedYScaleCount <- scale_y_discrete(breaks=seq(0, 100, 10))
reducedYScaleDis <- scale_y_discrete(breaks=seq(0, 1, 0.2))
reducedYScaleSize <- scale_y_discrete(breaks=seq(0, 100, 10))

hs.makeggtitle <- function(text, vars) {
  paramString <- printParams(vars)
  return(ggtitle(paste0(text,paramString)))
}

# TIME EVOLUTION

#count
p <- ggplot(clu, aes(t, count, group=sigma, colour=sigma))
p <- p + geom_smooth() + ylab('number of clusters') + plotScaleCount
p + hs.makeggtitle("Cluster as function of noise\n\n", c("sigma", "tau"))

#size
p <- ggplot(clu, aes(t, size.avg, group=sigma, colour=sigma))
p + geom_smooth() + ylab('cluster size') + plotScaleSize +
  hs.makeggtitle("Avg cluster size as function of noise\n\n", c("sigma", "tau"))

#from truth
p <- ggplot(clu, aes(t, fromtruth.avg, group=sigma, colour=sigma))
p + geom_smooth() + ylab('distance from truth') + plotScaleDis + hs.makeggtitle("Distance from truth as function of noise\n\n", c("sigma", "tau"))

# BOXPLOTS

#count
p <- ggplot(clu, aes(sigma, count, group=sigma, colour=sigma))
p + geom_boxplot() + ylab('number of clusters') + plotScaleCount +
  hs.makeggtitle("Cluster as function of noise\n\n", c("sigma", "tau"))

#size
p <- ggplot(clu, aes(sigma, size.avg, group=sigma, colour=sigma))
p + geom_boxplot() + ylab('cluster size') + plotScaleSize +
  hs.makeggtitle("Avg cluster size as function of noise\n\n", c("sigma", "tau"))

#from truth
p <- ggplot(clu, aes(sigma, fromtruth.avg, group=sigma, colour=sigma))
p + geom_boxplot() + ylab('distance from truth') + plotScaleDis + 
  hs.makeggtitle("Distance from truth as function of noise\n\n", c("sigma", "tau"))


# POINTS

p.conv <- ggplot(clu[clu$tau==1,], aes(t, fromtruth.avg))
p.conv + geom_point(aes(colour = sigma)) + geom_jitter(aes(colour = sigma)) + plotScaleDis

## FACETS


#count ts
p <- ggplot(clu, aes(t, count, group=sigma, colour=sigma))
p.facets <- p + geom_smooth() + facet_grid(sigma ~ tau, margins = T) + plotScaleCount + reducedXScale
p.facets <- p.facets + ggtitle("Number of clusters as a function of sigma and tau")
p.facets 

#count boxplot
p <- ggplot(clu, aes(t, count, group=sigma, colour=sigma))
p.facets <- p + geom_boxplot() + facet_grid(sigma ~ tau, margins = T) + plotScaleCount + reducedXScale
p.facets <- p.facets + ggtitle("Number of clusters as a function of sigma and tau")
p.facets

#size ts
p <- ggplot(clu, aes(t, size.avg, group=sigma, colour=sigma))
p.facets <- p + geom_smooth() + facet_grid(sigma ~ tau, margins = T) + plotScaleCount + reducedXScale
p.facets <- p.facets + ggtitle("Number of clusters as a function of sigma and tau")
p.facets 

#size boxplot
p <- ggplot(clu, aes(t, size.avg, group=sigma, colour=sigma))
p.facets <- p + geom_boxplot() + facet_grid(sigma ~ tau, margins = T) + plotScaleCount + reducedXScale
p.facets <- p.facets + ggtitle("Number of clusters as a function of sigma and tau")
p.facets

#fromtruth ts
p <- ggplot(clu, aes(t, fromtruth.avg, group=sigma, colour=sigma))
p + geom_smooth() + facet_grid(tau ~ sigma, margins = T) +
  reducedYScaleDis + reducedXScale +
  ylab('distance from truth') +
  hs.makeggtitle("Distance from truth as function of noise and tau\n\n", c("sigma", "tau"))


#fromtruth boxplot
p <- ggplot(clu, aes(tau, fromtruth.avg, group=sigma, colour=sigma))
p + geom_boxplot() + facet_grid(tau ~ sigma, margins = T) +
  reducedYScaleDis +
  ylab('distance from truth') +
  hs.makeggtitle("Distance from truth as function of noise and tau\n\n", c("sigma", "tau"))
  


