# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/alpha-k-2013-3-8-9-59/"
PATH = "/opt/MATLAB_WORKSPACE/hs/dump/alpha-A-B-2013-3-6-23-16/"
PATH = "/opt/MATLAB_WORKSPACE/hs/dump/sigma_tau-2013-3-6-10-36/"


PATH = "/opt/MATLAB_WORKSPACE/hs/dump/R-alpha-noA-noB-2013-3-6-20-8/"

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/tau-sigma-by-alpha-R-noA-noB-truth_middle-2013-3-9-18-42/"

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/A-B-alpha-R-tau-2013-3-10-0-29/"

setwd(PATH)
IMGPATH <- paste0(PATH, "img/");
params <- read.table('params.csv', head=TRUE, sep=",")
clusters <- read.table('clusters.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)
clusters$simname <- as.factor(clusters$simname)
clusters$simcount <- as.factor(clusters$simcount)
clusters$run <- as.factor(clusters$run)
# Transforms params in factors
clu <- merge(params, clusters, by=c("simname","simcount","run"))
for (n in names(clu[1:23])) {
  clu[, n] <- as.factor(clu[, n])      
}
clu$t <- as.factor(clu$t)
clu$fromtruth.avg.cut <- cut(clu$fromtruth.avg, seq(0,1,0.1))
clu$size.avg.cut <- cut(clu$size.avg, seq(0,100,5))
clu$count.cut <- cut(clu$count, seq(0,100,5))

# START



allPlots("A","R")

## TODO CHANGE YSCALE for DIS in facets

# POINTS
#p.conv <- ggplot(clu[clu$tau==1,], aes(t, fromtruth.avg))
#p.conv + geom_point(aes(colour = sigma)) + geom_jitter(aes(colour = sigma)) + plotScaleDis

v1="sigma"
v2="tau"
v3="R"
v4="alpha"
facetFormula <- as.formula(sprintf('%s~%s~%s~%s', v2, v1, v3, v4))
title <- paste0("Convergence levels in time by ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x="t", y="fromtruth.avg", group=v1, colour=v1))
  p <- p + geom_smooth()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedXScale + yLabDis
  p <- p  + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(TRUE, p, "STE2", IMGPATH)
#  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)


p <- ggplot(clu, aes(x=alpha, y=count))
p <- p + geom_smooth(aes(group=R))
p


plot.ts(clu$alpha, clu$count)
