# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/alpha-k-2013-3-8-9-59/"
PATH = "/opt/MATLAB_WORKSPACE/hs/dump/alpha-A-B-2013-3-6-23-16/"




DIR = "R-alpha-noA-noB-2013-3-6-20-8/"
DIR = "A-alpha-R-2013-3-10-10-25/"
DIR = "sigma_tau-2013-3-6-10-36/"
DIR = "tau-sigma-by-alpha-R-noA-noB-truth_middle-2013-3-9-18-42/"

DIR = "A-B-alpha-R-tau-2013-3-10-0-29/"


DIR = "R-alpha-99-2013-3-10-21-37/"
DIR = "R-alpha-01-2013-3-10-21-4/"

DIR = "vscaling-2013-3-10-22-3/"

DUMPDIR = "/opt/MATLAB_WORKSPACE/hs/dump/"
PATH = paste0(DUMPDIR,DIR)
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
#clu$t <- as.factor(clu$t)
clu$fromtruth.avg.cut <- cut(clu$fromtruth.avg, seq(0,1,0.1))
clu$size.avg.cut <- cut(clu$size.avg, seq(0,100,5))
clu$count.cut <- cut(clu$count, seq(0,100,5))

# START

clu.orig <- clu

clu <- clu[clu$t == 31,]

heatmap2by2Detail("init.vscaling","t")

# Spinning 3d Scatterplot
library(rgl)

plot3d(clu$init.vscaling, clu$t, clu$count, col="red", size=3) 

plot3d(clu$init.vscaling, clu$t, clu$fromtruth.avg, col="red", size=3) 

clu2 <- clu
clu2$init.vscaling <- as.numeric(clu2$init.vscaling)

scatter3d(clu2$init.vscaling, clu2$fromtruth.avg, clu2$t) 


# 3D Scatterplot with Coloring and Vertical Lines
# and Regression Plane
library(scatterplot3d)

s3d <-scatterplot3d(clu$init.vscaling, clu$t, clu$fromtruth.avg, pch=16, highlight.3d=TRUE, type="h", main="3D Scatterplot")


fit <- lm(mpg ~ wt+disp)
s3d$plane3d(fit)
##

attach(mtcars)
s3d <-scatterplot3d(wt,disp,mpg, pch=16, highlight.3d=TRUE,
  type="h", main="3D Scatterplot")
fit <- lm(mpg ~ wt+disp)
s3d$plane3d(fit)
detach(mtcars)
                    
plotTS2by2("init.vscaling", "R")

boxplot2by2("init.vscaling")

boxplot(clu$fromtruth.avg ~ clu$init.vscaling)

boxplot(clu$count ~ clu$init.vscaling)

boxplot(clu$size.avg ~ clu$init.vscaling)

clu0 <- clu[clu$alpha == 0,]

plot.ts(clu0$R, clu0$fromtruth.avg)


# when B fucks it up
clu <- clu[clu$B == 0,]



allPlots("R","sigma")

## TODO CHANGE YSCALE for DIS in facets

# POINTS
#p.conv <- ggplot(clu[clu$tau==1,], aes(t, fromtruth.avg))
#p.conv + geom_point(aes(colour = sigma)) + geom_jitter(aes(colour = sigma)) + plotScaleDis


cluB <- clu[clu$B != 0,]

R <- summarySE(clu, c("fromtruth.avg"), c("R","sigma"))
sort(R,partial=c(4))

v1="R"
v2="alpha"

v3="R"
v4="alpha"

facetFormula <- as.formula(sprintf('%s~%s', v1, v2))
title <- paste0("Convergence levels in time by ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x="t", y="fromtruth.avg", group=v1, colour=v1))
  p <- p + geom_smooth()
  p <- p + facet_grid(facetFormula)
  p <- p + reducedXScale + yLabDis
  p <- p  + hs.makeggtitle(title, c(v1, v2))

p

  saveOrPlot(TRUE, p, "R1_alpha", IMGPATH)
#  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)


p <- ggplot(clu, aes(x=alpha, y=count, colour=B, group=B))
p <- p + geom_smooth()
p


cl <- clu[clu$B == 0 & clu$A == 0.1 & clu$tau == 9 & clu$R == 0,]
plot.ts(cl$alpha, cl$count)
