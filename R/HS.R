# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/alpha-k-2013-3-8-9-59/"
PATH = "/opt/MATLAB_WORKSPACE/hs/dump/alpha-A-B-2013-3-6-23-16/"




DIR = "R-alpha-noA-noB-2013-3-6-20-8/"
DIR = "A-alpha-R-2013-3-10-10-25/"

DIR = "tau-sigma-by-alpha-R-noA-noB-truth_middle-2013-3-9-18-42/"

DIR = "A-B-alpha-R-tau-2013-3-10-0-29/"

DIR = "sigma_tau-2013-3-6-10-36/"



DIR = "vscaling-2013-3-10-22-3/"
DIR = "vscaling-tau-easy-to-converge-2013-3-11-9-31/"

DIR = "R-alpha-01-2013-3-10-21-4/"

DIR = "R-alpha-99-2013-3-10-21-37/"

DIR ="bottom-R-alpha-2013-3-11-10-33/"

DIR = "upper-R-alpha-2013-3-11-17-45/"

DIR = "limited_sigma_R/"

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
#for (n in names(clu[1:23])) {
#  clu[, n] <- as.factor(clu[, n])      
}
#clu$t <- as.factor(clu$t)
clu$fromtruth.avg.cut <- cut(clu$fromtruth.avg, seq(0,1,0.1))
clu$size.avg.cut <- cut(clu$size.avg, seq(0,100,5))
clu$count.cut <- cut(clu$count, seq(0,100,5))

#lagging

clu$fromtruth.lag <- lagg(clu$fromtruth.avg)
clu$truthdiff <- (clu$fromtruth.avg - clu$fromtruth.lag)
clu$ratediff <- clu$truthdiff / clu$fromtruth.avg
clu$Rjump <- clu$truthdiff > 0.02


# START

allPlots("R","sigma")

image(clu$R, clu$sigma, clu$fromtruth.avg)

heatmap2by2Detail("R","sigma")

### Print Convergence by R
# We need to have sigma and t as numeric, not as factors

AA <- clu[clu$t == 21,]

selected <- AA[AA$sigma == 0, ]
p <- ggplot(AA, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p
ggsave(filename=paste0(PATH,"convergence_by_r_alpha=.6_full.jpg"), plot=p)

selected <- AA[AA$sigma == 0 & AA$R < 0.6 & AA$R > 0.3,]
p <- ggplot(selected, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p
ggsave(filename=paste0(PATH,"convergence_by_r_alpha=.6_zoom.jpg"), plot=p)


selected <- AA[AA$sigma == 0.2 AA$R < 0.6 & AA$R > 0.3,]
p <- ggplot(selected, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p

### End



clu2 <- clu
clu2$init.vscaling <- as.numeric(clu2$init.vscaling)
clu2$tau <- as.numeric(clu2$tau)
clu2$t <- as.numeric(clu2$t)
clu2$alpha <- as.numeric(clu2$alpha)
clu2$R <- as.numeric(clu2$R)
clu2$run <- as.numeric(clu2$run)
clu2$simcount <- as.numeric(clu2$simcount)


cluB <- clu2[clu2$alpha == 0.1 | clu2$alpha == 0.2 | clu2$alpha == 0.05,]

plot.ts(cluB$t, cluB$fromtruth.avg)




cluB <- clu2[clu2$t == 31 & clu2$run==1 & clu2$alpha == 0.2,]

plot(cluB$R, cluB$fromtruth.avg)

heatmap2by2Detail("R","alpha", data=cluB)

# Spinning 3d Scatterplot
library(rgl)

plot3d(clu$init.vscaling, clu$t, clu$count, col="red", size=3) 
plot3d(clu$init.vscaling, clu$t, clu$fromtruth.avg, col="red", size=3) 


plot3d(clu2$R, clu2$alpha, clu2$count, col="red", size=3)

plot3d(clu2$R, clu2$alpha, clu2$fromtruth.avg, col="red", size=3) 


scatter3d(clu2$R, clu2$fromtruth.avg, clu2$alpha)

scatter3d(cluB$R, cluB$fromtruth.avg, cluB$alpha) 


scatter3d(clu2$init.vscaling, clu2$fromtruth.avg, clu2$tau) 

## JUST FOR A FEW POINTS
s3d <-scatterplot3d(clu2$init.vscaling, clu2$tau, clu2$fromtruth.avg, highlight.3d=TRUE,
type="h", main="3D Scatterplot")
fit <- lm(fromtruth.avg ~ tau+init.vscaling, data=clu2)
s3d$plane3d(fit)
##

                    
plotTS2by2("init.vscaling", "R")

boxplot2by2("init.vscaling")

boxplot(clu$fromtruth.avg ~ clu$init.vscaling)

boxplot(clu$count ~ clu$init.vscaling)

boxplot(clu$size.avg ~ clu$init.vscaling)

clu0 <- clu[clu$alpha == 0,]




# when B fucks it up
clu <- clu[clu$B == 0,]







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




R <- summarySE(clu, c("fromtruth.avg"), c("R","alpha", "sigma"))


last <- clu[clu$t == 21,]
R <- summarySE(last, c("fromtruth.avg"), c("R","alpha","sigma"))

A <- R[order(R$fromtruth.avg),]

A$ratio = as.numeric(A$alpha) / as.numeric(A$R)

p <- ggplot(A, aes(R, sigma))
p + geom_tile(aes(fill=fromtruth.avg)) + scale_fill_continuous(limits=c(0,max(clu$fromtruth.avg)), guide="legend", breaks= seq(0,1,0.05), low='lightblue',high='red')

obs=nrow(R)



plot.ts(1:nrow(R), R$truthdiff)

p <- ggplot(A[A$sigma == 0,], aes(R, fromtruth.avg, group=sigma, colour=sigma))
p + geom_line() + scale_fill_continuous(limits=c(0,max(clu$fromtruth.avg)), guide="legend", breaks= seq(0,1,0.05), low='lightblue',high='red')


+ geom_text(data = labeled.dat, aes(R, fromtruth.avg, label = R), size=4)


jumpR <- R[R$truthdiff > 0.02 & R$sigma == 0,]$R




plot.ts(1:nrow(A), A$R)


B <- A[sample(nrow(A)),]

#B$ratio = as.numeric(B$alpha) / as.numeric(B$R)

plot.ts(1:nrow(B), B$R)
