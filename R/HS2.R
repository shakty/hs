# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")


DIR = "attrExpo_nv_rndseq_tm_Rleft"

DUMPDIR = "/opt/MATLAB_WORKSPACE/hs/dump/"
PATH = paste0(DUMPDIR,DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");
params <- read.table('params.csv', head=TRUE, sep=",")

clusters <- read.table('clusters_macro.csv', head=TRUE, sep=",")
clusters.micro <- read.table('clusters_micro.csv', head=TRUE, sep=",")

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
# clu$t <- as.factor(clu$t)
clu$fromtruth.avg.cut <- cut(clu$fromtruth.avg, seq(0,1,0.1))
clu$size.avg.cut <- cut(clu$size.avg, seq(0,100,5))
clu$count.cut <- cut(clu$count, seq(0,100,5))

# lagging

clu$fromtruth.lag <- lagg(clu$fromtruth.avg)
clu$truthdiff <- (clu$fromtruth.avg - clu$fromtruth.lag)
clu$ratediff <- clu$truthdiff / clu$fromtruth.avg
clu$Rjump <- clu$truthdiff > 0.02

# SimCount

a <- params[params$alpha == 0.4 & params$sigma == 0 & params$R == 0.07,]
as.numeric(a$simcount)       

a <- clu[clu$fromtruth.avg > 0.6,]

max.dist.from.truth <- which(clu$fromtruth.avg == max(clu$fromtruth.avg))

clu[max.dist.from.truth,]


# Cluster Micro

CUTOFF <- 1
START.AFTER.ROUND <- 50
clusters.micro.sub <- clusters.micro[clusters.micro$size > CUTOFF & clusters$t > START.AFTER.ROUND, ]

p <- ggplot(clusters.micro.sub, aes(size, speed))
p <- p + geom_point(aes(colour = as.factor(simcount)))
p

p <- ggplot(clusters.micro.sub, aes(size, fromtruth))
p <- p + geom_point()
p

p <- ggplot(clusters.micro.sub, aes(size, move))
p <- p + geom_point()
p

nrow(a)

# Cluster Macro in time

cl <- clusters

title = "Evolution of cluster count and size"
p.count <- ggplot(cl, aes(t))
p.count <- p.count + geom_jitter(aes(y = size.avg), alpha=.2)
p.count <- p.count + geom_smooth(aes(y = count, colour="count"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.avg, colour="size"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.sd, colour="std. size"), size=2)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents")
p.count

title = "Average and cumulative space exploration"
p.explo <- ggplot(cl, aes(t))
p.explo <- p.explo + geom_jitter(aes(y = coverage), alpha=.2)
p.explo <- p.explo + geom_smooth(aes(y = coverage, colour="avg"), size=2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.cum, colour="cum"), size=2)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
p.explo

title = "Mean and std. agents speed"
p.speed <- ggplot(cl, aes(t))
p.speed <- p.speed + geom_jitter(aes(y = speed.avg), alpha=.2)
p.speed <- p.speed + geom_smooth(aes(y = speed.avg, colour="avg"), size=2)
p.speed <- p.speed + geom_smooth(aes(y = speed.sd, colour="sd"), size=2)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
p.speed

title = "Mean and std. agents movements"
p.move <- ggplot(cl, aes(t))
p.move <- p.move + geom_jitter(aes(y = move.avg), alpha=.2)
p.move <- p.move + geom_smooth(aes(y = move.avg, colour="avg"), size=2)
p.move <- p.move + geom_smooth(aes(y = move.sd, colour="sd"), size=2)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Displacement")
p.move

title = "Mean and std. agents distance from truth"
p.truth <- ggplot(cl, aes(t))
p.truth <- p.truth + geom_jitter(aes(y = fromtruth.avg), alpha=.2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.avg, colour="avg"), size=2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.sd, colour="sd"), size=2)
p.truth <- p.truth + ggtitle(title) + xlab("Rounds") + ylab("Distance")
p.truth

multiplot(p.explo, p.count, p.speed, p.move, p.truth, cols=3)

# START


v1 <- "R"
v2 <- "alpha"
v3 <- "sigma"
#data <- clu
data <- clu[clu$t == 21,]
paramsData <- params
heatmapFacets(v1,v2,v3, data)


OLDPATH = IMGPATH
for (S in unique(params$sigma)) {
  curDir <-  paste0("sigma_",S,"/")
  dir.create(file.path(OLDPATH, curDir), showWarnings = FALSE)
  IMGPATH <- paste0(OLDPATH, curDir)
  heatmap2by2Detail(v1,v2, data = clu[clu$sigma == S,], paramsData = params[params$sigma == S,])
}
IMGPATH = OLDPATH


#image(clu$R, clu$sigma, clu$fromtruth.avg)

# ALL PLOTS
#allPlots("R","sigma")

### Print Convergence by R
# We need to have sigma and t as numeric, not as factors

AA <- clu[clu$t == 21,]

selected <- AA[AA$sigma == 0, ]
p <- ggplot(selected, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p

ggsave(filename=paste0(PATH,"convergence_by_r_alpha=.6_full.jpg"), plot=p)

selected <- AA[AA$sigma == 0 & AA$R < 0.6 & AA$R > 0.3,]
p <- ggplot(selected, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p

ggsave(filename=paste0(PATH,"convergence_by_r_alpha=.6_zoom.jpg"), plot=p)




selected <- AA #[AA$sigma == 0.1,]
selected$sigma <- as.factor(selected$sigma)
p <- ggplot(selected, aes(R, fromtruth.avg, group=sigma, colour=sigma))
p <- p + geom_line() + facet_grid(sigma~.)
p

selected$sigma <- as.factor(selected$sigma)
p <- ggplot(selected, aes(R, fromtruth.avg, group=sigma, colour=sigma))
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
