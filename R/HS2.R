# HS analysis
source("./init.R")

DIR = "test_t-2013-6-4-12-14"

DIR = "attrExpo_nv_rndseq_tm_Rleft" # ADD NEW

DUMPDIR = "/cluster/home/gess/balistef/matlab/hsnew/dump/"

PATH = paste0(DUMPDIR,"/",DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "/img/");
params <- read.table('params.csv', head=TRUE, sep=",")

# NORMAL READ.TABLE
#clusters <- read.table('clusters_macro.csv', head=TRUE, sep=",")

library(sqldf)
tic()
f1 <- file('clusters_macro_smaller.csv')
toc()

tic()
clusters <- sqldf("select * from f1", dbname = tempfile(), file.format = list(header = T, row.names = F))
toc()

# Micro commented for the moment
#clusters.micro <- read.table('clusters_micro.csv', head=TRUE, sep=",")




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
#clu$fromtruth.avg.cut <- cut(clu$fromtruth.avg, seq(0,1,0.1))
#clu$size.avg.cut <- cut(clu$size.avg, seq(0,100,5))
#clu$count.cut <- cut(clu$count, seq(0,100,5))

# lagging

#clu$fromtruth.lag <- lagg(clu$fromtruth.avg)
#clu$truthdiff <- (clu$fromtruth.avg - clu$fromtruth.lag)
#clu$ratediff <- clu$truthdiff / clu$fromtruth.avg
#clu$Rjump <- clu$truthdiff > 0.02

# SimCount

#a <- params[params$alpha == 0.4 & params$sigma == 0 & params$R == 0.07,]
#as.numeric(a$simcount)       

#a <- clu[clu$fromtruth.avg > 0.6,]

#max.dist.from.truth <- which(clu$fromtruth.avg == max(clu$fromtruth.avg))

#clu[max.dist.from.truth,]


# Cluster Macro in time

library(biglm)
library(bigmemory)
library(biganalytics)


x <- read.big.matrix("clusters_micro_cutoff.csv", header=TRUE,
                     backingfile="micro_cutoff.bin",
                     descriptorfile="micro_cutoff.desc",
                     col.names = c("simname", "simcount", "run", "t", "size", "speed", "move", "fromtruth"))


xdesc <- dget("micro_cutoff.desc")
x <- attach.big.matrix(xdesc)


tic()
xfit <- biglm.big.matrix(speed ~ size + t + size*t, data = x)
toc()

xcut <- cut(x[, "size"], c(3,5,10,20,50,75,100))
x <- cbind(x, xcut)

xsplit <- split(1:nrow(x), x[, "size"]);


plot(x[,"size"],x[,"speed"])

cl <- clusters

title = "Evolution of cluster count and size"
p.count <- ggplot(cl, aes(t))
p.count <- p.count + geom_jitter(aes(y = size.avg), alpha=.2)
p.count <- p.count + geom_smooth(aes(y = count, colour="count"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.avg, colour="size"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.sd, colour="std. size"), size=2)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
#p.count

title = "Average and cumulative space exploration"
p.explo <- ggplot(cl, aes(t))
p.explo <- p.explo + geom_jitter(aes(y = coverage), alpha=.2)
p.explo <- p.explo + geom_smooth(aes(y = coverage, colour="avg"), size=2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.cum, colour="cum"), size=2)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
#p.explo

title = "Mean and std. agents speed"
p.speed <- ggplot(cl, aes(t))
p.speed <- p.speed + geom_jitter(aes(y = speed.avg), alpha=.2)
p.speed <- p.speed + geom_smooth(aes(y = speed.avg, colour="avg"), size=2)
p.speed <- p.speed + geom_smooth(aes(y = speed.sd, colour="sd"), size=2)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
#p.speed

title = "Mean and std. agents movements"
p.move <- ggplot(cl, aes(t))
p.move <- p.move + geom_jitter(aes(y = move.avg), alpha=.2)
p.move <- p.move + geom_smooth(aes(y = move.avg, colour="avg"), size=2)
p.move <- p.move + geom_smooth(aes(y = move.sd, colour="sd"), size=2)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Displacement")
#p.move

title = "Mean and std. distance from truth"
p.truth <- ggplot(cl, aes(t))
p.truth <- p.truth + geom_jitter(aes(y = fromtruth.avg), alpha=.2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.avg, colour="avg"), size=2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.sd, colour="sd"), size=2)
p.truth <- p.truth + ggtitle(title) + xlab("Rounds") + ylab("Distance")
#print(p.truth)

jpeg(paste0(IMGPATH, "mycoolplot.jpeg"), width=800, height=480)
p <- grid.arrange(p.count, p.explo, p.speed, p.truth, p.move, ncol=3,
            main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()

stop("Message")
warning("Message")


# Cluster Micro

clusters.micro <- read.table('clusters_micro.csv', head=TRUE, sep=",")

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
  
