# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

DIR = "test_t-2013-6-4-12-14"

DIR = "attrExpo_nv_rndseq_tm_Rleft" # ADD NEW


DUMPDIR = "/opt/MATLAB_WORKSPACE/hs/dump/NEW/"
PATH = paste0(DUMPDIR,DIR)
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

# Cluster Macro in time

cl <- clusters

title = "Evolution of cluster count and size"
p.count <- ggplot(cl, aes(t))
p.count <- p.count + geom_jitter(aes(y = size.avg), alpha=.2)
p.count <- p.count + geom_smooth(aes(y = count, colour="count"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.avg, colour="size"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.sd, colour="std. size"), size=2)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
#p.count
ggsave("count_size.jpg", plot=p.count)

#title = "Average and cumulative space exploration"
#p.explo <- ggplot(cl, aes(t))
#p.explo <- p.explo + geom_jitter(aes(y = coverage), alpha=.2)
#p.explo <- p.explo + geom_smooth(aes(y = coverage, colour="avg"), size=2)
#p.explo <- p.explo + geom_smooth(aes(y = coverage.cum, colour="cum"), size=2)
#p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
#p.explo
#
#title = "Mean and std. agents speed"
#p.speed <- ggplot(cl, aes(t))
#p.speed <- p.speed + geom_jitter(aes(y = speed.avg), alpha=.2)
#p.speed <- p.speed + geom_smooth(aes(y = speed.avg, colour="avg"), size=2)
#p.speed <- p.speed + geom_smooth(aes(y = speed.sd, colour="sd"), size=2)
#p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
#p.speed
#
#title = "Mean and std. agents movements"
#p.move <- ggplot(cl, aes(t))
#p.move <- p.move + geom_jitter(aes(y = move.avg), alpha=.2)
#p.move <- p.move + geom_smooth(aes(y = move.avg, colour="avg"), size=2)
#p.move <- p.move + geom_smooth(aes(y = move.sd, colour="sd"), size=2)
#p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Displacement")
#p.move
#
#title = "Mean and std. distance from truth"
#p.truth <- ggplot(cl, aes(t))
#p.truth <- p.truth + geom_jitter(aes(y = fromtruth.avg), alpha=.2)
#p.truth <- p.truth + geom_smooth(aes(y = fromtruth.avg, colour="avg"), size=2)
#p.truth <- p.truth + geom_smooth(aes(y = fromtruth.sd, colour="sd"), size=2)
#p.truth <- p.truth + ggtitle(title) + xlab("Rounds") + ylab("Distance")
#print(p.truth)
#
#jpeg(paste0(IMGPATH, "mycoolplot.jpeg"), width=800, height=480)
#p <- grid.arrange(p.count, p.explo, p.speed, p.truth, p.move, ncol=3,
#            main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
#dev.off()


