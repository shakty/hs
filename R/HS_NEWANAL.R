# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR VELOCITY
DUMPDIR = "/home/stefano/hs/test/"
DIR = "NEWTEST-2013-12-8-17-49/"

DUMPDIR = "/home/stefano/"
DIR = ""

INTERACTIVE = FALSE
PATH = paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");

# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/img/"))
}

######################################
## CLUSTER MACRO Average per round per simulation: Loading
######################################

macro <- read.table('clusters_macro_avg_all_split.csv', head=TRUE, sep=",")
macro$simname <- as.character(macro$simname)
macro$simname <- substr(macro$simname, nchar(macro$simname)-1, nchar(macro$simname))
macro$simname <- as.factor(macro$simname)
macro$simcount <- as.factor(macro$simcount)
macro$t <- as.factor(macro$t)

cl <- macro
#cl <- subset(macro, t %% 100 == 0)
#cl$t <- as.factor(cl$t)

# SIZE with SD
title = "Evolution of average cluster size (+Std.Dev) by sigma"
p.size <- ggplot(cl, aes(t,group=simname))
p.size <- p.size + geom_area(aes(y = meansize.avg, colour=simname, group=simname))
p.size <- p.size + geom_area(aes(y = meansize.sd, colour=simname, group=simname, fill=simname))
p.size <- p.size + facet_grid(simname~.,margins=F)
p.size <- p.size + ggtitle(title) + xlab("Rounds") + ylab("Average cluster size")
if (INTERACTIVE) {
  p.size
}

# SIZE with SMOOTH (loess)
title = "Evolution of average cluster size (Loess + Std.Errs) by sigma"
p.size.se <- ggplot(cl, aes(t,group=simname))
p.size.se <- p.size.se + geom_jitter(aes(y = meansize.avg, colour=simname, group=simname), alpha=0.2)
p.size.se <- p.size.se + geom_smooth(aes(y = meansize.avg, colour=simname, group=simname))
p.size.se <- p.size.se + ggtitle(title) + xlab("Rounds") + ylab("Average cluster size")
if (INTERACTIVE) {
  p.size.se
}

# SIZE MAX with SD
title = "Evolution of average size of the biggest cluster (+Std.Dev) by sigma"
p.size.max <- ggplot(cl, aes(t,group=simname))
p.size.max <- p.size.max + geom_area(aes(y = maxsize.avg, colour=simname, group=simname))
p.size.max <- p.size.max + geom_area(aes(y = maxsize.sd, colour=simname, group=simname, fill=simname))
p.size.max <- p.size.max + facet_grid(simname~.,margins=F)
p.size.max <- p.size.max + ggtitle(title) + xlab("Rounds") + ylab("Agents in biggest cluster")
if (INTERACTIVE) {
  p.size.max
}

# SIZE MAX points with SMOOTH (loess)
title = "Evolution of average size of the biggest (Loess + Std.Errs) by sigma"
p.size.max.se <- ggplot(cl, aes(t,group=simname))
p.size.max.se <- p.size.max.se + geom_jitter(aes(y = maxsize.avg, colour=simname, group=simname), alpha=0.2)
p.size.max.se <- p.size.max.se + geom_smooth(aes(y = maxsize.avg, colour=simname, group=simname))
p.size.max.se <- p.size.max.se + facet_grid(simname~.,margins=F)
p.size.max.se <- p.size.max.se + ggtitle(title) + xlab("Rounds") + ylab("Agents in biggest cluster")
if (INTERACTIVE) {
  p.size.max.se
}

# COUNT with SD
title = "Evolution of cluster counts (+Std.Dev) by sigma"
p.count <- ggplot(cl, aes(t,group=simname))
p.count <- p.count + geom_area(aes(y = count.avg, colour=simname, group=simname))
p.count <- p.count + geom_area(aes(y = count.sd, colour=simname, group=simname, fill=simname))
p.count <- p.count + facet_grid(simname~.,margins=F)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents in biggest cluster")
if (INTERACTIVE) {
  p.count
}

# COUNT points with SMOOTH (loess)
title = "Evolution of cluster counts (Loess + Std.Errs) by sigma"
p.count.se <- ggplot(cl, aes(t, group=simname))
p.count.se <- p.count.se + geom_jitter(aes(y = count.avg, colour=simname, group=simname), alpha=0.2)
p.count.se <- p.count.se + geom_smooth(aes(y = count.avg, colour=simname, group=simname))
p.count.se <- p.count.se + ggtitle(title) + xlab("Rounds") + ylab("Agents in biggest cluster")
if (INTERACTIVE) {
  p.count.se
}

# WITHIN-CLUSTER CONSENSUS with SD
title = "Evolution of within-cluster consensus (+Std.Dev) by sigma"
p.bigc.consensus <- ggplot(cl, aes(t,group=simname))
p.bigc.consensus <- p.bigc.consensus + geom_area(aes(y = bigc.pdist.avg, colour=simname, group=simname))
p.bigc.consensus <- p.bigc.consensus + geom_area(aes(y = bigc.pdist.sd, colour=simname, group=simname, fill=simname))
p.bigc.consensus <- p.bigc.consensus + facet_grid(simname~.,margins=F)
p.bigc.consensus <- p.bigc.consensus + ggtitle(title) + xlab("Rounds") + ylab("Average pair-wise distance in biggest cluster")
if (INTERACTIVE) {
  p.bigc.consensus
}

#  WITHIN-CLUSTER CONSENSUS with SMOOTH (loess)
title = "Evolution of within-cluster consensus (Loess + Std.Errs) by sigma"
p.bigc.consensus.se <- ggplot(cl, aes(t, group = simname))
p.bigc.consensus.se <- p.bigc.consensus.se + geom_jitter(aes(y = bigc.pdist.avg, colour=simname, group=simname), alpha=0.2)
p.bigc.consensus.se <- p.bigc.consensus.se + geom_smooth(aes(y = bigc.pdist.avg, colour=simname, group=simname))
p.bigc.consensus.se <- p.bigc.consensus.se + ggtitle(title) + xlab("Rounds") + ylab("Average pair-wise distance in biggest cluster")
if (INTERACTIVE) {
  p.bigc.consensus.se
}


### PLOTTING ALL OF THEM TOGETHER!!
jpeg(paste0(IMGPATH, "t_clusters_macro_overview_split.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.size, p.size.se,
                  p.size.max, p.size.max.se,
                  p.count, p.count.se,
                  p.bigc.consensus, p.bigc.consensus.se,
                  ncol=4,
                  main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()


######################################
## AGENTS Average per round per simulation: Loading
######################################

agents <- read.table('agents_avg_all_split.csv', head=TRUE, sep=",")
agents$simname <- as.character(agents$simname)
agents$simname <- substr(agents$simname, nchar(agents$simname)-1, nchar(agents$simname))
agents$simname <- as.factor(agents$simname)
agents$simcount <- as.factor(agents$N)
agents$t <- as.factor(agents$t)

ag <- agents

# CUM EXPLORATION + SD
title = "Cumulative exploration (+Std.Dev.) by sigma"
p.explo.cum <- ggplot(ag, aes(t))
p.explo.cum <- p.explo.cum + geom_area(aes(y = coverage.cum.avg, colour=simname, group=simname))
p.explo.cum <- p.explo.cum + geom_area(aes(y = coverage.cum.sd, colour=simname, group=simname, fill=simname))
p.explo.cum <- p.explo.cum + facet_grid(simname~.,margins=F)
p.explo.cum <- p.explo.cum + ggtitle(title) + xlab("Rounds") + ylab("Cumulative % of explored space")
if (INTERACTIVE) {
  p.explo.cum
}

# CUM EXPLORATION + SE
title = "Evolution of cumulative (Loess + Std.Errs) exploration by sigma"
p.explo.cum.se <- ggplot(ag, aes(t, group=simname))
p.explo.cum.se <- p.explo.cum.se + geom_jitter(aes(y = coverage.cum.avg, colour=simname, group=simname), alpha=0.2)
p.explo.cum.se <- p.explo.cum.se + geom_smooth(aes(y = coverage.cum.avg, colour=simname, group=simname))
p.explo.cum.se <- p.explo.cum.se + ggtitle(title) + xlab("Rounds") + ylab("Cumulative % of explored space")
if (INTERACTIVE) {
  p.explo.cum.se
}

# POINT EXPLORATION + SD
title = "Point exploration (+Std.Dev.) by sigma"
p.explo <- ggplot(ag, aes(t))
p.explo <- p.explo + geom_area(aes(y = coverage.avg, colour=simname, group=simname))
p.explo <- p.explo + geom_area(aes(y = coverage.sd, colour=simname, group=simname, fill=simname))
p.explo <- p.explo + facet_grid(simname~.,margins=F)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("% of explored space at a given round")
if (INTERACTIVE) {
  p.explo
}

# POINT EXPLORATION + SE
title = "Evolution of point exploration (Loess + Std.Errs) by sigma"
p.explo.se <- ggplot(ag, aes(t,group=simname))
p.explo.se <- p.explo.se + geom_jitter(aes(y = coverage.avg, colour=simname, group=simname), alpha=0.2)
p.explo.se <- p.explo.se + geom_smooth(aes(y = coverage.avg, colour=simname, group=simname))
p.explo.se <- p.explo.se + ggtitle(title) + xlab("Rounds") + ylab("% of explored space at a given round")
if (INTERACTIVE) {
  p.explo.se
}

# SPEED + SD
title = "Mean and std. agents speed"
p.speed <- ggplot(ag, aes(t))
p.speed <- p.speed + geom_area(aes(y = speed.avg, colour=simname, group=simname))
p.speed <- p.speed + geom_area(aes(y = speed.sd, colour=simname, group=simname, fill=simname))
p.speed <- p.speed + facet_grid(simname~.,margins=F)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed
}

# SPEED + SE
title = "Evolution of speed (Loess + Std.Errs) by sigma"
p.speed.se <- ggplot(ag, aes(t,group=simname))
p.speed.se <- p.speed.se + geom_jitter(aes(y = speed.avg, colour=simname, group=simname), alpha=0.2)
p.speed.se <- p.speed.se + geom_smooth(aes(y = speed.avg, colour=simname, group=simname))
p.speed.se <- p.speed.se + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed.se
}

#MOVEMENTS + SD
title = "Evolution of movements (+Std.Dev.) by sigma"
p.move <- ggplot(ag, aes(t))
p.move <- p.move + geom_area(aes(y = move.avg, colour=simname, group=simname))
p.move <- p.move + facet_grid(simname~.,margins=F)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Move")
if (INTERACTIVE) {
  p.move
}

# MOVEMENTS + SE
title = "Evolution of movements (Loess + Std.Errs) by sigma"
p.move.se <- ggplot(ag, aes(t, group=simname))
p.move.se <- p.move.se + geom_jitter(aes(y = move.avg, colour=simname, group=simname))
p.move.se <- p.move.se + geom_smooth(aes(y = move.avg, colour=simname, group=simname))
p.move.se <- p.move.se + ggtitle(title) + xlab("Rounds") + ylab("Move")
if (INTERACTIVE) {
  p.move.se
}


# FROM TRUTH + SD
title = "Evolution of distance from truth (+Std.Dev.) by sigma"
p.ft <- ggplot(ag, aes(t))
p.ft <- p.ft + geom_area(aes(y = fromtruth.avg, colour=simname, group=simname))
p.ft <- p.ft + geom_area(aes(y = fromtruth.sd, colour=simname, group=simname, fill=simname))
p.ft <- p.ft + facet_grid(simname~.,margins=F)
p.ft <- p.ft + ggtitle(title) + xlab("Rounds") + ylab("Distance from truth")
if (INTERACTIVE) {
  p.ft
}

# FROM TRUTH + SE
title = "Evolution of distance from truth (Loess + Std.Errs) by sigma"
p.ft.se <- ggplot(ag, aes(t,group=simname))
p.ft.se <- p.ft.se + geom_jitter(aes(y = fromtruth.avg, colour=simname, group=simname))
p.ft.se <- p.ft.se + geom_smooth(aes(y = fromtruth.avg, colour=simname, group=simname))
p.ft.se <- p.ft.se + ggtitle(title) + xlab("Rounds") + ylab("Distance from truth")
if (INTERACTIVE) {
  p.ft.se
}

# WITHIN-CLUSTER CONSENSUS with SD
title = "Evolution of between agents consensus (+Std.Dev) by sigma"
p.consensus <- ggplot(ag, aes(t,group=simname))
p.consensus <- p.consensus + geom_area(aes(y = pdist.avg, colour=simname, group=simname))
p.consensus <- p.consensus + geom_area(aes(y = pdist.sd, colour=simname, group=simname, fill=simname))
p.consensus <- p.consensus + facet_grid(simname~.,margins=F)
p.consensus <- p.consensus + ggtitle(title) + xlab("Rounds") + ylab("Average pair-wise distance")
if (INTERACTIVE) {
  p.consensus
}

#  WITHIN-CLUSTER CONSENSUS with SMOOTH (loess)
title = "Evolution of between agents consensus (Loess + Std.Errs) by sigma"
p.consensus.se <- ggplot(ag, aes(t, group = simname))
p.consensus.se <- p.consensus.se + geom_jitter(aes(y = pdist.avg, colour=simname, group=simname), alpha=0.2)
p.consensus.se <- p.consensus.se + geom_smooth(aes(y = pdist.avg, colour=simname, group=simname))
p.consensus.se <- p.consensus.se + ggtitle(title) + xlab("Rounds") + ylab("Average pair-wise distance")
if (INTERACTIVE) {
  p.consensus.se
}

### PLOTTING ALL OF THEM TOGETHER!!
jpeg(paste0(IMGPATH, "t_agents_overview_split.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.explo.cum, p.explo.cum.se,
                  p.explo, p.explo.se,
                  p.speed, p.speed.se,
                  p.move,p.move.se,
                  p.ft, p.ft.se,
                  p.consensus, p.consensus.se,
                  ncol=4,
                  main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()


######################################
## TRUTHRADIUS Average per round per simulation: Loading
######################################

tr <- read.table('truthradius_avg_all_split.csv', head=TRUE, sep=",")
tr$simname <- as.character(tr$simname)
tr$simname <- substr(tr$simname, nchar(tr$simname)-1, nchar(tr$simname))
tr$simname <- as.factor(tr$simname)
tr$simcount <- as.factor(tr$N)
tr$t <- as.factor(tr$t)


# Agents in radiuses
title = "Agents count in each radius"
p.radiuses <- ggplot(tr, aes(t))
p.radiuses <- p.radiuses + geom_area(aes(y = r_01_mean, group=simname, fill="0.01"))
p.radiuses <- p.radiuses + geom_area(aes(y = r_05_mean, group=simname, fill="0.05"))
p.radiuses <- p.radiuses + geom_area(aes(y = r_1_mean, group=simname, fill="0.1"))
p.radiuses <- p.radiuses + geom_area(aes(y = r_25_mean, group=simname, fill="0.25"))
p.radiuses <- p.radiuses + geom_area(aes(y = r_4_mean, group=simname, fill="0.4"))
p.radiuses <- p.radiuses + geom_area(aes(y = r_out_mean, group=simname, fill="out"))
p.radiuses <- p.radiuses + facet_grid(simname~.,margins=F)
p.radiuses <- p.radiuses + ggtitle(title) + xlab("Rounds") + ylab("Cumulative % of explored space")
if (INTERACTIVE) {
  p.radiuses
}

# Agents in radiuses SMOOTH
title = "Agents count in each radius"
p.radiuses.smooth <- ggplot(tr, aes(t))
p.radiuses.smooth <- p.radiuses.smooth + geom_smooth(aes(y = r_01_mean, group=simname, color="0.01"))
p.radiuses.smooth <- p.radiuses.smooth + geom_smooth(aes(y = r_05_mean, group=simname, color="0.05"))
p.radiuses.smooth <- p.radiuses.smooth + geom_smooth(aes(y = r_1_mean, group=simname, color="0.1"))
p.radiuses.smooth <- p.radiuses.smooth + geom_smooth(aes(y = r_25_mean, group=simname, color="0.25"))
p.radiuses.smooth <- p.radiuses.smooth + geom_smooth(aes(y = r_4_mean, group=simname, color="0.4"))
p.radiuses.smooth <- p.radiuses.smooth + geom_smooth(aes(y = r_out_mean, group=simname, color="out"))
p.radiuses.smooth <- p.radiuses.smooth + facet_grid(simname~.,margins=F)
p.radiuses.smooth <- p.radiuses.smooth + ggtitle(title) + xlab("Rounds") + ylab("Cumulative % of explored space")
if (INTERACTIVE) {
  p.radiuses.smooth
}

# Agents reaches consensus (share of simulations where agents reached consensus).
# A certain share of agents (default 2/3) must stay in the inner radius (default 0.01) for a long enough number of time steps (default 20).
title = "Consensus reached at time t."
p.consensus.reach <- ggplot(tr, aes(t))
p.consensus.reach <- p.consensus.reach + geom_jitter(aes(y = consensus.avg, colour=simname, group=simname), alpha=0.2)
p.consensus.reach <- p.consensus.reach + geom_smooth(aes(y = consensus.avg, group=simname, color=simname))
p.consensus.reach <- p.consensus.reach + facet_grid(simname ~ ., margins=F)
p.consensus.reach <- p.consensus.reach + ggtitle(title) + xlab("Rounds") + ylab("Share of simulations that reached consensus at time t.")
if (INTERACTIVE) {
  p.consensus.reach
}

### PLOTTING ALL OF THEM TOGETHER!!
jpeg(paste0(IMGPATH, "t_truthradius_overview_split.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.radiuses,
                  p.radiuses.smooth,
                  p.consensus.reach,
                  ncol=3,
                  main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()

################################
## HEATMAPS: 1 res every X time steps (default 100) per simulation (default 20 x 2970)
################################

params <- read.table('params.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)

params <- subset(params, select=-c(seed, attr_on_v, attrtype, noisetype,
                                   truth.x, truth.y, init.clusterradio,
                                   init.nclusters, init.vscaling, tau,
                                   d1, B, d0, A, k, spacesize, spacedim,
                                   nagents, t.end, dt, timestamp))

macro <- read.table('clusters_macro.csv', head=TRUE, sep=",")
macro$simname <- as.factor(macro$simname)
macro$simcount <- as.factor(macro$simcount)
macro$run <- as.factor(macro$run)

macro <- subset(macro, select=-c(move.avg,
                                     speed.avg,
                                     move.sd,
                                     speed.sd,
                                     fromtruth.avg,
                                     fromtruth.sd))

agents <- read.table('agents.csv', head=TRUE, sep=",")
agents$simname <- as.factor(agents$simname)
agents$simcount <- as.factor(agents$simcount)
agents$run <- as.factor(agents$run)

tr <- read.table('truthradius.csv', head=TRUE, sep=",")
agents$simname <- as.factor(agents$simname)
agents$simcount <- as.factor(agents$simcount) # or N?
agents$run <- as.factor(agents$run)

## MACRO: Heatmap

v1 <- "R"
v2 <- "alpha"
v3 <- "sigma"
paramsData <- params

if (!file.exists(paste0(IMGPATH, "ft/"))) {
  dir.create(paste0(IMGPATH, "ft/"))
}
if (!file.exists(paste0(IMGPATH, "count/"))) {
  dir.create(paste0(IMGPATH, "count/"))
}   
if (!file.exists(paste0(IMGPATH, "size/"))) {
  dir.create(file.path(IMGPATH, "size/"))
}
if (!file.exists(paste0(IMGPATH, "sizemax/"))) {
  dir.create(file.path(IMGPATH, "sizemax/"))
}
if (!file.exists(paste0(IMGPATH, "speed/"))) {
  dir.create(file.path(IMGPATH, "speed/"))
}
if (!file.exists(paste0(IMGPATH, "move/"))) {
  dir.create(file.path(IMGPATH, "move/"))
}    
if (!file.exists(paste0(IMGPATH, "cumcov/"))) {
  dir.create(file.path(IMGPATH, "cumcov/"))
}
if (!file.exists(paste0(IMGPATH, "cov/"))) {
  dir.create(file.path(IMGPATH, "cov/"))
} 
if (!file.exists(paste0(IMGPATH, "pdist/"))) {
  dir.create(file.path(IMGPATH, "pdist/"))
}
if (!file.exists(paste0(IMGPATH, "bigcpdist/"))) {
  dir.create(file.path(IMGPATH, "bigcpdist/"))
}
if (!file.exists(paste0(IMGPATH, "consensus/"))) {
  dir.create(file.path(IMGPATH, "consensus/"))
}

# Visualize the clusters.

# Merging clusters macro and params: clu
clu <- merge(params, macro, by=c("simname","simcount","run"))
# Factorising
for (n in names(clu[1:7])) {
  clu[, n] <- as.factor(clu[, n])      
}

data <- clu
idx = 1;
for (t in unique(clu$t)) {
  data <- clu[clu$t == t,]
  # Count
  pt.c <- heatmapFacets_count(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"count/count_",sprintf("%04d",idx),".jpg"),plot=pt.c$p)
  # Size Mean
  pt.s <- heatmapFacets_size(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"size/size_",sprintf("%04d",idx),".jpg"),plot=pt.s$p)
  # Size Max
  pt.sm <- heatmapFacets_sizemax(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"sizemax/sizemax_",sprintf("%04d",idx),".jpg"),plot=pt.sm$p)
  # Pairwise-distance in biggest cluster
  pt.bigcpdist <- heatmapFacets_bigcpdist(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"size/size_",sprintf("%04d",idx),".jpg"),plot=pt.bigcpdist$p) 
  idx = idx + 1
}

# Merging truthradius and params: clu
clu <- merge(params, tr, by=c("simname","simcount","run"))
# Factorising
for (n in names(clu[1:7])) {
  clu[, n] <- as.factor(clu[, n])      
}

data <- clu
idx = 1;
for (t in unique(clu$t)) {
  data <- clu[clu$t == t,]
  # Consensus
  pt.consensus <- heatmapFacets_consensus(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"consensus/consensus_", sprintf("%04d",idx),".jpg"), plot=pt.consensus$p)
  idx = idx + 1
}

# Merging agents and params: clu.
clu <- merge(params, agents, by=c("simname","simcount","run"))   
# Factorising
for (n in names(clu[1:7])) {
  clu[, n] <- as.factor(clu[, n])      
}
data <- clu
idx = 1;

for (t in unique(clu$t)) {
  data <- clu[clu$t == t,]
  # From truth
  pt.ft <- heatmapFacets_fromtruth(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"ft/ft_",sprintf("%04d",idx),".jpg"),plot=pt.ft$p)
  # Move
  pt.m <- heatmapFacets_move(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"move/move_",sprintf("%04d",idx),".jpg"),plot=pt.m$p)
  # Speed
  pt.sp <- heatmapFacets_speed(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"speed/speed_",sprintf("%04d",idx),".jpg"),plot=pt.sp$p)
  # Cum Cov
  pt.ccov <- heatmapFacets_cumcoverage(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"cumcov/cumcov_",sprintf("%04d",idx),".jpg"),plot=pt.ccov$p)
  # Cov Instantaneous
  pt.cov <- heatmapFacets_coverage(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"cov/cov_",sprintf("%04d",idx),".jpg"),plot=pt.cov$p)
  # PDist
  pt.pdist <- heatmapFacets_pdist(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"pdist/pdist_",sprintf("%04d",idx),".jpg"),plot=pt.pdist$p)
  idx = idx + 1
}

# Creating an overview of the last frame of all indexes
jpeg(paste0(IMGPATH, "ht_overview.jpeg"), width=2048, height=768)
p <- grid.arrange(pt.s$p, pt.sm$p, pt.c$p, pt.bigcpdist$p,
                  pt.m$p, pt.sp$p, pt.ccov$p, pt.cov$p,
                  pt.ft$p, pt.pdist$p, pt.consensus$p,
                  ncol=4,
                  main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()

# Making a video out of the images
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/ft/ft_%04d.jpg img/ft/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/speed/speed_%04d.jpg img/speed/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/move/move_%04d.jpg img/move/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/size/size_%04d.jpg img/size/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/cumcov/cumcov_%04d.jpg img/cumcov/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/cov/cov_%04d.jpg img/cov/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/count/count_%04d.jpg img/count/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/sizemax/sizemax_%04d.jpg img/sizemax/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/pdist/pdist_%04d.jpg img/pdist/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/bigcpdist/bigcpdist_%04d.jpg img/bigcpdist/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/consensus/consensus_%04d.jpg img/consensus/movie.avi')


################
## MACRO AVG ALL
################

# if heatmap, needs to rename count to count.avg, maybe others too (???)

macro <- read.table('clusters_macro_avg_all.csv', head=TRUE, sep=",")
agents <- read.table('agents_avg_all.csv', head=TRUE, sep=",")
tr <- read.table('truthradius_avg_all.csv', head=TRUE, sep=",")

# attach(macro.avg.all)

# COUNT AND SIZE
title = "Evolution of cluster count and size (Loess + Std.Errs)"
p.count <- ggplot(macro, aes(t))
p.count <- p.count + geom_smooth(aes(y = meansize.avg, colour="mean size"), size=2)
p.count <- p.count + geom_smooth(aes(y = maxsize.avg, colour="max size"), size=2)
p.count <- p.count + geom_smooth(aes(y = count.avg, colour="count"), size=2)
p.count <- p.count + geom_smooth(aes(y = meansize.sd, colour="std. mean size"), size=2)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
if (INTERACTIVE) {
  p.count
}

#COVERAGE
title = "Average and cumulative space exploration"
p.explo <- ggplot(agents, aes(t))
p.explo <- p.explo + geom_jitter(aes(y = coverage.avg), alpha=.2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.avg, colour="avg"), size=2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.cum.avg, colour="cum"), size=2)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
if (INTERACTIVE) {
  p.explo
}

# SPEED
title = "Mean and std. agents speed"
p.speed <- ggplot(agents, aes(t))
p.speed <- p.speed + geom_jitter(aes(y = speed.avg), alpha=.2)
p.speed <- p.speed + geom_smooth(aes(y = speed.avg, colour="avg"), size=2)
p.speed <- p.speed + geom_smooth(aes(y = speed.sd, colour="sd"), size=2)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed
}

#MOVEMENTS
title = "Mean and std. agents movements"
p.move <- ggplot(agents, aes(t))
p.move <- p.move + geom_jitter(aes(y = move.avg), alpha=.2)
p.move <- p.move + geom_smooth(aes(y = move.avg, colour="avg"), size=2)
p.move <- p.move + geom_smooth(aes(y = move.sd, colour="sd"), size=2)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Displacement")
if (INTERACTIVE) {
  p.move
}

# FROM TRUTH
title = "Mean and std. distance from truth"
p.truth <- ggplot(agents, aes(t))
p.truth <- p.truth + geom_jitter(aes(y = fromtruth.avg), alpha=.2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.avg, colour="avg"), size=2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.sd, colour="sd"), size=2)
p.truth <- p.truth + ggtitle(title) + xlab("Rounds") + ylab("Distance")
if (INTERACTIVE) {
  p.truth
}

# PDIST
title = "Mean and std. consensus (pair-wise distance agents)"
p.pdist <- ggplot(agents, aes(t))
p.pdist <- p.pdist + geom_jitter(aes(y = pdist.avg), alpha=.2)
p.pdist <- p.pdist + geom_smooth(aes(y = pdist.avg, colour="avg"), size=2)
p.pdist <- p.pdist + geom_smooth(aes(y = pdist.sd, colour="sd"), size=2)
p.pdist <- p.pdist + ggtitle(title) + xlab("Rounds") + ylab("Pair-wise Distance")
if (INTERACTIVE) {
  p.pdist
}

# BIGC.PDIST
title = "Mean and std. consensus big cluster (pair-wise distance agents)"
p.bigc.pdist <- ggplot(macro, aes(t))
p.bigc.pdist <- p.bigc.pdist + geom_jitter(aes(y = bigc.pdist.avg), alpha=.2)
p.bigc.pdist <- p.bigc.pdist + geom_smooth(aes(y = bigc.pdist.avg, colour="avg"), size=2)
p.bigc.pdist <- p.bigc.pdist + geom_smooth(aes(y = bigc.pdist.sd, colour="sd"), size=2)
p.bigc.pdist <- p.bigc.pdist + ggtitle(title) + xlab("Rounds") + ylab("Pair-wise Distance")
if (INTERACTIVE) {
  p.bigc.pdist
}

# CONSENSUS on TRUTH
title = "Mean and std. consensus big cluster (pair-wise distance agents)"
p.consensus <- ggplot(tr, aes(t))
p.consensus <- p.consensus + geom_jitter(aes(y = consensus.avg), alpha=.2)
p.consensus <- p.consensus + geom_smooth(aes(y = consensus.avg, colour="avg"), size=2)
p.consensus <- p.consensus + geom_smooth(aes(y = consensus.sd, colour="sd"), size=2)
p.consensus <- p.consensus + ggtitle(title) + xlab("Rounds") + ylab("Pair-wise Distance")
if (INTERACTIVE) {
  p.consensus
}

# grid.Arrange
jpeg(paste0(IMGPATH, "t_avg_overview.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.count, p.explo, p.speed, p.truth,
                  p.move, p.pdist, p.bigc.pdist, p.consensus,
                  ncol=4,
                  main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()

