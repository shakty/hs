# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR VELOCITY

DUMPDIR = "/home/stefano/hs/test/"
DIR = "NEWTEST-2013-12-8-17-49/"

DUMPDIR = "/home/stefano/"
DIR = ""


INTERACTIVE = FALSE
PATH = paste0(DUMPDIR, DIR, "aggr_velocities/")
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");

# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/img/"))
}


##############################
# Explanation of loaded files:
##############################
#
# - *_avg_all_split.csv contains the averages for each time step of all
#   the simulations diveded by each sigma level (~2000 * every sigma level).
#
# - *_avg_all.csv contains the averages for each time step (~2000) 
#
# - *.csv contains each ~20 snapshots of all simulations (every ~100 time step)
#   (20 * nCombinations of parameter sweep * levels of sigma)
#
#############################


##########
# PARAMS #
##########
params <- read.table('params.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)


params <- subset(params, select=-c(seed, attr_on_v, attrtype, noisetype,
                                   truth.x, truth.y, init.clusterradio,
                                   init.nclusters, tau,
                                   d1, B, d0, A, k, spacesize, spacedim,
                                   nagents, t.end, dt, timestamp))

################################
## HEATMAPS: 1 res every X time steps (default 100) per simulation (default 20 x 2970)
################################

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
tr$simname <- as.factor(tr$simname)
tr$simcount <- as.factor(tr$simcount) # or N?
tr$run <- as.factor(tr$run)

## MACRO: Heatmap

v1 <- "R"
v2 <- "alpha"
v3 <- "sigma"
paramsData <- params


if (!file.exists(paste0(IMGPATH, "velocity/"))) {
  dir.create(file.path(IMGPATH, "velocity/"))
}

# Visualize the clusters.

# Merging clusters macro and params: clu
########################################
clu <- merge(params, macro, by=c("simname","simcount","run"))
clu$simname <- as.character(clu$simname)
clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
clu$simname <- as.factor(clu$simname)

# Factorising
for (n in names(clu[1:7])) {
  clu[, n] <- as.factor(clu[, n])      
}

cl <- clu

# VELOCITY ANALYSIS
################

# SIZE.MAX distribution by sigma and epsilon
title = "Distributions of size of the biggest cluster at the end of simulation by sigma and epsilon"
p <- ggplot(cl[cl$t == 2000,], aes(size.max))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Velocities") + ylab("Sigmas")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/sizemax_distributions.jpg"), plot=p)

# SIZE.MAX heatmaps + R and Alpha
title = "Distributions of size of the biggest cluster at the end of simulation by sigma, epsilon (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=size.max), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Velocitie") + ylab("Sigmas")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/sizemax_r_alpha_ht.jpg"), plot=p)

# BIGC PDIST distribution by sigma and epsilon
title = "Distributions of pairwise distances within biggest cluster at the end of simulation by sigma and epsilon"
p <- ggplot(cl[cl$t == 2000,], aes(bigc.pdist.mean))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Velocities") + ylab("Sigmas")
if (INTERACTIVE) {
  p
}

ggsave(filename=paste0(IMGPATH, "noise/bigcpdist_distributions.jpg"), plot=p)

# BIGC PDIST heatmaps + R and Alpha
title = "Distributions of pairwise distances within biggest cluster at the end of simulation by sigma, epsilon (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=bigc.pdist.mean), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Velocities") + ylab("Sigmas")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/bigcpdist_r_alpha_ht.jpg"), plot=p)


# Merging truthradius and params: clu
#####################################
clu <- merge(params, tr, by=c("simname","simcount","run"))
# Factorising
for (n in names(clu[1:7])) {
  clu[, n] <- as.factor(clu[, n])      
}

cl <- clu

# CONSENSUS distribution by sigma and epsilon
title = "Consensus on truth at the end of simulation by sigma and epsilon"
p <- ggplot(cl[cl$t == 2000,], aes(consensus))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Epsilons")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/consensus_on_t_distributions.jpg"), plot=p)

# CONSENSUS heatmaps + R and Alpha
title = "Consensus on truth at the end of simulation by sigma, velocity (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=consensus), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/consensus_on_t_r_alpha_ht.jpg"), plot=p)

# RADIUSES distribution by sigma and velocity
cl.melted <- melt(cl,
                  measure.vars=c("r_01","r_05","r_1","r_25", "r_4","r_out"),
                  variable_name="condition")

a <- givemeSummary(cl.melted, "value", c("condition","init.vscaling","sigma"), TRUE)

title = "Average agent counts in each range (+Std.Err.) by sigma and velocity"
limits <- aes(ymax = value + value.se, ymin=value - value.se)
p <- ggplot(a, aes(x = condition, y=value))
p <- p + geom_bar(aes(color=condition, fill=condition))
p <- p + geom_errorbar(limits, width=0.2)
p <- p + facet_grid(init.vscaling ~ sigma)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/radiuses_hist.jpg"), plot=p)


# Merging agents and params: clu.
##################################
clu <- merge(params, agents, by=c("simname","simcount","run"))   
# Factorising
for (n in names(clu[1:7])) {
  clu[, n] <- as.factor(clu[, n])      
}

# NOISE ANALYSIS
################

cl <- clu

# CUM.COVERAGE distribution by sigma and velocity
title = "Distributions of cumulative space exploration the end of simulation by sigma and velocity"
p <- ggplot(cl[cl$t == 2000,], aes(coverage.cum))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/cumcov_distributions.jpg"), plot=p)

# CUM.COVERAGE heatmaps + R and Alpha
title = "Distributions of cumulative space exploration at the end of simulation by sigma, velocity (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=coverage.cum), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/cumcov_r_alpha_ht.jpg"), plot=p)


# FROMTRUTH distribution by sigma and velocity
title = "Distributions of average distance from truth the end of simulation by sigma and velocity"
p <- ggplot(cl[cl$t == 2000,], aes(fromtruth.avg))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/fromtruth_distributions.jpg"), plot=p)

# FROMTRUTH heatmaps + R and Alpha
title = "Distributions of average distance from truth at the end of simulation by sigma, velocity (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=fromtruth.avg), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/fromtruth_r_alpha_ht.jpg"), plot=p)


# MOVS distribution by sigma and velocity
title = "Distributions of average movements at the end of simulation by sigma and velocity"
p <- ggplot(cl[cl$t == 2000 & cl$init.vscaling != 100,], aes(move.avg))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/move_distributions.jpg"), plot=p)

# MOVS heatmaps + R and Alpha
title = "Distributions of average movements at the end of simulation by sigma, velocity (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=move.avg), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/move_r_alpha_ht.jpg"), plot=p)

# SPEED distribution by sigma and velocity
title = "Distributions of average speed at the end of simulation by sigma and velocity"
#p <- ggplot(cl[cl$t == 2000 & cl$init.vscaling != 100,], aes(speed.avg))
p <- ggplot(cl[cl$t == 2000,], aes(speed.avg))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/speed_distributions.jpg"), plot=p)

# SPEED heatmaps + R and Alpha
title = "Distributions of average speed at the end of simulation by sigma, velocity (outer), and alpha and R (inner)"
#p <- ggplot(cl[cl$t == 2000 & cl$init.vscaling != 100,], aes(x=R, y=alpha))
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=speed.avg), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/speed_r_alpha_ht.jpg"), plot=p)

# MEAN PDIST distribution by sigma and velocity
title = "Distributions of average pair-wise distance at the end of simulation by sigma and velocity"
p <- ggplot(cl[cl$t == 2000,], aes(pdist.mean))
p <- p + geom_bar(aes(color=sigma))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "noise/pdist_distributions.jpg"), plot=p)

# MEAN PDIST heatmaps + R and Alpha
title = "Distributions of average speed at the end of simulation by sigma, velocity (outer), and alpha and R (inner)"
p <- ggplot(cl[cl$t == 2000,], aes(x=R, y=alpha))
p <- p + geom_tile(aes(fill=pdist.mean), color="white")
p <- p + scale_fill_continuous(low='lightblue',high='red')
p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))
p <- p + facet_grid(sigma ~ init.vscaling)
p <- p + ggtitle(title) + xlab("Sigmas") + ylab("Velocities")
if (INTERACTIVE) {
  p
}

ggsave(filename=paste0(IMGPATH, "noise/pdist_r_alpha_ht.jpg"), plot=p)
