# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR 
DUMPDIR = "/mnt/tmp/dump/NAVNP/"

# Linear
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon/"
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v/"

DUMPDIR <- '/home/stefano/HS/'
DIR = 'aggrLinear_epsilon/'
DIR = 'aggrLinear_RFull_epsilon_v/'


INTERACTIVE = FALSE
PATH = paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");

# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/img/"))
}
if (!file.exists(paste0(IMGPATH, "new/"))) {
  dir.create(file.path(PATH, "/img/new/"))
}
if (!file.exists(paste0(IMGPATH, "newpdist/"))) {
  dir.create(file.path(PATH, "/img/newpdist/"))
}


##############################
# Explanation of loaded files:
##############################
#
# - *_avg_all_split.csv contains the averages for each time step of all
#   the simulations split by each sigma level (~2000 * every sigma level).
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


params <- subset(params, select=-c(seed, run, attr_on_v, attrtype, noisetype,
                                   truth.x, truth.y, init.clusterradio,
                                   init.nclusters, tau,
                                   d1, B, d0, A, k, spacesize, spacedim,
                                   nagents, t.end, dt, timestamp))




macro <- read.table('clusters_macro.csv', head=TRUE, sep=",")

clu <- merge(params, macro, by=c("simname","simcount"))



clu$simname <- as.character(clu$simname)
clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
clu$simname <- as.factor(clu$simname)
clu$simcount <- as.factor(clu$simcount)
clu$t <- as.factor(clu$t)

# Loaded from another dir.
cl.R.full <- clu

cl <- clu[clu$t == 100 & clu$R < 1.4,]

p <- ggplot(cl, aes(R, count))
p <- p +  geom_bar(stat = "identity", position="dodge")
p

p <- ggplot(cl.R.full, aes(R, count))
p <- p +  geom_bar(stat = "identity", position="dodge")
p



p <- ggplot(cl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge")
p

p <- ggplot(clu[clu$t == 2000,], aes(alpha, count))
p <- p + geom_bar(stat = "identity", position="dodge")
p

p <- ggplot(clu[clu$t == 2000,], aes(epsilon, count))
p <- p + geom_bar(stat = "identity", position="dodge")
p

p <- ggplot(clu[clu$t == 2000,], aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge")
p





p <- p + facet_grid(sigma ~.)
p
