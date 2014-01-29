# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR 
DUMPDIR = "/mnt/tmp/dump/NAVNP/"

# Linear
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon/"
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v/"

DUMPDIR <- '/home/stefano/HS/'

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


loadData <- function(DUMPDIR, DIR) {

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


  ##########
  # PARAMS #
  ##########
  params <- read.table('params.csv', head=TRUE, sep=",")

  params$simname <- as.factor(params$simname)
  params$simcount <- as.factor(params$simcount)
  params$run <- as.factor(params$run)
  

  params <- subset(params, select=-c(seed, attr_on_v, attrtype, noisetype,
                                     truth.x, truth.y, init.clusterradio,
                                     init.nclusters,
                                     d1, B, d0, A, k, spacesize, spacedim,
                                     nagents, t.end, dt, timestamp))




  macro <- read.table('clusters_macro.csv', head=TRUE, sep=",")
  
  clu <- merge(params, macro, by=c("simname","simcount", "run"))
  


  clu$simname <- as.character(clu$simname)
  clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
  clu$simname <- as.factor(clu$simname)
  clu$simcount <- as.factor(clu$simcount)
  clu$t <- as.factor(clu$t)

  cl <- clu[clu$t == 2000,]

  return(clu)
}

theme_white <- function() {
  theme_update(panel.background = element_blank())
}



myLabeller <- function(var, value){
  value <- as.character(value)
  if (var == "R") {
    value[value== 0.03] <- "Small radius (R = 0.03)"
    value[value== 0.3] <- "Large Radius (R = 0.3)"
  } 
  return(value)
}

## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

## R ##

cl <- loadData(DUMPDIR, 'final_R/')

theme_set(theme_bw(base_size = 18))
theme_white()
title <- 'Cluster counts by radius of influence'
p <- ggplot(cl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + xlab('Radius of Influence') + ylab('Cluster counts')
p <- p + ggtitle(title)


p.mini <- ggplot(cl[cl$R <= 0.1,], aes(R, count))
p.mini <- p.mini + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p.mini <- p.mini + theme(legend.position = "none") 
p.mini

vp <- viewport(width = 0.4, height = 0.4,
               x = 0.6, y = 0.7)

full <- function() {
     print(p)
     theme_set(theme_bw(base_size = 8))
     theme_white()
     print(p.mini, vp = vp)
     theme_set(theme_bw(base_size = 18))
     theme_white()
}

jpeg(filename = paste0(IMGPATH, "scan_R.jpg"))
full()
dev.off()

## ALPHA ##

cl <- loadData(DUMPDIR, 'final_alpha/')

title <- 'Cluster counts by strength of influence'
p <- ggplot(cl, aes(alpha, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Strength of influence') + ylab('Cluster counts')
p <- p + ggtitle(title)

ggsave(filename = paste0(IMGPATH, "scan_alpha.jpg"), plot = p)
       
## VSCALING ##

cl <- loadData(DUMPDIR, 'final_vscaling/')

title <- 'Cluster counts by initial velocity'
p <- ggplot(cl, aes(init.vscaling, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Cluster counts')
p <- p + ggtitle(title)

ggsave(filename = paste0(IMGPATH, "scan_vscaling.jpg"), plot = p)

title <- 'Distance from truth by initial velocity'
p <- ggplot(cl[cl$t == 2000,], aes(init.vscaling, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Cluster counts')
p <- p + ggtitle(title)
p

ggsave(filename = paste0(IMGPATH, "scan_vscaling_fromtruth.jpg"), plot = p)

cl$vbr <- cut(cl$init.vscaling,  breaks=c(0,0.5,1, 1.5,2,10))
cl$vbr <- as.factor(cl$vbr)

title <- 'Cluster counts and distance from truth'
p <- ggplot(cl[cl$t == 2000 & cl$init.vscaling > 8,], aes(count, fromtruth.avg))
p <- p + geom_jitter(aes(size=init.vscaling, color=as.factor(R)), alpha=0.5)
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
p

ggsave(filename = paste0(IMGPATH, "scatter_count_fromtruth.jpg"), plot = p)

SCATTERPATH <- paste0(IMGPATH, 'scatter_count_fromtruth/')
idx = 1;
for (t in unique(cl$t)) {
  data <- cl[cl$t == t,]
  title <- paste0('Cluster counts and distance from truth - T', t) 
  p <- ggplot(data, aes(count, fromtruth.avg))
  p <- p + geom_jitter(aes(size=init.vscaling, color=as.factor(R)), alpha=0.5)
  p <- p + xlab('Number of clusters') + ylab('Distance from truth')
  p <- p + xlim(0,20) + ylim(0,0.5)
  p <- p + ggtitle(title)
  ggsave(filename = paste0(SCATTERPATH, "img_",sprintf("%04d",idx),".jpg"), plot = p)
  idx = idx + 1
}

system(paste0('ffmpeg -qscale 1 -r 1 -b 9600 -y -i ', SCATTERPATH, 'img_%04d.jpg ', SCATTERPATH, 'scatter_count_fromtruth.avi'))


## NOISES ##

cl <- loadData(DUMPDIR, 'final_noises/')

title <- 'Cluster counts by individualization noise'
p <- ggplot(cl, aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Individualization noise') + ylab('Cluster counts')
p <- p + ggtitle(title)

ggsave(filename = paste0(IMGPATH, "scan_noises.jpg"), plot = p)
