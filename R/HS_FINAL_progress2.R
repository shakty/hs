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
  if (!file.exists(paste0(IMGPATH, "scatter_count_fromtruth/"))) {
    dir.create(file.path(IMGPATH, "scatter_count_fromtruth/"))
  }
  if (!file.exists(paste0(IMGPATH, "tau_fromtruth/"))) {
    dir.create(file.path(IMGPATH, "tau_fromtruth/"))
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
  macro$fromtruth.avg <- NULL
  macro$fromtruth.sd <- NULL
  agents <- read.table('agents.csv', head=TRUE, sep=",")
  
  clu <- merge(params, macro, by=c("simname","simcount", "run"))
  clu <- merge(clu, agents,  by=c("simname","simcount", "run", "t"))

  clu$simname <- as.character(clu$simname)
  clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
  clu$simname <- as.factor(clu$simname)
  clu$simcount <- as.factor(clu$simcount)
  #clu$t <- as.factor(clu$t)

  # cl <- clu[clu$t == 2000,]

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

limits <- aes(ymax = count + se, ymin = count - se)
limits <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)

theme_set(theme_bw(base_size = 18))
theme_white()

## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

## R ##

cl <- loadData(DUMPDIR, 'final_R/')

summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("R"), na.rm=TRUE)

title <- 'Distance from truth by radius of influence'
p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limits)
p <- p + xlab('Radius of Influence') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none") 
p <- p + ggtitle(title)


p.mini <- ggplot(summaryFt[summaryFt$R <= 0.1,], aes(R, fromtruth.avg))
p.mini <- p.mini + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p.mini <- p.mini + geom_errorbar(limits)
p.mini <- p.mini + theme(legend.position = "none") 
p.mini

vp <- viewport(width = 0.4, height = 0.4,
               x = 0.7, y = 0.7)

full <- function() {
     print(p)
     theme_set(theme_bw(base_size = 8))
     theme_white()
     print(p.mini, vp = vp)
     theme_set(theme_bw(base_size = 18))
     theme_white()
}

jpeg(filename = paste0(IMGPATH, "progress_scan_R.jpg"))
full()
dev.off()

## ALPHA ##

cl <- loadData(DUMPDIR, 'final_alpha/')

cl$tbr <- cut(cl$t,  breaks=seq(0,20000,1000))
cl$tbr <- as.factor(cl$tbr)

summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("alpha", "R", "tbr"), na.rm=TRUE)

title <- 'Distance from truth by social influence'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Strength of social influence') + ylab('Distance from truth')
p <- p + ylim(0,0.3)
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none") 
p <- p + ggtitle(title)
p

ggsave(filename = paste0(IMGPATH, "progress_scan_alpha.jpg"), plot = p)

# t = 20000
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)
title <- 'Distance from truth by social influence (20k iters)'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
# p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
#p <- p + geom_errorbar(limits)
p <- p + geom_pointrange(aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se, group=as.factor(R), color=as.factor(R)), size=0.5)


data <- summarySE(cl[cl$t == 20000 | cl$t == 2000,], c("fromtruth.avg"), c("alpha", "R","t"), na.rm=TRUE)
p <- ggplot(data, aes((1 - alpha), fromtruth.avg))
p <- p + geom_point(aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se, group=as.factor(t), color=as.factor(t), shape=as.factor(R)))
p

p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Strength of social influence') + ylab('Distance from truth')
p <- p + ylim(0,0.3)
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none") 
p <- p + ggtitle(title) # + annotate("text", x = 0.55, y = 0.25, label = "(After 20k iterations)")
p

ggsave(filename = paste0(IMGPATH, "progress_scan_alpha_20000.jpg"), plot = p)
       
## VSCALING ##

cl <- loadData(DUMPDIR, 'final_vscaling/')
cl$vbr <- cut(cl$init.vscaling,  breaks=c(0,0.5,1, 1.5,2,10))
cl$vbr <- as.factor(cl$vbr)

summaryFt <- summarySE(cl, c("fromtruth.avg"), c("init.vscaling", "R", "t"), na.rm=TRUE)

title <- 'Distance from truth by initial velocity'
p <- ggplot(summaryFt[summaryCl$t == 2000,], aes(init.vscaling, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none") 
p <- p + ggtitle(title)
p

ggsave(filename = paste0(IMGPATH, "progress_scan_vscaling.jpg"), plot = p)

## NOISES ##

cl <- loadData(DUMPDIR, 'final_noises/')

summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("sigma", "epsilon", "R"), na.rm=TRUE)

title <- 'Distance from truth by individualization and \nmeasurament noise'
p <- ggplot(summaryFt, aes(sigma, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Individualization noise') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none")
# p <- p + scale_y_continuous(breaks=c(0,0.1,0.15))
p <- p + ggtitle(title)

ggsave(filename = paste0(IMGPATH, "progress_scan_noises.jpg"), plot = p)


## TAU ##

cl <- loadData(DUMPDIR, 'final_tau/')

summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("tau", "R"), na.rm=TRUE)

title <- 'Distance from truth by strength of the truth'
p <- ggplot(summaryFt, aes((100 - tau), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none")
p <- p + ggtitle(title)

ggsave(filename = paste0(IMGPATH, "progress_scan_tau.jpg"), plot = p)

summaryFt <- summarySE(cl, c("fromtruth.avg"), c("tau", "R", "t"), na.rm=TRUE)
SCATTERPATH <- paste0(IMGPATH, 'tau_fromtruth/')
idx = 1;
for (t in sort(unique(summaryFt$t))) {
  data <- summaryFt[summaryFt$t == t,]
  title <- paste0('Distance from truth by strength of the truth - ', t)
  p <- ggplot(data, aes((100 - tau), fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
  p <- p + geom_errorbar(limits)
  p <- p + ylim(0,0.55)
  p <- p + facet_grid(.~ R, labeller = myLabeller)
  p <- p + xlab('Truth strength in percentage') + ylab('Distance from truth')
  p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none")
  p <- p + ggtitle(title)
  ggsave(filename = paste0(SCATTERPATH, "img_",sprintf("%04d",idx),".jpg"), plot = p)
  idx = idx + 1
}

system(paste0('ffmpeg -qscale 1 -r 1 -b 9600 -y -i ', SCATTERPATH, 'img_%04d.jpg ', SCATTERPATH, 'tau_fromtruth.avi'))
