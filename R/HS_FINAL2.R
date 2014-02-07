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
  if (!file.exists(paste0(IMGPATH, "scatter_count_fromtruth_tau/"))) {
    dir.create(file.path(IMGPATH, "scatter_count_fromtruth_tau/"))
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
  #clu$t <- as.factor(clu$t)

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

limits <- aes(ymax = count + se, ymin = count - se)

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

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("R"), na.rm=TRUE)

title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_pointrange(limits)
p <- p + xlab('Radius of Influence') + ylab('Cluster counts')
p <- p + ggtitle(title) + theme(legend.position = "none") 


p.mini <- ggplot(summaryCl[summaryCl$R <= 0.1,], aes(R, count))
# p.mini <- p.mini + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p.mini <- p.mini + geom_bar(stat = "identity", position="dodge")
p.mini <- p.mini + geom_errorbar(limits)
p.mini <- p.mini + theme(legend.position = "none") 
p.mini

vp <- viewport(width = 0.4, height = 0.4,
               x = 0.75, y = 0.7)

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

cl$tbr <- cut(cl$t,  breaks=seq(0,20000,1000))
cl$tbr <- as.factor(cl$tbr)

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("alpha", "R", "tbr"), na.rm=TRUE)

title <- 'Cluster counts by social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + ylim(0,28)
p <- p + xlab('Strength of social influence') + ylab('Cluster counts')
p <- p + ggtitle(title) + theme(legend.position = "none")
p

ggsave(filename = paste0(IMGPATH, "scan_alpha.jpg"), plot = p)

# t = 20000
summaryCl <- summarySE(cl[cl$t == 20000,], c("count"), c("alpha", "R"), na.rm=TRUE)

title <- 'Cluster counts by social influence (20k iters)'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + ylim(0,28)
p <- p + xlab('Strength of social influence') + ylab('Cluster counts') + theme(legend.position = "none")
p <- p + ggtitle(title) # + annotate("text", x = 0.55, y = 0.25, label = "(After 20k iterations)")

ggsave(filename = paste0(IMGPATH, "scan_alpha_20000.jpg"), plot = p)
       
## VSCALING ##

cl <- loadData(DUMPDIR, 'final_vscaling/')
cl$vbr <- cut(cl$init.vscaling,  breaks=c(0,0.5,1, 1.5,2,10))
cl$vbr <- as.factor(cl$vbr)

summaryCl <- summarySE(cl, c("count"), c("init.vscaling", "R", "t"), na.rm=TRUE)

title <- 'Cluster counts by initial velocity'
p <- ggplot(summaryCl[summaryCl$t == 2000,], aes(init.vscaling, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Cluster counts')
p <- p + ggtitle(title) + theme(legend.position = "none")
p

ggsave(filename = paste0(IMGPATH, "scan_vscaling.jpg"), plot = p)

# From truth
summaryCl2 <- summarySE(cl, c("fromtruth.avg"), c("init.vscaling", "R", "t"), na.rm=TRUE)

title <- 'Distance from truth by initial velocity'
p <- ggplot(summaryCl2[summaryCl2$t == 2000,], aes(init.vscaling, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se))
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Distance from truth')
p <- p + ggtitle(title) + theme(legend.position = "none")
p

ggsave(filename = paste0(IMGPATH, "scan_vscaling_fromtruth.jpg"), plot = p)

# Scatter fromtruth count

summaryCl2 <- rename(summaryCl2, c("se" = "se.fromtruth.avg"))
summaryCl2$sd <- NULL
summaryCl2$ci <- NULL
summaryCl2$N <- NULL

summaryCl3 <- merge(summaryCl, summaryCl2, by=c("init.vscaling","R","t"))


title <- 'Cluster counts and distance from truth'
p <- ggplot(summaryCl3[summaryCl3$t == 2000,], aes(count, fromtruth.avg))
p <- p + geom_jitter(aes(size=init.vscaling, color=as.factor(R)), alpha=0.5)
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Initial velocity\nintervals")
p

ggsave(filename = paste0(IMGPATH, "scatter_count_fromtruth.jpg"), plot = p)

SCATTERPATH <- paste0(IMGPATH, 'scatter_count_fromtruth/')
idx = 1;
for (t in sort(unique(summaryCl3$t))) {
  data <- summaryCl3[summaryCl3$t == t,]
  title <- paste0('Cluster counts and distance from truth - T', t) 
  p <- ggplot(data, aes(count, fromtruth.avg))
  p <- p + geom_jitter(aes(size=init.vscaling, color=as.factor(R)), alpha=0.5)
  p <- p + xlab('Number of clusters') + ylab('Distance from truth')
  p <- p + xlim(0,20) + ylim(0,0.5)
  p <- p + ggtitle(title)
  p <- p + scale_color_hue(name="Influence\nradius size")
  p <- p + scale_size_continuous(name="Initial velocity\nintervals")
  ggsave(filename = paste0(SCATTERPATH, "img_",sprintf("%04d",idx),".jpg"), plot = p)
  idx = idx + 1
}

system(paste0('ffmpeg -qscale 1 -r 1 -b 9600 -y -i ', SCATTERPATH, 'img_%04d.jpg ', SCATTERPATH, 'scatter_count_fromtruth.avi'))

## NOISES ##

cl <- loadData(DUMPDIR, 'final_noises/')

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("sigma", "epsilon", "R"), na.rm=TRUE)

title <- 'Cluster counts by individualization and \nmeasurament noise'
p <- ggplot(summaryCl, aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Individualization noise') + ylab('Cluster counts')
p <- p + ggtitle(title) + theme(legend.position = "none")

ggsave(filename = paste0(IMGPATH, "scan_noises.jpg"), plot = p)


## TAU ##

cl <- loadData(DUMPDIR, 'final_tau/')

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("tau", "R"), na.rm=TRUE)

title <- 'Cluster counts by strength of the truth'
p <- ggplot(summaryCl, aes((100 - tau), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + ggtitle(title) + theme(legend.position = "none")

ggsave(filename = paste0(IMGPATH, "scan_tau.jpg"), plot = p)



# Scatter fromtruth count by tau

summaryCl <- summarySE(cl, c("count"), c("tau", "R", "t"), na.rm=TRUE)
summaryCl2 <- summarySE(cl, c("fromtruth.avg"), c("tau", "R", "t"), na.rm=TRUE)
summaryCl2 <- rename(summaryCl2, c("se" = "se.fromtruth.avg"))
summaryCl2$sd <- NULL
summaryCl2$ci <- NULL
summaryCl2$N <- NULL
summaryCl3 <- merge(summaryCl, summaryCl2, by=c("tau","R","t"))


title <- 'Cluster counts and distance from truth'
p <- ggplot(summaryCl3[summaryCl3$t == 2000,], aes(count, fromtruth.avg))
p <- p + geom_jitter(aes(size=(1/tau), color=as.factor(R)), alpha=0.8)
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Truth\nstrength")
p

ggsave(filename = paste0(IMGPATH, "scatter_count_fromtruth_tau.jpg"), plot = p)

SCATTERPATH <- paste0(IMGPATH, 'scatter_count_fromtruth_tau/')
idx = 1;
for (t in sort(unique(summaryCl3$t))) {
  data <- summaryCl3[summaryCl3$t == t,]
  title <- paste0('Cluster counts and distance from truth - T', t) 
  p <- ggplot(data, aes(count, fromtruth.avg))
  p <- p + geom_jitter(aes(size=(1/tau), color=as.factor(R)), alpha=0.8)
  p <- p + xlab('Number of clusters') + ylab('Distance from truth')
  p <- p + xlim(0,40) + ylim(0,0.55)
  p <- p + ggtitle(title)
  p <- p + scale_color_hue(name="Influence\nradius size")
  p <- p + scale_size_continuous(name="Truth\nstrength")
  ggsave(filename = paste0(SCATTERPATH, "img_",sprintf("%04d",idx),".jpg"), plot = p)
  idx = idx + 1
}

system(paste0('ffmpeg -qscale 1 -r 1 -b 9600 -y -i ', SCATTERPATH, 'img_%04d.jpg ', SCATTERPATH, 'scatter_count_fromtruth_tau.avi'))
