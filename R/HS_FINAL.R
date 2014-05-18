# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR 
DUMPDIR = "/mnt/tmp/dump/NAVNP/"

# Linear
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon/"
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v/"

DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'

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


loadData <- function(DUMPDIR, DIR, TWO_THOUSANDS = 0) {

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


  if (TWO_THOUSANDS) {
    macro <- macro[macro$t == 2000,]
  }
  
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

myThemeMod <- theme(legend.position = "none",
                    axis.title.x = element_text(vjust=-1, size=24),
                    axis.title.y = element_text(vjust=-0.1, size=24),
                    plot.margin=unit(c(10,10,10,10),"mm"),                    
                    plot.title = element_text(vjust=3, size=24,face="bold")
                    )

limits <- aes(ymax = count + se, ymin = count - se)
limitsCI <- aes(ymax = count + ci, ymin = count - ci)

theme_set(theme_bw(base_size = 30))
theme_white()

XINTERCEPT <- 0.15

## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

## R ##

cl <- loadData(DUMPDIR, 'final_R/')

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("R"), na.rm=TRUE)


# Breaks for background rectangles
# rects <- data.frame(xstart = c(-Inf, XINTERCEPT), xend = c(XINTERCEPT, Inf), col = letters[1:2])

# color <- "#FFEDD1"

title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryCl, aes(R, count))
# p <- p + geom_rect(aes(xmin = -Inf, xmax = XINTERCEPT, ymin = -Inf, ymax = Inf), fill = "#E8FCFF")
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
#p <- p + geom_point()
#p <- p + geom_line()
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
#p <- p + annotate("text", x = 0.05, y = 30, label = "Clusters", size=8)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
#p <- p + geom_segment(aes(x = 0, xend = 0.1, y = 29.4, yend = 29.4))
#p <- p + geom_segment(aes(x = 0.475, xend = 0.625, y = 29.4, yend = 29.4))
p <- p + xlab('Radius of Influence') + ylab('Number of Clusters')
p <- p  + myThemeMod # + ggtitle(title)
p

ggsave(filename = paste0(IMGPATH, "scan_R3.svg"),
       plot = p, width=10, height=5, dpi=300)


p.mini <- ggplot(summaryCl[summaryCl$R <= 0.1,], aes(R, count))
p.mini <- p.mini + geom_bar(stat = "identity", position="dodge", aes(fill=count))
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

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("alpha", "R", "tbr"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + ylim(0,29)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab('Strength of Social Influence') + ylab('Number of Clusters')
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "scan_alpha.svg"),
       plot = p, width=10, height=5, dpi=300)


# t = 20000
summaryCl <- summarySE(cl[cl$t == 20000 & cl$R == 0.03,], c("count"), c("alpha", "R"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
# p <- p + facet_grid(~ R, labeller = myLabeller)
# p <- p + ylim(0,28)
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab(xlabText) + ylab('Cluster counts')
labelText <- "Results after 20.000 iterations\nInteraction Radius = 0.03"
p <- p + annotate("text", x = 0.74, y = 13.5, label = labelText, size = 8)
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "scan_alpha_20000.jpg"),
       plot = p, width=10, height=5, dpi=300)
       
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

gsummaryCl2 <- rename(summaryCl2, c("se" = "se.fromtruth.avg"))
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

title <- 'Cluster counts vs Angular noise and \nPosition noise'
p <- ggplot(summaryCl, aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Angular noise') + ylab('Cluster counts')
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05", "0.075", "0.1"))
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "scan_noises.svg"),
       plot = p, width=10, height=10, dpi=300)


## TAU ##

cl <- loadData(DUMPDIR, 'final_tau/')

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("tau", "R", "init.vscaling"), na.rm=TRUE)

title <- 'Cluster counts by strength of the truth'
p <- ggplot(summaryCl[summaryCl$init.vscaling == 1,], aes((100 - tau), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=1))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + myThemeMod + theme(strip.background = element_blank())
p

# To be taken from TAU 20000
ggsave(filename = paste0(IMGPATH, "scan_tau_upto_10.jpg"),
       plot = p, width=10, height=5, dpi=300)


# Scatter fromtruth count by tau

summaryCl <- summarySE(cl, c("count"), c("tau", "R", "t", "init.vscaling"), na.rm=TRUE)
summaryCl2 <- summarySE(cl, c("fromtruth.avg"), c("tau", "R", "t", "init.vscaling"), na.rm=TRUE)
summaryCl2 <- rename(summaryCl2, c("se" = "se.fromtruth.avg"))
summaryCl2$sd <- NULL
summaryCl2$ci <- NULL
summaryCl2$N <- NULL
summaryCl3 <- merge(summaryCl, summaryCl2, by=c("tau","R","t","init.vscaling"))


title <- 'Cluster counts and distance from truth'
p <- ggplot(summaryCl3[summaryCl3$t == 2000 & summaryCl3$init.vscaling == 1,], aes(fromtruth.avg, count))
p <- p + geom_jitter(aes(size=(1/tau), color=as.factor(R)))
p <- p + xlab('Distance from truth') + ylab('Number of clusters')
p <- p + ggtitle(title)
#p <- p + facet_grid(init.vscaling ~ .)
p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Truth\nstrength")
p

# 2 POINTS, AVG
title <- 'Cluster counts and distance from truth'
p <- ggplot(summaryCl3[summaryCl3$t == 2000 & summaryCl3$tau == 1 & summaryCl3$init.vscaling == 1,], aes(count, fromtruth.avg))
p <- p + geom_point(aes(color=as.factor(R)))
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Truth\nstrength")
p


# POINTS ALL TAU = 1 INITV = 1
title <- 'Cluster counts and distance from truth'
p <- ggplot(cl[cl$t == 2000 & cl$init.vscaling == 1 & cl$tau == 1,], aes(count, fromtruth.avg))
p <- p + geom_jitter(aes(size=tau, color=as.factor(R))))
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Truth\nstrength")
p

# ALL POINTS. VINIT = 1
title <- 'Cluster counts and distance from truth'
p <- ggplot(cl[cl$t == 2000 & cl$init.vscaling == 1,], aes(count, fromtruth.avg))
p <- p + geom_jitter(aes(size=tau, color=as.factor(R)))
#p <- p + geom_jitter(aes(size=(1/tau), color=as.factor(R)))
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
#p <- p + facet_grid(init.vscaling ~ .)
p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Truth\nstrength")
p


title <- 'Cluster counts and distance from truth'
p <- ggplot(summaryCl3[summaryCl3$t == 2000,], aes(count, fromtruth.avg))
p <- p + geom_jitter(aes(color=(tau)))
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)
p <- p + facet_grid( R ~ init.vscaling, margins=T)
#p <- p + scale_color_hue(name="Influence\nradius size")
#p <- p + scale_size_continuous(name="Truth\nstrength")
p



title <- 'Cluster counts and distance from truth'
p <- ggplot(summaryCl3[summaryCl3$t == 2000,], aes(init.vscaling, tau))
p <- p + geom_jitter(aes(color=fromtruth.avg))
p <- p + xlab('Number of clusters') + ylab('Distance from truth')
p <- p + ggtitle(title)

3p <- p + scale_color_hue(name="Influence\nradius size")
p <- p + scale_size_continuous(name="Truth\nstrength")
p

9
p <- p + myThemeMod
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


## TAU 20000

cl <- loadData(DUMPDIR, 'final_tau_20000/')

summaryCl <- summarySE(cl[cl$t == 2000 & cl$tau <= 10,], c("count"), c("tau", "R"), na.rm=TRUE)

title <- 'Cluster counts by strength of the truth'
p <- ggplot(summaryCl, aes((11 - tau), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=1))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
xlabText <- expression(paste('Strength ', tau, ' of attraction by ground-truth'))
p <- p + xlab(xlabText) + ylab('Number of Clusters')
p <- p + scale_x_continuous(breaks=c(2,4,6,8,10))
p <- p + myThemeMod + theme(strip.background = element_blank())
p

# To be taken from TAU 20000
ggsave(filename = paste0(IMGPATH, "scan_tau_upto_10.svg"),
       plot = p, width=10, height=5, dpi=300)

#ggsave(filename = paste0(IMGPATH, "scan_tau_20000.jpg"),
#       plot = p, width=10, height=5, dpi=600)


## TAU 20000 NO BOUND

cl <- loadData(DUMPDIR, 'final_tau_20000_nobound/', 1)

summaryCl <- summarySE(cl[cl$t == 20000,], c("count"), c("tau", "R"), na.rm=TRUE)

title <- 'Cluster counts by strength of the truth'
p <- ggplot(summaryCl, aes((100 - tau), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=1))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "scan_tau_20000_nobound.jpg"),
       plot = p, width=10, height=5, dpi=600)


## ALPHA TAU

## get also tau = 1
clTau1 <- loadData(DUMPDIR, 'final_alpha/', 0)
clTau1 <- clTau1[clTau1$t == 2000,]

# might take long...
cl <- loadData(DUMPDIR, 'alpha_tau/', 1)

cl <- rbind(cl,clTau1)

mycl <- cl[cl$tau <= 20,]


summaryCl <- summarySE(mycl[mycl$R == 0.3,], c("count"), c("alpha"), na.rm=TRUE)

# Create plot
title <- paste0('Ccount: tau vs alpha - Big R')
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + ylim(0,40)
p <- p + xlab(xlabText) + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p


title <- "Tau vs Alpha vs Clusters Count"
p <- ggplot(cl, aes(alpha, tau, fill = count))
p <- p + geom_tile()
p <- p + ggtitle(title) + myThemeMod + theme(legend.position = "right")
p <- p + facet_grid(~R, labeller = myLabeller)
p <- p + scale_y_continuous(breaks=c(1,25,50,75,100))
p

ggsave(filename = paste0(IMGPATH, "heatmap_tau_alpha_ccount.jpg"),
       plot = p, width=14)


title <- "Tau vs Alpha vs Distance from Truth"
p <- ggplot(cl, aes(alpha, tau, fill = fromtruth.avg))
p <- p + geom_tile()
p <- p + ggtitle(title) + myThemeMod + theme(legend.position = "right")
p <- p + scale_fill_continuous(name="Dist.\nfrom\nTruth")
p <- p + facet_grid(~R, labeller = myLabeller)
p <- p + scale_y_continuous(breaks=c(1,25,50,75,100))
p

ggsave(filename = paste0(IMGPATH, "heatmap_tau_alpha_dist.jpg"),
       plot = p, width = 14)



mycl$taubr <- cut(mycl$tau,  breaks=c(0,5,10,15,20,25,50,100))
summaryCl <- summarySE(mycl, c("count"), c("alpha", "R", "taubr"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(aes(ymin = count -se, ymax = count + se))
p <- p + facet_grid(taubr ~ R, labeller = myLabeller)
#p <- p + ylim(0,28)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab(xlabText) + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p


ggsave(filename = paste0(IMGPATH, "MYCL_tau_alpha_brk_count.jpg"),
       plot = p, width = 14, height = 14)


#ggsave(filename = paste0(IMGPATH, "tau_alpha_brk_count.jpg"),
#       plot = p, width = 14, height = 14)

mycl2 <- mycl[mycl$tau == 1 | mycl$tau == 10 | mycl$tau == 20,]  

summaryFt <- summarySE(mycl2, c("fromtruth.avg"), c("alpha", "R", "tau"), na.rm=TRUE)

title <- 'Distance from Truth vs Strength of social influence'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg, group = tau, color=as.factor(tau)))
#p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(aes(ymin=fromtruth.avg - se, ymax = fromtruth.avg + se))
p <- p + facet_grid(. ~ R, labeller = myLabeller)
#p <- p + ylim(0,28)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "MYCL_tau_alpha_brk_dist.jpg"),
       plot = p, width = 14, height = 14)

# ggsave(filename = paste0(IMGPATH, "tau_alpha_brk_dist.jpg"),
#       plot = p, width = 14, height = 14)


taus <- 1:100
for (tau in taus) {
# Get data
summaryCl <- summarySE(cl[cl$tau == tau & cl$R == 0.3,], c("count"), c("alpha"), na.rm=TRUE)
# Create plot
title <- paste0(tau, ' - Ccount: tau vs alpha - Big R')
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
                                        #p <- p + geom_point()
                                        #p <- p + geom_line()
p <- p + geom_errorbar(limits)
                                        #p <- p + ylim(0,28)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + ylim(0,40)
p <- p + xlab(xlabText) + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
#p
ggsave(filename = paste0(IMGPATH, "tau_alpha_ccount/tau_alpha_brk_count_", tau, ".jpg"))
}


taus <- 1:100
for (tau in taus) {
# Get data
summaryFt <- summarySE(cl[cl$tau == tau & cl$R == 0.3,], c("fromtruth.avg"), c("alpha"), na.rm=TRUE)
# Create plot
title <- paste0(tau, ' - Dist: tau vs alpha - Big R')
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(aes(ymin = fromtruth.avg - se, ymax = fromtruth.avg + se))
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + ylim(0,0.5)
p <- p + xlab(xlabText) + ylab('Avg. Dist from Truth')
p <- p + ggtitle(title) + myThemeMod
#p
ggsave(filename = paste0(IMGPATH, "tau_alpha_dist/big_r/tau_alpha_dist_", tau, ".jpg"))
}

taus <- 1:100
for (tau in taus) {
# Get data
summaryFt <- summarySE(cl[cl$tau == tau & cl$R == 0.03,], c("fromtruth.avg"), c("alpha"), na.rm=TRUE)
# Create plot
title <- paste0(tau, ' - Dist: tau vs alpha - Small R')
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(aes(ymin = fromtruth.avg - se, ymax = fromtruth.avg + se))
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + ylim(0,0.55)
p <- p + xlab(xlabText) + ylab('Avg. Dist from Truth')
p <- p + ggtitle(title) + myThemeMod
#p
ggsave(filename = paste0(IMGPATH, "tau_alpha_dist/small_r/tau_alpha_dist_", tau, ".jpg"))
}


tau = 100
summaryCl <- summarySE(cl[cl$tau == tau & cl$R == 0.3,], c("count"), c("alpha"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
#p <- p + geom_point()
#p <- p + geom_line()
p <- p + geom_errorbar(limits)
#p <- p + ylim(0,28)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab(xlabText) + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p


tau = 100
summaryFt <- summarySE(cl[cl$tau == tau & cl$R == 0.03,], c("fromtruth.avg"), c("alpha"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(aes(ymin=fromtruth.avg - se, ymax = fromtruth.avg + se))
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab(xlabText) + ylab('Distance from Truth')
p <- p + ggtitle(title) + myThemeMod
p


# R TAU

## get also tau = 1
clTau1 <- loadData(DUMPDIR, 'final_R/', 0)
clTau1 <- clTau1[clTau1$t == 2000,]

# might take long...
cl <- loadData(DUMPDIR, 'R_tau/', 1)

cl <- rbind(cl,clTau1)

summaryCl <- summarySE(cl, c("count"), c("R"), na.rm=TRUE)

# Create plot
title <- paste0('Ccount: tau vs R')
p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + ylim(0,40)
p <- p + xlab("R") + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p




taus <- 1:100
for (tau in taus) {
# Get data
summaryCl <- summarySE(cl[cl$tau == tau,], c("count"), c("R"), na.rm=TRUE)
# Create plot
title <- paste0(tau, ' - Ccount: tau vs R')
p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + ylim(0,40)
p <- p + xlab("R") + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
#p
ggsave(filename = paste0(IMGPATH, "tau_R/count/tau_R_count_", tau, ".jpg"))
}


taus <- 1:100
for (tau in taus) {
# Get data
summaryFt <- summarySE(cl[cl$tau == tau,], c("fromtruth.avg"), c("R"), na.rm=TRUE)
# Create plot
title <- paste0(tau, ' - Dist: tau vs R')
p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(aes(ymin = fromtruth.avg - se, ymax = fromtruth.avg + se))
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + ylim(0,0.5)
p <- p + xlab("R") + ylab('Avg. Dist from Truth')
p <- p + ggtitle(title) + myThemeMod
p
ggsave(filename = paste0(IMGPATH, "tau_R/dist/tau_R_dist_", tau, ".jpg"))
}


cl$rbr <- cut(cl$R, c(0,0.05,0.1,0.15,0.2,0.3,1))

cl$taubr <- cut(cl$tau,  breaks=c(0,5,10,15,20,25,50,100))
summaryCl <- summarySE(cl, c("count"), c("R", "taubr"), na.rm=TRUE)

title <- paste0('Ccount: tau vs R')
p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + ylim(0,40)
p <- p + xlab("R") + ylab('Cluster counts')
p <- p + facet_grid(taubr ~ .)
p <- p + ggtitle(title) + myThemeMod
p

summaryFt <- summarySE(cl, c("fromtruth.avg"), c("R", "taubr"), na.rm=TRUE)

title <- paste0('Dist: tau vs R')
p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(aes(ymin = fromtruth.avg - se, ymax = fromtruth.avg + se))
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + ylim(0,0.5)
p <- p + xlab("R") + ylab('Dist')
p <- p + facet_grid(taubr ~ .)
p <- p + ggtitle(title) + myThemeMod
p


## ALL in one


cl <- loadData(DUMPDIR, 'final_R/')
clall <- cl[cl$t == 2000,]

cl <- loadData(DUMPDIR, 'final_alpha/')
clall <- rbind(clall, cl[cl$t == 2000,])

cl <- loadData(DUMPDIR, 'final_tau/')
clall2 <- rbind(clall, cl[cl$t == 2000,])

cl <- loadData(DUMPDIR, 'final_noises/')
clall <- rbind(clall, cl[cl$t == 2000,])



p <- ggplot(clall2[clall2$tau == 1,])
p <- p + geom_jitter(aes(fromtruth.avg, count))
p <- p + geom_smooth(aes(fromtruth.avg, count), method="lm") #glm, gam, loess,rlm
#p <- p + facet_grid(~tau)
#p <- p + facet_grid(~init.vscaling)
#p <- p + geom_smooth(data = cl[cl$t == 500,], aes(fromtruth.avg, count, color="red"))
#p <- p + geom_smooth(data = cl[cl$t == 1000,], aes(fromtruth.avg, count, color="blue"))
#p <- p + geom_smooth(data = cl[cl$t == 1500,], aes(fromtruth.avg, count, color="green"))
#p <- p + geom_smooth(data = cl[cl$t == 2000,], aes(fromtruth.avg, count, color="yellow"))
p


p <- ggplot(clall)
p <- p + geom_jitter(aes(fromtruth.avg, count))
p <- p + geom_smooth(aes(fromtruth.avg, count), method="loess") #glm, gam, loess,rlm
p <- p + facet_grid(~R)
#p <- p + facet_grid(~init.vscaling)
#p <- p + geom_smooth(data = cl[cl$t == 500,], aes(fromtruth.avg, count, color="red"))
#p <- p + geom_smooth(data = cl[cl$t == 1000,], aes(fromtruth.avg, count, color="blue"))
#p <- p + geom_smooth(data = cl[cl$t == 1500,], aes(fromtruth.avg, count, color="green"))
#p <- p + geom_smooth(data = cl[cl$t == 2000,], aes(fromtruth.avg, count, color="yellow"))
p




p <- ggplot(cl[cl$t == 2000,])
p <- p + geom_jitter(aes(fromtruth.avg, count, color = tau))
p <- p + geom_smooth(aes(fromtruth.avg, count), method="loess") #glm, gam, loess,rlm
p <- p + facet_grid(R~init.vscaling)
p
