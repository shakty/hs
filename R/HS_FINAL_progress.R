# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'
#DUMPDIR <- '/home/stefano/HS/'

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
  macro <- subset(macro, select=-c(move.avg,
                                   speed.avg,
                                   move.sd,
                                   speed.sd,
                                   fromtruth.avg,
                                   fromtruth.sd))
 

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

myFormatter <- function() {
  function(x) format(x, nsmall = 1, scientific = FALSE)
}

myThemeMod <- theme(legend.position = "none",
                    axis.title.x = element_text(vjust=-1, size=24),
                    axis.title.y = element_text(vjust=-0.1, size=24),
                    plot.margin=unit(c(10,10,10,10),"mm"),                    
                    plot.title = element_text(vjust=3, size=24,face="bold")
                    )

limits <- aes(ymax = count + se, ymin = count - se)
limits <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)
limitsFt <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)

theme_set(theme_bw(base_size = 30))
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
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Distance from Truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod # + ggtitle(title)
p

ggsave(filename = paste0(IMGPATH, "progress_scan_R2b.svg"),
       plot = p, width=10, height=5, dpi=300)

p.mini <- ggplot(summaryFt[summaryFt$R <= 0.1,], aes(R, fromtruth.avg))
p.mini <- p.mini + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p.mini <- p.mini + geom_errorbar(limits) + ylab('Distance from truth')
p.mini <- p.mini + myThemeMod
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


## R (truth-aside) ##

cl <- loadData(DUMPDIR, 'truth_aside_R_alpha/')

summaryFt <- summarySE(cl[cl$t == 21,], c("fromtruth.avg"), c("R"), na.rm=TRUE)

title <- 'Distance from truth by radius of influence'
p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "progress_scan_R2.jpg"),
       plot = p, width=10, height=5, dpi=600)

p.mini <- ggplot(summaryFt[summaryFt$R <= 0.1,], aes(R, fromtruth.avg))
p.mini <- p.mini + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p.mini <- p.mini + geom_errorbar(limits) + ylab('Distance from truth')
p.mini <- p.mini + myThemeMod
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

# & cl$alpha > 0.8
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("alpha", "R", "tbr"), na.rm=TRUE)

title <- 'Distance from truth vs Strength of social influence'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(~ R, labeller = myLabeller)
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab("Strength of social influence") + ylab('Avg. Distance from Truth')
p <- p + ylim(0,0.3)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "scan_alpha_dist.svg"),
       plot = p, width=10, height=5, dpi=300)

# t = 20000

summaryFt <- summarySE(cl[cl$t == 20000 & cl$R == 0.03,], c("fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)

title <- 'Distance from truth vs Strength of social influence'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab(xlabText) + ylab('Distance from truth')
labelText <- "Results after 20.000 iterations\nInteraction Radius = 0.03"
p <- p + annotate("text", x = 0.74, y = 0.155, label = labelText, size = 8)
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "progress_scan_alpha_20000.jpg"),
       plot = p, width=10, height=5, dpi=300)


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

title <- 'Distance from truth vs Angular noise and \nPosition noise'
p <- ggplot(summaryFt, aes(sigma, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Angular noise') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
#p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05", "0.075", "0.1"))
#p <- p + scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15))
p <- p + ggtitle(title) + myThemeMod + theme(panel.margin = unit(c(5),"mm"))
#p <- p + theme(axis.text.y = element_text(size=18))
p

ggsave(filename = paste0(IMGPATH, "progress_scan_noises.jpg"),
       plot = p, width=10, height=10, dpi=300)

## TAU ##

cl <- loadData(DUMPDIR, 'final_tau/')

summaryFt <- summarySE(cl[cl$t == 2000 & cl$tau <= 10,], c("fromtruth.avg"), c("tau", "R","init.vscaling"), na.rm=TRUE)

#options(warn = 2, error = recover)

vline_frame <- data.frame(intercept=89, R = 0.3)
ann_text <- data.frame(tau = 75, fromtruth.avg = 0.5, R = 0.3)
ann_rect <- data.frame(tau = 75, fromtruth.avg = 0.2, R = 0.3)
ann_arrow <- data.frame(x = 48, y = 0.48, R = 0.3)

title <- 'Distance from truth vs Strength of the truth\'s signal'
p <- ggplot(summaryFt[summaryFt$init.vscaling == 1,], aes((100 - tau), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=1))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(.~ R, labeller = myLabeller)
#p <- p + geom_vline(data = vline_frame, aes(xintercept = intercept), colour="red", linetype = "longdash", size = 1)
#p <- p + geom_text(data = ann_text, label = "Convergence far away\n from ground-truth", size=8)
#p <- p + geom_rect(data=ann_rect, aes(xmin = tau,
#                     xmax = tau + 20, ymin = fromtruth.avg,
#                     ymax = fromtruth.avg + 0.15), alpha = .2)
#p <- p + geom_segment(data=ann_arrow, aes(x = x, y = y, xend = x + 34, yend = y - 0.14), arrow = arrow())
p <- p + xlab('Truth strength in percentage') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "progress_scan_tau_upto_10.jpg"),
       plot = p, width=10, height=5, dpi=300)

ggsave(filename = paste0(IMGPATH, "progress_scan_tau.jpg"),
       plot = p, width=10, height=5, dpi=300)

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


## TAU 20000

cl <- loadData(DUMPDIR, 'final_tau_20000/')

summaryFt <- summarySE(cl[cl$t == 2000 & cl$tau <= 10,], c("fromtruth.avg"), c("tau", "R"), na.rm=TRUE)

#options(warn = 2, error = recover)

title <- 'Distance from truth by strength of the truth'
p <- ggplot(summaryFt, aes((11 - tau), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=1))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(.~ R, labeller = myLabeller)
xlabText <- expression(paste('Strength ', tau, ' of Attraction by Ground-Truth'))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_x_continuous(breaks=c(2,4,6,8,10))
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "progress_scan_tau_upto_10.svg"),
       plot = p, width=10, height=5, dpi=300)

#ggsave(filename = paste0(IMGPATH, "progress_scan_tau_20000.jpg"),
#       plot = p, width=10, height=5, dpi=600)

summaryCl <- summarySE(cl[cl$t == 20000,], c("coverage.cum"), c("tau", "R"), na.rm=TRUE)

title <- 'Cohesiveness biggest cluster by strength of the truth'
p <- ggplot(summaryCl, aes((100 - tau), coverage.cum))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=coverage.cum, width=1))
p <- p + geom_errorbar(aes(ymax = coverage.cum + se, ymin = coverage.cum - se))
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p

summaryCl <- summarySE(cl[cl$t == 20000,], c("pdist.sd"), c("tau", "R"), na.rm=TRUE)

title <- 'Cohesiveness biggest cluster by strength of the truth'
p <- ggplot(summaryCl, aes((100 - tau), pdist.sd))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=pdist.sd, width=1))
p <- p + geom_errorbar(aes(ymax = pdist.sd + se, ymin = pdist.sd - se))
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p

summaryCl <- summarySE(cl[cl$t == 20000,], c("pdist.mean"), c("tau", "R"), na.rm=TRUE)

title <- 'Cohesiveness biggest cluster by strength of the truth'
p <- ggplot(summaryCl, aes((100 - tau), pdist.mean))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=pdist.mean, width=1))
p <- p + geom_errorbar(aes(ymax = pdist.mean + se, ymin = pdist.mean - se))
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod
p


## TAU 20000 NO BOUND


cl <- loadData(DUMPDIR, 'final_tau_20000_nobound/')

summaryFt <- summarySE(cl[cl$t == 20000,], c("fromtruth.avg"), c("tau", "R"), na.rm=TRUE)

#options(warn = 2, error = recover)

title <- 'Distance from truth by strength of the truth'
p <- ggplot(summaryFt, aes((100 - tau), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=1))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + ggtitle(title) + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "progress_scan_tau_20000_nobound.jpg"),
       plot = p, width=10, height=5, dpi=300)
