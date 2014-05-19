# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR 
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
limitsFt <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)
limitsCI <- aes(ymax = count + ci, ymin = count - ci)

theme_set(theme_bw(base_size = 30))
theme_white()

XINTERCEPT <- 0.15

## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/NOBOUND/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

## R ##
#######
cl <- loadData(DUMPDIR, 'nobound_R_tau/', 1)

# CL
summaryCl <- summarySE(cl[cl$t == 2000 & cl$tau == 2,], c("count"), c("R"), na.rm=TRUE)


title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Number of Clusters')
p <- p  + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "nobound_R_count.svg"),
       plot = p, width=10, height=5, dpi=300)


# FT
summaryFt <- summarySE(cl[cl$t == 2000 & cl$tau == 100,], c("fromtruth.avg"), c("R"), na.rm=TRUE)

title <- 'Distance from truth by radius of influence'
p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Distance from Truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod 
p

ggsave(filename = paste0(IMGPATH, "nobound_R_dist.svg"),
       plot = p, width=10, height=5, dpi=300)

## ALPHA ##
###########

cl <- loadData(DUMPDIR, 'final_alpha/')

cl$tbr <- cut(cl$t,  breaks=seq(0,20000,1000))

# CL
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

# FT
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

## NOISES ##
############

cl <- loadData(DUMPDIR, 'final_noises/')

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("sigma", "epsilon", "R"), na.rm=TRUE)

title <- 'Cluster counts vs Angular noise and \nPosition noise'
p <- ggplot(summaryCl, aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Angular Noise') + ylab('Number of Clusters')
#p <- p + scale_x_continuous(labels = c("0", "0.02", "0.04", "0.06", "0.08", "0.1"),
#                            breaks = seq(0,0.1,0.02))
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
p <- p + myThemeMod + theme(strip.background = element_blank())
p


# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[22]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0"))
grob[["grobs"]][[23]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.1"))
grob[["grobs"]][[24]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.2"))
grob[["grobs"]][[25]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.3"))
grob[["grobs"]][[26]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.4"))
grob[["grobs"]][[27]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.5"))

svg(filename = paste0(IMGPATH, "scan_noises.svg"),
     width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()


# FT
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("sigma", "epsilon", "R"), na.rm=TRUE)

title <- 'Distance from truth vs Angular noise and \nPosition noise'
p <- ggplot(summaryFt, aes(sigma, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Angular Noise') + ylab('Avg. Distance from Truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
#p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05", "0.075", "0.1"))
#p <- p + scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15))
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
p <- p + myThemeMod +  theme(strip.background = element_blank())
#p <- p + theme(axis.text.y = element_text(size=18))
p

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[22]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0"))
grob[["grobs"]][[23]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.1"))
grob[["grobs"]][[24]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.2"))
grob[["grobs"]][[25]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.3"))
grob[["grobs"]][[26]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.4"))
grob[["grobs"]][[27]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.5"))

svg(filename = paste0(IMGPATH, "scan_noises_progress.svg"),
     width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()

## TAU ##
#########

cl <- loadData(DUMPDIR, 'final_tau/')

# CL
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


# FT
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
