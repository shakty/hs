# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR VELOCITY
DUMPDIR = "/mnt/tmp/dump/RND_SEED/"

DIR = "attrLinear_nv_rndseed_rndseq_tm_Rleft_n100_fv0/"
DIR = "attrZero_nav_rndseed_rndseq_tm_Rleft_n100_fv0/"

DIR = "attrZero_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0/"

INTERACTIVE = FALSE
PATH = paste0(DUMPDIR,DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");

# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/img/"))
} 

######################################
## MACRO AVG_SPLIT Aggregated: Loading
######################################

params <- read.table('params_all.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)

macro <- read.table('clusters_macro_avg_split.csv', head=TRUE, sep=",")
macro$simname <- as.character(macro$simname)
macro$simname <- substr(macro$simname, nchar(macro$simname)-1, nchar(macro$simname))
macro$simname <- as.factor(macro$simname)
macro$simcount <- as.factor(macro$simcount)
macro$t <- as.factor(macro$t)

cl <- macro
#cl <- subset(macro, t %% 100 == 0)
#cl$t <- as.factor(cl$t)

# SIZE with SD
title = "Evolution of cluster size (+Std.Dev) by sigma"
p.size <- ggplot(cl, aes(t,group=simname))
p.size <- p.size + geom_area(aes(y = size.avg, colour=simname, group=simname))
p.size <- p.size + geom_area(aes(y = size.sd, colour=simname, group=simname, fill=simname))
p.size <- p.size + facet_grid(simname~.,margins=F)
p.size <- p.size + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
if (INTERACTIVE) {
  p.size
}
# SIZE with SE
limits <- aes(ymax = cl$size.avg + cl$size.se, ymin=cl$size.avg - cl$size.se)
title = "Evolution of cluster size (+Std.Err.) by sigma"
p.size.se <- ggplot(cl, aes(t,group=simname))
#p.size.se <- p.size.se + geom_line(aes(y = size.avg, colour=simname, group=simname))
#p.size.se <- p.size.se + geom_errorbar(limits, width=0.0002)
p.size.se <- p.size.se + geom_smooth(aes(y = size.avg, colour=simname, group=simname))
p.size.se <- p.size.se + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
#p.size.se <- p.size.se + scale_x_log10()
if (INTERACTIVE) {
  p.size.se
}

# CUM EXPLORATION + SD
title = "Cumulative exploration (+Std.Dev.) by sigma"
p.explo.cum <- ggplot(cl, aes(t))
p.explo.cum <- p.explo.cum + geom_area(aes(y = coverage.cum.avg, colour=simname, group=simname))
p.explo.cum <- p.explo.cum + geom_area(aes(y = coverage.cum.sd, colour=simname, group=simname, fill=simname))
p.explo.cum <- p.explo.cum + facet_grid(simname~.,margins=F)
p.explo.cum <- p.explo.cum + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
if (INTERACTIVE) {
  p.explo.cum
}

# CUM EXPLORATION + SE
limits <- aes(ymax = cl$coverage.cum.avg + cl$coverage.cum.se, ymin=cl$coverage.cum.avg - cl$coverage.cum.se)
title = "Evolution of cumulative exploration (+Std.Err.) by sigma"
p.explo.cum.se <- ggplot(cl, aes(t,group=simname))
#p.explo.cum.se <- p.explo.cum.se + geom_line(aes(y = coverage.cum.avg, colour=simname, group=simname))
#p.explo.cum.se <- p.explo.cum.se + geom_errorbar(limits, width=0.2)
p.explo.cum.se <- p.explo.cum.se + geom_point(aes(y = coverage.cum.avg, colour=simname, group=simname))
p.explo.cum.se <- p.explo.cum.se + ggtitle(title) + xlab("Rounds") + ylab("% of space")
if (INTERACTIVE) {
  p.explo.cum.se
}

# POINT EXPLORATION + SD
title = "Point exploration (+Std.Dev.) by sigma"
p.explo <- ggplot(cl, aes(t))
p.explo <- p.explo + geom_area(aes(y = coverage.avg, colour=simname, group=simname))
p.explo <- p.explo + geom_area(aes(y = coverage.sd, colour=simname, group=simname, fill=simname))
p.explo <- p.explo + facet_grid(simname~.,margins=F)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
if (INTERACTIVE) {
  p.explo
}

# POINT EXPLORATION + SE
limits <- aes(ymax = cl$coverage.avg + cl$coverage.se, ymin=cl$coverage.avg - cl$coverage.se)
title = "Evolution of point exploration (+Std.Err.) by sigma"
p.explo.se <- ggplot(cl, aes(t,group=simname))
#p.explo.se <- p.explo.se + geom_line(aes(y = coverage.avg, colour=simname, group=simname))
#p.explo.se <- p.explo.se + geom_errorbar(limits, width=0.2)
#p.explo.se <- p.explo.se + geom_point(aes(y = coverage.avg, colour=simname, group=simname), alpha=0.2)
p.explo.se <- p.explo.se + geom_smooth(aes(y = coverage.avg, colour=simname, group=simname))
p.explo.se <- p.explo.se + ggtitle(title) + xlab("Rounds") + ylab("% of space")
if (INTERACTIVE) {
  p.explo.se
}

# SPEED + SD
title = "Mean and std. agents speed"
p.speed <- ggplot(cl, aes(t))
p.speed <- p.speed + geom_area(aes(y = speed.avg, colour=simname, group=simname))
p.speed <- p.speed + geom_area(aes(y = speed.sd, colour=simname, group=simname, fill=simname))
p.speed <- p.speed + facet_grid(simname~.,margins=F)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed
}

# SPEED + SE
limits <- aes(ymax = cl$speed.avg + cl$speed.se, ymin=cl$speed.avg - cl$speed.se)
title = "Evolution of speed (+Std.Err.) by sigma"
p.speed.se <- ggplot(cl, aes(t,group=simname))
p.speed.se <- p.speed.se + geom_smooth(aes(y = speed.avg, colour=simname, group=simname))
#p.speed.se <- p.speed.se + geom_line(aes(y = speed.avg, colour=simname, group=simname))
#p.speed.se <- p.speed.se + geom_errorbar(limits, width=0.2)
p.speed.se <- p.speed.se + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed.se
}

#MOVEMENTS + SD
title = "Evolution of movements (+Std.Dev.) by sigma"
p.move <- ggplot(cl, aes(t))
p.move <- p.move + geom_area(aes(y = move.avg, colour=simname, group=simname))
#p.move <- p.move + geom_area(aes(y = move.sd, colour=simname, group=simname, fill=simname))
p.move <- p.move + facet_grid(simname~.,margins=F)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Move")
if (INTERACTIVE) {
  p.move
}

# MOVEMENTS + SE
limits <- aes(ymax = cl$move.avg + cl$move.se, ymin=cl$move.avg - cl$move.se)
title = "Evolution of movements (+Std.Err.) by sigma"
p.move.se <- ggplot(cl, aes(t,group=simname))
p.move.se <- p.move.se + geom_point(aes(y = move.avg, colour=simname, group=simname))
#p.move.se <- p.move.se + geom_line(aes(y = move.avg, colour=simname, group=simname))
#p.move.se <- p.move.se + geom_errorbar(limits, width=0.2)
p.move.se <- p.move.se + ggtitle(title) + xlab("Rounds") + ylab("Move")
if (INTERACTIVE) {
  p.move.se
}

# FROM TRUTH + SD
title = "Evolution of distance from truth (+Std.Dev.) by sigma"
p.ft <- ggplot(cl, aes(t))
p.ft <- p.ft + geom_area(aes(y = fromtruth.avg, colour=simname, group=simname))
p.ft <- p.ft + geom_area(aes(y = fromtruth.sd, colour=simname, group=simname, fill=simname))
p.ft <- p.ft + facet_grid(simname~.,margins=F)
p.ft <- p.ft + ggtitle(title) + xlab("Rounds") + ylab("Distance from truth")
if (INTERACTIVE) {
  p.ft
}

# FROM TRUTH + SE
limits <- aes(ymax = cl$fromtruth.avg + cl$fromtruth.se, ymin=cl$fromtruth.avg - cl$fromtruth.se)
title = "Evolution of distance from truth (+Std.Err.) by sigma"
p.ft.se <- ggplot(cl, aes(t,group=simname))
p.ft.se <- p.ft.se + geom_point(aes(y = fromtruth.avg, colour=simname, group=simname))
#p.ft.se <- p.ft.se + geom_line(aes(y = fromtruth.avg, colour=simname, group=simname))
#p.ft.se <- p.ft.se + geom_errorbar(limits, width=0.2)
p.ft.se <- p.ft.se + ggtitle(title) + xlab("Rounds") + ylab("Distance from truth")
if (INTERACTIVE) {
  p.ft.se
}

### PLOTTING ALL OF THEM TOGETHER!!
jpeg(paste0(IMGPATH, "t_avg_overview_split.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.size, p.size.se,
                  p.explo.cum, p.explo.cum.se,
                  p.explo, p.explo.se,
                  p.speed, p.speed.se,
                  p.move,p.move.se,
                  p.ft, p.ft.se,
                  ncol=4,
                  main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()

################################
## MACRO AVG Aggregated: Loading
################################

params <- read.table('params_all.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)

macro <- read.table('clusters_macro_all.csv', head=TRUE, sep=",")
macro$simname <- as.factor(macro$simname)
macro$simcount <- as.factor(macro$simcount)
macro$run <- as.factor(macro$run)
# Merging macro and params: clu
clu <- merge(params, macro, by=c("simname","simcount","run"))
# Factorising
for (n in names(clu[1:25])) {
  clu[, n] <- as.factor(clu[, n])      
}

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
      
data <- clu
maxT <- max(clu$t)
idx = 1;
for (t in unique(clu$t)) {
  data <- clu[clu$t == t,]
  # From truth
  pt.ft <- heatmapFacets_fromtruth(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"ft/ft_",sprintf("%04d",idx),".jpg"),plot=pt.ft$p)
  # Count
  pt.c <- heatmapFacets_count(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"count/count_",sprintf("%04d",idx),".jpg"),plot=pt.c$p)
  # Size
  pt.s <- heatmapFacets_size(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"size/size_",sprintf("%04d",idx),".jpg"),plot=pt.s$p)
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
  # Creating an overview of the last frame of all indexes
  if (t == maxT) {
    jpeg(paste0(IMGPATH, "ht_overview.jpeg"), width=2048, height=768)
    p <- grid.arrange(pt.s$p, pt.c$p, pt.ft$p, pt.m$p, pt.sp$p, pt.ccov$p, pt.cov$p, ncol=4,
            main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
    dev.off()
  }
  idx = idx + 1
}


# Making a video out of the images
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/ft/ft_%04d.jpg img/ft/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/speed/speed_%04d.jpg img/speed/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/move/move_%04d.jpg img/move/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/size/size_%04d.jpg img/size/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/cumcov/cumcov_%04d.jpg img/cumcov/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/cov/cov_%04d.jpg img/cov/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/count/count_%04d.jpg img/count/movie.avi')


################
## MACRO AVG ALL
################

# if heatmap, needs to rename count to count.avg, maybe others too (???)

macro.avg.all <- read.table('clusters_macro_avg_all.csv', head=TRUE, sep=",")

# attach(macro.avg.all)

cl <- macro.avg.all

# COUNT AND SIZE
limits <- aes(ymax = cl$size.avg + cl$size.sd, ymin=cl$size.avg - cl$size.sd, colour="Std.Dev.")
limits2 <- aes(ymax = cl$size.avg + cl$size.se, ymin=cl$size.avg - cl$size.se, colour="Std.Err.")
title = "Evolution of cluster count and size"
p.count <- ggplot(cl, aes(t))
#p.count <- p.count + geom_line(aes(y = size.avg, colour="Avg"))
#p.count <- p.count + geom_point(aes(y = size.avg), alpha=.2, size=2)
#p.count <- p.count + geom_errorbar(limits, width=0.2)
#p.count <- p.count + geom_errorbar(limits2, width=0.2)
p.count <- p.count + geom_smooth(aes(y = size.avg, colour="size"), size=2)
p.count <- p.count + geom_smooth(aes(y = count.avg, colour="count"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.sd, colour="std. size"), size=2)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
if (INTERACTIVE) {
  p.count
}
title = "Evolution of speed (+Std.Err.) by sigma"
p.speed.se <- ggplot(cl, aes(t,group=simname))
p.speed.se <- p.speed.se + geom_line(aes(y = speed.avg, colour=simname, group=simname))
#p.speed.se <- p.speed.se + geom_errorbar(limits, width=0.2)
p.speed.se <- p.speed.se + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed.se
}

#COVERAGE
title = "Average and cumulative space exploration"
p.explo <- ggplot(cl, aes(t))
p.explo <- p.explo + geom_jitter(aes(y = coverage.avg), alpha=.2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.avg, colour="avg"), size=2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.cum.avg, colour="cum"), size=2)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
if (INTERACTIVE) {
  p.explo
}

# SPEED
title = "Mean and std. agents speed"
p.speed <- ggplot(cl, aes(t))
p.speed <- p.speed + geom_jitter(aes(y = speed.avg), alpha=.2)
p.speed <- p.speed + geom_smooth(aes(y = speed.avg, colour="avg"), size=2)
p.speed <- p.speed + geom_smooth(aes(y = speed.sd, colour="sd"), size=2)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
if (INTERACTIVE) {
  p.speed
}

#MOVEMENTS
title = "Mean and std. agents movements"
p.move <- ggplot(cl, aes(t))
p.move <- p.move + geom_jitter(aes(y = move.avg), alpha=.2)
p.move <- p.move + geom_smooth(aes(y = move.avg, colour="avg"), size=2)
p.move <- p.move + geom_smooth(aes(y = move.sd, colour="sd"), size=2)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Displacement")
if (INTERACTIVE) {
  p.move
}

# FROM TRUTH
title = "Mean and std. distance from truth"
p.truth <- ggplot(cl, aes(t))
p.truth <- p.truth + geom_jitter(aes(y = fromtruth.avg), alpha=.2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.avg, colour="avg"), size=2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.sd, colour="sd"), size=2)
p.truth <- p.truth + ggtitle(title) + xlab("Rounds") + ylab("Distance")
if (INTERACTIVE) {
  p.truth
}

#print(p.truth)
# grid.Arrange
jpeg(paste0(IMGPATH, "t_avg_overview.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.count, p.explo, p.speed, p.truth, p.move, ncol=3,
            main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()

