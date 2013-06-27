# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

DIR = "test_t-2013-6-4-12-14"

DIR = "attrExpo_nv_rndseq_tm_Rleft_reduced/"
DIR = "attrExpo_nv_rndseq_tm_Rleft/"

DIR = "attrExpo_nv_rndseq_tm_Rleft/attrExpo_nv_rndseq_tm_Rleft_s5/" 


DUMPDIR = "/opt/MATLAB_WORKSPACE/hs/dump/NEW/"
PATH = paste0(DUMPDIR,DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");

# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/img/"))
}

## MACRO Aggregated: Loading
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

## MACRO: Loading
params <- read.table('params.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)

macro <- read.table('clusters_macro.csv', head=TRUE, sep=",")
macro$simname <- as.factor(macro$simname)
macro$simcount <- as.factor(macro$simcount)
macro$run <- as.factor(macro$run)
# Merging macro and params: clu
clu <- merge(params, macro, by=c("simname","simcount","run"))
# Factorising
for (n in names(clu[1:23])) {
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
for (t in unique(clu$t)) {
  data <- clu[clu$t == t,]
  # From truth
  pt.ft <- heatmapFacets_fromtruth(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"ft/ft_",sprintf("%04d",t),".jpg"),plot=pt.ft$p)
  # Count
  pt.c <- heatmapFacets_count(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"count/count_",sprintf("%04d",t),".jpg"),plot=pt.c$p)
  # Size
  pt.s <- heatmapFacets_size(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"size/size_",sprintf("%04d",t),".jpg"),plot=pt.s$p)
  # Move
  pt.m <- heatmapFacets_move(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"move/move_",sprintf("%04d",t),".jpg"),plot=pt.m$p)
  # Speed
  pt.sp <- heatmapFacets_speed(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"speed/speed_",sprintf("%04d",t),".jpg"),plot=pt.sp$p)
  # Cum Cov
  pt.ccov <- heatmapFacets_cumcoverage(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"cumcov/cumcov_",sprintf("%04d",t),".jpg"),plot=pt.ccov$p)
   # Cov Instantaneous
  pt.cov <- heatmapFacets_coverage(v1,v2,v3, data, t=t)
  ggsave(filename=paste0(IMGPATH,"cov/cov_",sprintf("%04d",t),".jpg"),plot=pt.cov$p)
  # Creating an overview of the last frame of all indexes
  if (t == maxT) {
    jpeg(paste0(IMGPATH, "ht_overview.jpg"), width=2048, height=768)
    p <- grid.arrange(pt.s$p, pt.c$p, pt.ft$p, pt.m$p, pt.sp$p, pt.ccov$p, pt.cov$p, ncol=4,
            main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
    dev.off()
  } 
}


# Making a video out of the images
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/ft/ft_%04d.jpg img/ft/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/speed/speed_%04d.jpg img/speed/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/move/move_%04d.jpg img/move/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/size/size_%04d.jpg img/size/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/cumcov/cumcov_%04d.jpg img/cumcov/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/cov/cov_%04d.jpg img/cov/movie.avi')
system('ffmpeg -qscale 1 -r 1 -b 9600 -y -i img/count/count_%04d.jpg img/count/movie.avi')

## MICRO

micro <- read.table('clusters_micro.csv', head=TRUE, sep=",")

CUTOFF <- 1
START.AFTER.ROUND <- 1
micro.sub <- micro[clusters.micro$size > CUTOFF & clusters$t > START.AFTER.ROUND, ]

micro.sub <- micro

p <- ggplot(micro.sub, aes(size, speed))
p <- p + geom_point()
p

p <- ggplot(micro.sub, aes(size, fromtruth))
p <- p + geom_point()
p

p <- ggplot(micro.sub, aes(size, move))
p <- p + geom_point()
p

p <- ggplot(micro.sub, aes(speed, move))
p <- p + geom_point()
p

p <- ggplot(micro.sub, aes(speed, fromtruth))
p <- p + geom_point()
p




## MACRO AVG ALL

# if heatmap, needs to rename count to count.avg, maybe others too

macro.avg.all <- read.table('clusters_macro_avg_all.csv', head=TRUE, sep=",")

attach(macro.avg.all)

cl <- macro.avg.all
# COUNT AND SIZE
title = "Evolution of cluster count and size"
p.count <- ggplot(cl, aes(t))
p.count <- p.count + geom_jitter(aes(y = size.avg), alpha=.2)
p.count <- p.count + geom_smooth(aes(y = count.avg, colour="count"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.avg, colour="size"), size=2)
p.count <- p.count + geom_smooth(aes(y = size.sd, colour="std. size"), size=2)
p.count <- p.count + ggtitle(title) + xlab("Rounds") + ylab("Agents per cluster")
#p.count
#COVERAGE
title = "Average and cumulative space exploration"
p.explo <- ggplot(cl, aes(t))
p.explo <- p.explo + geom_jitter(aes(y = coverage.avg), alpha=.2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.avg, colour="avg"), size=2)
p.explo <- p.explo + geom_smooth(aes(y = coverage.cum.avg, colour="cum"), size=2)
p.explo <- p.explo + ggtitle(title) + xlab("Rounds") + ylab("Percentage")
#p.explo
# SPEED
title = "Mean and std. agents speed"
p.speed <- ggplot(cl, aes(t))
p.speed <- p.speed + geom_jitter(aes(y = speed.avg), alpha=.2)
p.speed <- p.speed + geom_smooth(aes(y = speed.avg, colour="avg"), size=2)
p.speed <- p.speed + geom_smooth(aes(y = speed.sd, colour="sd"), size=2)
p.speed <- p.speed + ggtitle(title) + xlab("Rounds") + ylab("Speed")
#p.speed
#MOVEMENTS
title = "Mean and std. agents movements"
p.move <- ggplot(cl, aes(t))
p.move <- p.move + geom_jitter(aes(y = move.avg), alpha=.2)
p.move <- p.move + geom_smooth(aes(y = move.avg, colour="avg"), size=2)
p.move <- p.move + geom_smooth(aes(y = move.sd, colour="sd"), size=2)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Displacement")
#p.move
# FROM TRUTH
title = "Mean and std. distance from truth"
p.truth <- ggplot(cl, aes(t))
p.truth <- p.truth + geom_jitter(aes(y = fromtruth.avg), alpha=.2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.avg, colour="avg"), size=2)
p.truth <- p.truth + geom_smooth(aes(y = fromtruth.sd, colour="sd"), size=2)
p.truth <- p.truth + ggtitle(title) + xlab("Rounds") + ylab("Distance")
#print(p.truth)
# grid.Arrange
jpeg(paste0(IMGPATH, "t_avg_overview.jpeg"), width=2048, height=1024)
p <- grid.arrange(p.count, p.explo, p.speed, p.truth, p.move, ncol=3,
            main=textGrob(DIR, gp=gpar(cex=1.5, fontface="bold")))
dev.off()




# START


v1 <- "R"
v2 <- "alpha"
v3 <- "sigma"
#data <- clu
data <- clu[clu$t == 21,]
paramsData <- params
heatmapFacets(v1,v2,v3, data)


OLDPATH = IMGPATH
for (S in unique(params$sigma)) {
  curDir <-  paste0("sigma_",S,"/")
  dir.create(file.path(OLDPATH, curDir), showWarnings = FALSE)
  IMGPATH <- paste0(OLDPATH, curDir)
  heatmap2by2Detail(v1,v2, data = clu[clu$sigma == S,], paramsData = params[params$sigma == S,])
}
IMGPATH = OLDPATH


#image(clu$R, clu$sigma, clu$fromtruth.avg)

# ALL PLOTS
#allPlots("R","sigma")

### Print Convergence by R
# We need to have sigma and t as numeric, not as factors

AA <- clu[clu$t == 21,]

selected <- AA[AA$sigma == 0, ]
p <- ggplot(selected, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p

ggsave(filename=paste0(PATH,"convergence_by_r_alpha=.6_full.jpg"), plot=p)

selected <- AA[AA$sigma == 0 & AA$R < 0.6 & AA$R > 0.3,]
p <- ggplot(selected, aes(R, fromtruth.avg, group=alpha, colour=alpha))
p <- p + geom_line()
p

ggsave(filename=paste0(PATH,"convergence_by_r_alpha=.6_zoom.jpg"), plot=p)




selected <- AA #[AA$sigma == 0.1,]
selected$sigma <- as.factor(selected$sigma)
p <- ggplot(selected, aes(R, fromtruth.avg, group=sigma, colour=sigma))
p <- p + geom_line() + facet_grid(sigma~.)
p

selected$sigma <- as.factor(selected$sigma)
p <- ggplot(selected, aes(R, fromtruth.avg, group=sigma, colour=sigma))
p <- p + geom_line()
p


### End
  
