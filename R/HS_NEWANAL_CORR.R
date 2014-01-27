# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR 
DUMPDIR = "/mnt/tmp/dump/NAVNP/"

# Linear
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon/"
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v/"

# Zero
#DIR = "attrZero_navnp_RClean_n100_fv0_s1_epsilon/"
#DIR = "attrZero_navnp_RClean_n100_fv0_s1_epsilon_v/"

# Hard
#DIR = "attrHard_navnp_RClean_n100_fv0_s1_epsilon/"
#DIR = "attrHard_navnp_RClean_n100_fv0_s1_epsilon_v/"

# Millean
#DIR = "attrMillean_navnp_RClean_n100_fv0_s1_epsilon/"
#DIR = "attrMillean_navnp_RClean_n100_fv0_s1_epsilon_v/"


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

cluagents <- merge(params, agents, by=c("simname","simcount"))

subcluagents <- subset(cluagents, select = c(simname, simcount,t,
                                    pdist.sd,
                                    fromtruth.avg))

clu <- merge(params, macro, by=c("simname","simcount"))
subclu <- subset(clu, select = c(simname, simcount, alpha, R, sigma, init.vscaling, t,
                        count, size.avg, bigc.pdist.mean,
                        size.max))

cluall <-  merge(subclu, subcluagents, by=c("simname","simcount","t"))



cluallall$sigma <- as.factor(cluall$sigma)
cluall$tbr <- cut(macro$t, b = 4)
                      
cluall$truth <- agents$fromtruth.avg

cluall$tbr <- cut(cluall$t, breaks=c(0,500,1000, 1500, 2002))
cluall$tbr <- as.factor(cluall$tbr)

cluall$Rbr <- cut(cluall$R,  breaks=c(0,0.04,0.07,1))
cluall$Rbr <- as.factor(cluall$Rbr)

cluall$alphabr <- cut(cluall$alpha, breaks=c(0,0.25,0.5,0.75,1))
cluall$alphabr <- as.factor(cluall$alphabr)

# FROM TRUTH
title = "Scatterplot avg cluster count ~ distance from truth by velocity and time"
p <- ggplot(cluall, aes(x = count, y = fromtruth.avg))
p <- p + geom_point(aes(color=as.factor(Rbr)), size=1.5, alpha=0.2)                         
p <- p + facet_grid(init.vscaling ~ tbr)
p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
p <- p + ggtitle(title)
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "new/velocity_ccount_fromtruth.jpg"), plot=p)

title = "Scatterplot max cluster size ~ distance from truth by velocity and time"
p <- ggplot(cluall, aes(x = size.max, y = fromtruth.avg))
p <- p + geom_point(aes(color=sigma), size=1.5, alpha=0.2)                         
p <- p + facet_grid(init.vscaling ~ tbr)
p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
p <- p + ggtitle(title)
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "new/velocity_maxsize_fromtruth.jpg"), plot=p)

# color = velocity
for (s in unique(cluall$sigma)) {
  title = paste0("Scatterplot max cluster size ~ distance from truth by alpha and R\nsigma = ", s)
  p <- ggplot(cluall[cluall$t == 2000 & cluall$sigma == s,], aes(x = size.max, y = fromtruth.avg))
  p <- p + geom_point(aes(color=as.factor(init.vscaling)), alpha=0.5)                         
  p <- p + facet_grid(Rbr ~ alphabr)
  p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
  p <- p + ggtitle(title)
  if (INTERACTIVE) {
    p
  }
  ggsave(filename=paste0(IMGPATH, paste0("new/velocity_maxsize_fromtruth_t2000_s", s, ".jpg")), plot=p)
}

# color = sigma, by velocity => ALPHA does not have an impact on the final distribution.
for (s in unique(cluall$init.vscaling)) {
  title = paste0("Scatterplot max cluster size ~ distance from truth by alpha and R\ninitial velocity = ",s)
  p <- ggplot(cluall[cluall$t == 2000 & cluall$init.vscaling == s,], aes(x = size.max, y = fromtruth.avg))
  p <- p + geom_point(aes(color=as.factor(sigma)), alpha=0.5)
  p <- p + facet_grid(Rbr ~ alphabr)
  p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
  p <- p + ggtitle(title)
  if (INTERACTIVE) {
    p
  }
  ggsave(filename=paste0(IMGPATH, paste0("new/maxsize_fromtruth_by_v", s, ".jpg")), plot=p)
}

# Summary of above, without alpha
title = paste0("Scatterplot max cluster size ~ distance from truth \n by initial velocity and R")
p <- ggplot(cluall[cluall$t == 2000, ], aes(x = size.max, y = fromtruth.avg))
p <- p + geom_point(aes(color=as.factor(sigma)), alpha=0.5)
p <- p + facet_grid(. ~ init.vscaling)
p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
p <- p + ggtitle(title)
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, paste0("new/maxsize_fromtruth_by_v_and_r.jpg")), plot=p)



fit <- lm(100*fromtruth.avg ~ size.max, data = cluall[cluall$t == 2000,])
summary(fit)


## Testing...
title = paste0("Scatterplot max cluster size ~ distance from truth \n by initial velocity and R")
p <- ggplot(cluall[cluall$t == 2000, ], aes(x = count, y = fromtruth.avg))
p <- p + geom_jitter(aes(group = as.factor(sigma), color=as.factor(sigma)))
p <- p + facet_grid(sigma ~ init.vscaling)
#p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
p <- p + ggtitle(title)
p


p <- ggplot(cluall[cluall$t == 2000 & cluall$count < 40, ], aes(x = count, y = fromtruth.avg))
p <- p + geom_smooth(method="lm")
p



p <- ggplot(cluall[cluall$t == 2000,], aes(x = init.vscaling, y = fromtruth.avg))
p <- p + geom_bar(aes(group = as.factor(sigma), color=as.factor(sigma)))
p


p <- ggplot(cluall, aes(x = size.avg, y = fromtruth.avg))
p <- p + geom_point(aes(color=as.factor(t)))
p <- p + facet_grid(. ~ init.vscaling)
p



title = "Distributions of size of the biggest cluster at the end of simulation by sigma, epsilon (outer), and alpha and R (inner)"

p <- ggplot(cluall[cluall$t == 2000,], aes(x=init.vscaling, y=fromtruth.avg))
p <- p + geom_density2d(aes(fill = fromtruth.avg), geom="polygon")
p

p <- p + scale_fill_continuous(low='lightblue',high='red')


p <- p + scale_x_discrete(breaks = seq(0, 1, 0.02))
p <- p + scale_y_discrete(breaks = seq(0, 1, 0.05))


p <- p + facet_grid(sigma ~ epsilon)
p <- p + ggtitle(title) + xlab("Epsilons") + ylab("Sigmas")




# Generate data library(reshape2) # for melt
volcano3d <- melt(volcano)
names(volcano3d) <- c("x", "y", "z")
# Basic plot
v <- ggplot(volcano3d, aes(x, y, z = z))
v + stat_contour() 




# Spinning 3d Scatterplot
library(rgl)

clulast <- cluall[cluall$t == 2000,]

plot3d(clulast$R, clulast$alpha, clulast$fromtruth.avg, type = "p")



data(volcano)
z <- 3 * volcano # Exaggerate the relief
x <- 10 * (1:nrow(z)) # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z)) # 10 meter spacing (E to W)
zlim <- range(y)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- terrain.colors(zlen) # height color lookup table
col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
open3d()

surface3d(x, y, z, color=col, back="lines")

C <- cluall[cluall$sigma == 0.01 & cluall$t == 2000 & cluall$init.vscaling == 1,]


open3d()
surface3d(cl, y, z, color=col, back="lines")


levelplot(fromtruth.avg ~ alpha * R, data=cluall[cluall$t == 2000,],
          shade = TRUE,
          # aspect = c(61/87, 0.4),
          # drape = TRUE, colorkey = TRUE,
          pretty = TRUE,
          strip = FALSE,
          #screen = list(x = 45, y = 90, z = 135)
          #light.source = c(10,0,10)
          )


p <- wireframe(fromtruth.avg ~ sigma * R, data=cluall[cluall$t == 2000,], shade = TRUE)
npanel <- c(4, 2)
rotx <- c(-50, -80)
rotz <- seq(30, 300, length = npanel[1]+1)
update(p[rep(1, prod(npanel))], layout = npanel,
    panel = function(..., screen) {
        panel.wireframe(..., screen = list(z = rotz[current.column()],
                                           x = rotx[current.row()]))
    })


# STD PDIST

title = "Scatterplot max cluster size ~ st.d. pair-wise distance between agents \n by velocity and time"
p <- ggplot(cluall, aes(x = size.max, y = pdist.sd))
p <- p + geom_point(aes(color=sigma), size=1.5, alpha=0.2)
p <- p + facet_grid(init.vscaling ~ tbr)
p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
p <- p + ggtitle(title)
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, "newpdist/velocity_maxsize_sdpdist.jpg"), plot=p)


# color = velocity
for (s in unique(cluall$sigma)) {
  title = paste0("Scatterplot max cluster size ~ st.d. pair-wise distance \n between agents by alpha and R. sigma = ", s)
  p <- ggplot(cluall[cluall$t == 2000 & cluall$sigma == s,], aes(x = size.max, y = pdist.sd))
  p <- p + geom_point(aes(color=as.factor(init.vscaling)), alpha=0.5)                         
  p <- p + facet_grid(Rbr ~ alphabr)
  p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
  p <- p + ggtitle(title)
  if (INTERACTIVE) {
    p
  }
  ggsave(filename=paste0(IMGPATH, paste0("newpdist/velocity_maxsize_pdistsd_t2000_s", s, ".jpg")), plot=p)
}

# color = sigma, by velocity => ALPHA does not have an impact on the final distribution.
for (s in unique(cluall$init.vscaling)) {
  title = paste0("Scatterplot max cluster size ~ st.d. pair-wise distance by alpha and R\ninitial velocity = ",s)
  p <- ggplot(cluall[cluall$t == 2000 & cluall$init.vscaling == s,], aes(x = size.max, y = pdist.sd))
  p <- p + geom_point(aes(color=as.factor(sigma)), alpha=0.5)
  p <- p + facet_grid(Rbr ~ alphabr)
  p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
  p <- p + ggtitle(title)
  if (INTERACTIVE) {
    p
  }
  ggsave(filename=paste0(IMGPATH, paste0("newpdist/maxsize_pdistsd_by_v", s, ".jpg")), plot=p)
}

# Summary of above, without alpha
title = paste0("Scatterplot max cluster size ~ st.d. pair-wise distance \n by initial velocity and R")
p <- ggplot(cluall[cluall$t == 2000, ], aes(x = size.max, y = pdist.sd))
p <- p + geom_point(aes(color=as.factor(sigma)), alpha=0.5)
p <- p + facet_grid(Rbr ~ init.vscaling)
p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
p <- p + ggtitle(title)
if (INTERACTIVE) {
  p
}
ggsave(filename=paste0(IMGPATH, paste0("newpdist/maxsize_pdistsd_by_v_and_r.jpg")), plot=p)


# OLD


 
title = "Clustering and Progress"
p <- ggplot(macro, aes(x=fromtruth.avg, y=count.avg.rate))
p <- p + geom_point(aes(color=t))
p

title = "Clustering and Progress"
p <- ggplot(macro, aes(x=count.avg.rate, y=fromtruth.avg))
p <- p + geom_point(aes(color=t))
p <- p + facet_grid( . ~ tbr)
p

                         
title = "Clustering and Progress"
p <- ggplot(macro, aes(x = fromtruth.avg, y = count.avg.rate))
p <- p + geom_point(aes(color=fromtruth.avg))
p <- p + scale_colour_gradient2()
p <- p + facet_wrap(~ tbr)                         
p

title = "Clustering and Progress"
p <- ggplot(clu, aes(x = size.max, y = fromtruth.avg))
p <- p + geom_point(aes(color=R))                         
#p <- p + scale_colour_gradient2()
p <- p + facet_grid(init.vscaling ~ tbr)                         
p

                         

                         
title = "Clustering and Progress"
p <- ggplot(macro, aes(x=fromtruth.avg, y=count.avg.rate))
p <- p + geom_point(aes(color=t))
p


                         
title = "Average size biggest cluster and distance from truth"
p <- ggplot(macro, aes(x = maxsize.avg, y = fromtruth.avg))
p <- p + geom_point(aes(color=fromtruth.avg))
p <- p + scale_color_gradient(low = "lightblue", high = "red", space = "Lab", na.value = "grey50", 
                              guide = "colourbar")
p <- p + ggtitle(title)
p

title = "Average cluster size and distance from truth"
p <- ggplot(macro, aes(x = meansize.avg, y = fromtruth.avg))
p <- p + geom_point(aes(color=fromtruth.avg))
p <- p + scale_color_gradient(low = "lightblue", high = "red", space = "Lab", na.value = "grey50", 
                              guide = "colourbar")
p <- p + ggtitle(title)
p
                         
title = "Average cluster size and distance from truth"
p <- ggplot(macro, aes(x = count.avg, y = fromtruth.avg))
p <- p + geom_point(aes(color=fromtruth.avg))
p <- p + scale_color_gradient(low = "lightblue", high = "red", space = "Lab", na.value = "grey50", 
                              guide = "colourbar")
p <- p + ggtitle(title)
p  

# All together
                         
title = "Clustering and distance from truth"
p <- ggplot(macro, aes(y = fromtruth.avg))
p <- p + geom_point(aes(x = count.avg, color="count"))
p <- p + geom_point(aes(x = maxsize.avg, color="max size"))                        
p <- p + geom_point(aes(x = meansize.avg, color="mean size"))
p
                         
p <- p + scale_color_gradient(low = "lightblue", high = "red", space = "Lab", na.value = "grey50", 
                              guide = "colourbar")
p <- p + ggtitle(title)
p  
                  
