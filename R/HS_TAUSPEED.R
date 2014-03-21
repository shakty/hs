# TAU SPEED

source("/opt/MATLAB_WORKSPACE/hs/R/init.R")
library(lattice)


myLabeller <- function(var, value){
  value <- as.character(value)
  if (var == "R") {
    value[value== 0.03] <- "Small radius (R = 0.03)"
    value[value== 0.3] <- "Large Radius (R = 0.3)"
  }
#  else if
# Does not work at the moment.
#   else if (var == "alpha") {
#     value[value == 0.01] <- expression(alpha)
#     value[value == 0.5] <- "b" # expression(paste(alpha, ' 0.5'))
#     value[value == 0.99] <- "c" # expression(paste(alpha, ' 0.99'))
#   }
  return(value)
}


theme_set(theme_bw(base_size = 30))


myThemeMod <- theme(legend.position = "none",
                    axis.title.x = element_text(vjust=-1, size=24),
                    axis.title.y = element_text(vjust=-0.1, size=24),
                    plot.margin=unit(c(10,10,10,10),"mm"),                    
                    plot.title = element_text(vjust=3, size=24,face="bold")
                    )

## Data

SMALLRADIUS <- 1

DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'

if (SMALLRADIUS) {
  DIR <- 'final_tau_vs_speed/'
} else {
  DIR <- 'final_tau_vs_speed_largeR/'
}

PATH <- paste0(DUMPDIR, DIR)
IMGPATH <- paste0(DUMPDIR, "imgs/");
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/imgs/"))
}
setwd(PATH)


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


clu$tau10 <- cut(clu$tau, breaks=seq(0,1,0.02))
clu$v10 <- cut(clu$init.vscaling, breaks=seq(0,1,0.02))
clu$invtau <- 1 - clu$tau

mydata <- clu[clu$t == 2000,]


if (SMALLRADIUS) {
  title <- "Truth's signal strength vs Velocity multiplier - Small Radius (R = 0.03)"
  file <- "wireframe_tau_v_r003.jpeg"
} else {
  title <- "Truth's signal strength vs Velocity multiplier - Large Radius (R = 0.3)"
  file <- "wireframe_tau_v_r03.jpeg"
}

theseCol=heat.colors(150)

xlabstr <- list(expression(tau))
xlabstr <- "Signal's\n strength"

jpeg(filename=paste0(IMGPATH, file),
     res=300, quality=100, width=2000, height=2000)
print(      
      #oldpar <- par(mgp=c(10,10,10))
      wireframe(fromtruth.avg ~ invtau * init.vscaling, data = mydata,
                shade = TRUE,
                scales = list(arrows = FALSE),
                pretty = TRUE,
                main=title,
                zlim=c(0,0.25),
                zlab=list(cex=1.3, label="DIST"), xlab=list(cex=1.3, label=xlabstr),
                ylab=list(cex=1.3, label="Velocity\n multiplier", distance=3, at=10)
                )
      )
dev.off()



# npanel <- c(4, 2)
# rotx <- c(-50, -80)
# rotz <- seq(30, 300, length = npanel[1]+1)
# update(p[rep(1, prod(npanel))], layout = npanel,
#     panel = function(..., screen) {
#         panel.wireframe(..., screen = list(z = rotz[current.column()],
#                                            x = rotx[current.row()]))
#     })
# 
# 
# library(rgl)
# plot3d(mydata$invtau, mydata$init.vscaling, mydata$fromtruth.avg, type = "p")
