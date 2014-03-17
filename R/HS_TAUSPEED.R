# SPEED TEST: many vs 1 cluster

source("/opt/MATLAB_WORKSPACE/hs/R/init.R")



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

DUMPDIR <- '/home/stefano/HS/'
DIR <- 'final_tau_vs_speed/'

PATH <- paste0(DUMPDIR, DIR)
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
