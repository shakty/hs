# SPEED TEST: many vs 1 cluster

source("/opt/MATLAB_WORKSPACE/hs/R/init.R")


DUMPDIR <- '/mnt/tmp/dump/SPEEDTEST/'

DIR  <- 'attrLinear_navnp_RFull_SpeedTest_epsilon/'

DUMPDIR <- '/home/stefano/HS/'
DIR <- 'attrLinear_navnp_SpeedTest_fullscan/'

DIR <- 'final_R_alpha/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

data <- read.table('speedtest.csv', head = T, sep = ",")

data$cluster <- as.numeric(data$R < 0.06)

mean(data[data$cluster == 1,]$t)
mean(data[data$cluster == 0,]$t)

summaryData <- summarySE(data, c('t'), c('cluster','sigma'))
summaryData



