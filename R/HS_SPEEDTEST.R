# SPEED TEST: many vs 1 cluster

source("/opt/MATLAB_WORKSPACE/hs/R/init.R")


DUMPDIR <- '/mnt/tmp/dump/SPEEDTEST/'

DIR  <- 'attrLinear_navnp_RFull_SpeedTest_epsilon/'

DUMPDIR <- '/home/stefano/HS/'

DIR <- 'attrLinear_navnp_SpeedTest_fullscan/'
DIR <- 'speedtest_final_R_alpha/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

data <- read.table('speedtest.csv', head = T, sep = ",")

data$cluster <- as.numeric(data$R < 0.06)

# If consensus is not reached it has value -1. Replace with NA
data[] <- lapply(data, function(x){replace(x, x == -1, NA)}) 


mean(data[data$cluster == 1,]$t)
mean(data[data$cluster == 0,]$t)

summaryData <- summarySE(data, c('consensus75'), c('cluster','sigma'), na.rm=TRUE)
summaryData


fit <- lm(consensus75 ~ ccount50, data=data)
summary(fit)
