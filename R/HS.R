# HS analysis

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/few_big_groups-DIM-vs-ALPHA/"
PATH = "/opt/MATLAB_WORKSPACE/hs/dump/sigma_tau-2013-3-6-10-36/"

PATH = "/opt/MATLAB_WORKSPACE/hs/dump/R-alpha-noA-noB-2013-3-6-20-8/"
#PATH = "/opt/MATLAB_WORKSPACE/hs/dump/R-alpha-noA-noB-2013-3-6-20-8"


setwd(PATH)

params <- read.table('params.csv', head=TRUE, sep=",")
clusters <- read.table('clusters.csv', head=TRUE, sep=",")
params$simname <- as.factor(params$simname)
params$simcount <- as.factor(params$simcount)
params$run <- as.factor(params$run)
clusters$simname <- as.factor(clusters$simname)
clusters$simcount <- as.factor(clusters$simcount)
clusters$run <- as.factor(clusters$run)
# Transforms params in factors
clu <- merge(params, clusters, by=c("simname","simcount","run"))
for (n in names(clu[1:23])) {
  clu[, n] <- as.factor(clu[, n])      
}
clu$t <- as.factor(clu$t)




# START



IMGPATH <- paste0(PATH, "img/");


allPlots("R","alpha")



# POINTS
#p.conv <- ggplot(clu[clu$tau==1,], aes(t, fromtruth.avg))
#p.conv + geom_point(aes(colour = sigma)) + geom_jitter(aes(colour = sigma)) + plotScaleDis


