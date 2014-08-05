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
loadData <- function(DUMPDIR, DIR, LOAD_T = 0, T = 2000) {

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


  if (LOAD_T) {
    macro <- macro[macro$t == T,]
  }
  
  clu <- merge(params, macro, by=c("simname","simcount", "run"))
 

  clu$simname <- as.character(clu$simname)
  clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
  clu$simname <- as.factor(clu$simname)
  clu$simcount <- as.factor(clu$simcount)
  #clu$t <- as.factor(clu$t)

  cl <- clu[clu$t == T,]

  return(clu)
}
#
loadDataAgents <- function(DUMPDIR, DIR, LOAD_T = 0, T = 2000) {

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


  #########
  # Macro #
  #########
  macro <- read.table('clusters_macro.csv', head=TRUE, sep=",")

  if (LOAD_T) {
    macro <- macro[macro$t == T,]
  }

  # Merging
  clu <- merge(params, macro, by=c("simname","simcount", "run"))

  ##########
  # Agents #
  ##########
  agents <- read.table('agents.csv', head=TRUE, sep=",")

  if (LOAD_T) {
    agents <- agents[agents$t == T,]
  }

  agents <- subset(agents, select=-c(coverage, coverage.cum))
  colnames(agents) <- c("simname", "simcount", "run", "t", "a.speed.avg",
                        "a.speed.sd", "a.move.avg", "a.move.sd",
                        "a.fromtruth.avg", "a.fromtruth.sd", "a.pdist.mean",
                        "a.pdist.sd")

  clu <- merge(clu, agents, by=c("simname","simcount", "run", "t"))
  
  clu$simname <- as.character(clu$simname)
  clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
  clu$simname <- as.factor(clu$simname)
  clu$simcount <- as.factor(clu$simcount)
  #clu$t <- as.factor(clu$t)


  
  cl <- clu[clu$t == T,]

  return(clu)
}
#
theme_white <- function() {
  theme_update(panel.background = element_blank())
}

#
myLabeller <- function(var, value){
  value <- as.character(value)
  if (var == "R") {
    value[value == 0.03] <- "Small Radius R = 0.03"
    value[value == 0.3] <- "Large Radius R = 0.3"
  } else if (var == "clbr") {
    value[value == "(0,1]"] <- "1 cluster"
    value[value == "(1,5]"] <- "2-5 clusters"
    value[value == "(5,10]"] <- "6-10\nclusters"
    value[value == "(10,20]"] <- "11-20\nclusters"
    value[value == "(20,30]"] <- "21-30\nclusters"
  }
  return(value)
}


myLabeller2 <- function(var, value){
  value <- as.character(value)
  if (var == "clbr") {
    value[value == "(0,1]"] <- "1 cluster"
    value[value == "(1,5]"] <- "2-5 clusters"
    value[value == "(5,10]"] <- "6-10 clusters"
    value[value == "(10,20]"] <- "11-20 clusters"
    value[value == "(20,30]"] <- "21-30 clusters"
  }
  return(value)
}

myLabeller3 <- function(var, value){
  value <- as.character(value)
  if (var == "convZone") {
    value[value == 0] <- "R <= 1"
    value[value == 1] <- "R > 1"
  }
  return(value)
}

#
myThemeMod <- theme(legend.position = "none",
                    axis.title.x = element_text(vjust=-1, size=24),
                    axis.title.y = element_text(vjust=-0.1, size=24),
                    plot.margin=unit(c(10,10,10,10),"mm"),                    
                    plot.title = element_text(vjust=3, size=24,face="bold")
                    )
#
limits <- aes(ymax = count + se, ymin = count - se)
limitsFt <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)
limitsCI <- aes(ymax = count + ci, ymin = count - ci)
#
theme_set(theme_bw(base_size = 30))
theme_white()
#
XINTERCEPT <- 0.15
