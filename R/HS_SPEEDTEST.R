# SPEED TEST: many vs 1 cluster

source("/opt/MATLAB_WORKSPACE/hs/R/init.R")
library(texreg)
library(mediation)
library(grid)
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

DUMPDIR <- '/mnt/tmp/dump/SPEEDTEST/'

DIR  <- 'attrLinear_navnp_RFull_SpeedTest_epsilon/'

DUMPDIR <- '/home/stefano/HS/'

DIR <- 'attrLinear_navnp_SpeedTest_fullscan/'
DIR <- 'speedtest_final_R_alpha/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

data <- read.table('speedtest.csv', head = T, sep = ",")

data$smallR <- as.numeric(data$R <= 0.1)
data$bigR <- as.numeric(data$R > 0.1)

# If consensus is not reached it has value -1. Replace with NA
data[] <- lapply(data, function(x){replace(x, x == -1, NA)}) 
  
mean(data[data$smallR == 1,]$consensus100, na.rm=TRUE)
mean(data[data$smallR == 0,]$consensus100, na.rm=TRUE)

### Mediation and Moderation Analysis

mydata <- data[complete.cases(data[,c("consensus75","ccount50")]) & data$alpha ==0.5,]

# Step1 Regress INDEPENDENT on DEPENDENT: init.ccount on consensus75

fit1a <- lm(consensus75 ~ smallR, data = mydata)
summary(fit1a)

fit1aa <- lm(consensus75 ~ ccount50, data = mydata)
summary(fit1aa)


# Step2 Regress INDEPENDENT on MEDIATOR: R and alpha on consensus75

fit1b <- lm(ccount50 ~ smallR, data = mydata)
summary(fit1b)

# Step3 Regress MEDIATOR on DEPENDENT: R and alpha on consensus75

fit1c <- lm(consensus75 ~ smallR + ccount50, data = mydata)
summary(fit1c)


# Mediation Analysis: Preacher & Hayes (2004) Bootstrap Method


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="smallR", mediator="ccount50")

summary(contcont)

plot(contcont)




#######################
### Cluster vs Clusters
#######################

DIR <- 'clusterVsClusters/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

IMGPATH <- paste0(DUMPDIR, "imgs/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

data <- read.table('speedtest.csv', head = T, sep = ",")
# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
data$init.placement <- as.factor(data$init.placement)
# If consensus is not reached it has value -1. Replace with NA
data[] <- lapply(data, function(x){replace(x, x == -1, NA)})
data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)


# Log data
# data$consensus75 <- log(data$consensus75)

# Good fit - 3
mydata <- data[data$alpha == 0.01 & data$init.placement == 0.4 & data$R == 0.03,]
fit1 <- lm(consensus75 ~ init.ccount, data = mydata)
summary(fit1)

# Medium fit - 1
mydata <- data[data$init.placement == 0.4 & data$R == 0.03,]
fit2 <- lm(consensus75 ~ init.ccount, data = mydata)
summary(fit2)

# Not so good fit - 4
mydata <- data[data$alpha == 0.01 & data$init.placement == 0.4 & data$R == 0.3,]
fit3 <- lm(consensus75 ~ init.ccount, data = mydata)
summary(fit3)

# Even worse - 2
mydata <- data[data$init.placement == 0.4 & data$R == 0.3,]
fit4 <- lm(consensus75 ~ init.ccount, data = mydata)
summary(fit4)

#### Mediation and Moderation Analysis

mydata <- data[data$init.placement == 0.4,]
mydata <- mydata[complete.cases(mydata[,c("consensus75")]),]

mydata$R003 <- mydata$R == 0.03
mydata$R003A001 <- mydata$R == 0.03 & mydata$alpha == 0.01
mydata$R03 <- mydata$R == 0.3
mydata$R03A001 <- mydata$R == 0.3 & mydata$alpha == 0.01

# Step1 Regress INDEPENDENT on DEPENDENT: init.ccount on consensus75

# Good fit - 3
fit1a <- lm(consensus75 ~ R003, data = mydata)
summary(fit1a)


# Step2 Regress INDEPENDENT on MEDIATOR: R and alpha on consensus75

# Good fit - 3
fit1b <- lm(init.ccount ~ R003, data = mydata)
summary(fit1b)


# Step3 Regress MEDIATOR on DEPENDENT: R and alpha on consensus75

# Good fit - 3
fit1c <- lm(consensus75 ~ init.ccount + R003, data = mydata)
summary(fit1c)

# Mediation Analysis: Preacher & Hayes (2004) Bootstrap Method

# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="R003", treat.value=TRUE, mediator="init.ccount")

summary(contcont)

plot(contcont)


#### Tables and Plots

texreg(list(fit2, fit4, fit1, fit3))

texreg(list(fit2, fit4, fit1, fit3), booktabs = TRUE, dcolumn = TRUE)

p <- ggplot(mydata, aes(init.ccount, consensus75))
p <- p + geom_jitter(aes(color=init.ccount))
p

mydata <- data[data$init.placement == 0.4,]
title <- "Number of clusters and time to reach a consensus"
p <- ggplot(mydata, aes(init.ccount, consensus75, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)

p + facet_grid(. ~ R, labeller = myLabeller)

ggsave(filename = paste0(IMGPATH, "race_init-cc.jpg"),
       width=10, height=10, dpi=600)

pAlpha <- p + facet_grid(alpha ~ R, labeller = myLabeller)

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(pAlpha)
grob[["grobs"]][[13]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.01"))
grob[["grobs"]][[14]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.5"))
grob[["grobs"]][[15]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.99"))

jpeg(filename = paste0(IMGPATH, "race_init-cc_alpha.jpg"),
     quality=100, width=6000, height=6000, res=600)
grid.newpage()
grid.draw(grob)
dev.off()


# Distribution of time to consensus

# Add this if wanna see by cluster as well.
data$clbr <- cut(data$init.ccount, breaks=c(0,1,5, 10, 20, 30))

if (exists("summaryData")) {
  rm(summaryData)
}
if (exists("summaryNew")) {
  rm(summaryNew)
}
counter <- 1
for (c in seq(10,100,10)) {
  myVar <- c(paste0('consensus', c))
  # Remove clbr if not needed
  summaryNew <- summarySE(data[data$alpha != 0.99,], c(myVar), c('R', 'clbr'), na.rm=TRUE)
  names(summaryNew) <- sub(myVar, "time", names(summaryNew))
  summaryNew$share <- c
  if (exists("summaryData")) {
    newRowNums <- c(counter:(nrow(summaryNew)+counter-1)); newRowNums
    rownames(summaryNew) <- newRowNums
    summaryData <- rbind(summaryData, summaryNew)
  } else {
    summaryData <- summaryNew
  }
  counter <- counter + nrow(summaryNew)
}


title <- 'Temporal evolution of consensus building \n by number of initial clusters'
p <- ggplot(summaryData, aes(x = share, weight=time, group=R, fill=R))
p <- p + geom_bar(aes(y=time), stat = "identity", position = "dodge")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se), position = "dodge")
p <- p + xlab('Share of Consensus') + ylab('Time Passed')
# Remove clbr if not needed
p <- p + facet_grid(.~clbr)
p + ggtitle(title)  + myThemeMod + theme(legend.position = c(1.1, 0.5),
                                         plot.margin = unit(c(10,30,10,10),"mm"),
                                         axis.text.x = element_text(size=15))

ggsave(filename = paste0(IMGPATH, "race_distribution_init-cc.jpg"),
       width=10, height=5, dpi=600)


#########################
### Clusters vs Progress
#########################

DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'

DIR <- 'clusters_vs_progress/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

IMGPATH <- paste0(DUMPDIR, "imgs/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

data <- read.table('speedtest.csv', head = T, sep = ",")
# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
#data$init.placement <- as.factor(data$init.placement)
# If consensus is not reached it has value -1. Replace with NA
data[] <- lapply(data, function(x){replace(x, x == -1, NA)})
data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)


# Add this if wanna see by cluster as well.
data$clbr <- cut(data$init.ccount, breaks=c(0,1,5, 10, 20, 30))

mydata <- data
title <- "Number of clusters vs progress and time to reach a consensus"
p <- ggplot(mydata, aes(init.ccount, consensus75, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(init.placement  ~ R, labeller = myLabeller)
p

mydata <- data
title <- "Number of clusters vs progress and time to reach a consensus"
p <- ggplot(mydata, aes(init.placement, consensus75, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Clusters placed at distance") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p




levelplot(consensus75 ~ init.placement * init.ccount, data= data[data$init.placement > 0.2 & data$R == 0.03 & data$alpha == 0.5,],
          shade = TRUE,
          scales = list(arrows = FALSE))
          pretty = TRUE)
          #main=title,
          #zlim=c(0,0.25),
          #zlab=list(cex=1.3, label="DIST"), xlab=list(cex=1.3, label=xlabstr),
          #ylab=list(cex=1.3, label="Velocity\n multiplier", distance=3, at=10))

library(rgl)
mydata <- data[data$R == 0.03 & data$alpha == 0.5,]
plot3d(mydata$init.ccount, mydata$init.placement, mydata$consensus75)

mydata <- data[data$init.ccount == 30,]
title <- "Number of clusters vs initial progress"
p <- ggplot(mydata, aes(init.placement, ccount50, color = R))
p <- p + geom_jitter(size=2)
#p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p

mydata <- data[data$init.ccount == 30,]
title <- "Number of clusters vs initial progress"
p <- ggplot(mydata, aes(init.placement, ccount50, color = R))
p <- p + geom_jitter(size=2)
#p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p


mydata <- data
title <- "Initial progress vs Number of clusters"
p <- ggplot(mydata, aes(init.placement, ccount50, color = R))
p <- p + geom_jitter(size=2)
#p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Clusters placed at a distance from") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(clbr  ~ R, labeller = myLabeller)
p

mydata <- data[data$init.placement == 0.4,]
title <- "Number of clusters and time to reach a consensus"
p <- ggplot(mydata, aes(init.ccount, consensus75, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Clusters placed at a distance from") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(. ~ R, labeller = myLabeller)
p

ggsave(filename = paste0(IMGPATH, "race_init-cc.jpg"),
       width=10, height=10, dpi=600)

pAlpha <- p + facet_grid(alpha ~ R, labeller = myLabeller)

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(pAlpha)
grob[["grobs"]][[13]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.01"))
grob[["grobs"]][[14]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.5"))
grob[["grobs"]][[15]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.99"))

jpeg(filename = paste0(IMGPATH, "race_init-cc_alpha.jpg"),
     quality=100, width=6000, height=6000, res=600)
grid.newpage()
grid.draw(grob)
dev.off()


# Distribution of time to consensus

# Add this if wanna see by cluster as well.
data$clbr <- cut(data$init.ccount, breaks=c(0,1,5, 10, 20, 30))

if (exists("summaryData")) {
  rm(summaryData)
}
if (exists("summaryNew")) {
  rm(summaryNew)
}
counter <- 1
for (c in seq(10,100,10)) {
  myVar <- c(paste0('consensus', c))
  # Remove clbr if not needed
  summaryNew <- summarySE(data[data$alpha != 0.99,], c(myVar), c('R', 'clbr'), na.rm=TRUE)
  names(summaryNew) <- sub(myVar, "time", names(summaryNew))
  summaryNew$share <- c
  if (exists("summaryData")) {
    newRowNums <- c(counter:(nrow(summaryNew)+counter-1)); newRowNums
    rownames(summaryNew) <- newRowNums
    summaryData <- rbind(summaryData, summaryNew)
  } else {
    summaryData <- summaryNew
  }
  counter <- counter + nrow(summaryNew)
}


title <- 'Temporal evolution of consensus building \n by number of initial clusters'
p <- ggplot(summaryData, aes(x = share, weight=time, group=R, fill=R))
p <- p + geom_bar(aes(y=time), stat = "identity", position = "dodge")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se), position = "dodge")
p <- p + xlab('Share of Consensus') + ylab('Time Passed')
# Remove clbr if not needed
p <- p + facet_grid(.~clbr)
p + ggtitle(title)  + myThemeMod + theme(legend.position = c(1.1, 0.5),
                                         plot.margin = unit(c(10,30,10,10),"mm"),
                                         axis.text.x = element_text(size=15))

ggsave(filename = paste0(IMGPATH, "race_distribution_init-cc.jpg"),
       width=10, height=5, dpi=600)


#########################
### Clusters vs Progress: RBANDS
#########################

DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'
DIR <- 'clusters_vs_progress_rbands01/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

IMGPATH <- paste0(DUMPDIR, "imgs/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

data <- read.table('speedtest.csv', head = T, sep = ",")
# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
data$init.placement <- as.factor(data$init.placement)
# If consensus is not reached it has value -1. Replace with NA
data[] <- lapply(data, function(x){replace(x, x == -1, NA)})
data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)


mydata <- data
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount75, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p

mydata <- data
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(as.factor(init.band.i), ccount75, color = R))
p <- p + geom_boxplot(notch=TRUE)
#p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p



mydata <- data
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount10, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(alpha  ~ R, labeller = myLabeller)
p




