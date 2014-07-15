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


# IT WAS:
# If consensus is not reached it has value -1. Replace with NA
#data[] <- lapply(data, function(x){replace(x, x == -1, NA)})

# IT IS NOW:
# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
#data$init.placement <- as.factor(data$init.placement)
# If consensus is not reached it has value -1. Replace with NA
#data[ , 15:28 ][ data[ , 15:28 ] == -1 ] <- NA
# 20.000
data[ , 15:28 ][ data[ , 11:22 ] == -1 ] <- 20000
# If no consensus was reached replace ccounts with NA
data[ , 29:40 ][ data[ , 23:34 ] == -1 ] <- NA
# Replaces everything, not good
#data[] <- lapply(data, function(x){replace(x, x == -1, 20000)})
data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)

  
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
DIR <- 'clusters_vs_progress_hidef/'

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
data[] <- lapply(data, function(x){replace(x, x == -1, 20000)})
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


ggsave(filename = paste0(IMGPATH, "clusters_vs_progress.jpg"),
       width=10, height=12, dpi=300)



mydata <- data[data$R == 0.03,]
mydata$moreThan9 <- mydata$init.ccount > 9

title <- "Only R = 0.03. Smoothed functions."
p <- ggplot(mydata, aes(init.ccount, consensus75))
p <- p + geom_smooth(alpha=0.5, aes(color = as.factor(init.placement)))
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p

ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_r003.jpg"),
       width=10, height=12, dpi=300)

##########
# Added later in session with Michael
#########

mydata <- data[data$R == 0.03 & data$alpha == 0.5 &
               (data$init.placement == 0.1 |
                data$init.placement == 0.15 |
                data$init.placement == 0.2 |
                data$init.placement == 0.25 |
                data$init.placement == 0.3 |
                data$init.placement == 0.35 |
                data$init.placement == 0.4 |
                data$init.placement == 0.45 |
                data$init.placement == 0.5
                ),]
mydata$moreThan9 <- mydata$init.ccount > 9

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata, aes(init.ccount, consensus75))
p <- p + geom_smooth(se = FALSE, aes(color = as.factor(init.placement)), size=2)
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
#p <- p + facet_grid(alpha~.)
#p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p <- p + myThemeMod + theme(legend.position = c(0.8, 0.3),
                            legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(vjust=3, size=16,face="bold"),
                            legend.text = element_text(size=14),
                            legend.key.width = unit(1.5, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA))
                            )
p

mydata.summary <- summarySE(mydata, "consensus75", c("init.ccount", "init.placement"), na.rm = TRUE)

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary, aes(init.ccount, consensus75))
p <- p + geom_point(aes(color = as.factor(init.placement)), size=4)
p <- p + geom_line(aes(color = as.factor(init.placement)))
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
#p <- p + facet_grid(alpha~.)
#p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p <- p + myThemeMod + theme(legend.position = c(0.8, 0.3),
                            legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(vjust=3, size=16,face="bold"),
                            legend.text = element_text(size=14),
                            legend.key.width = unit(1.5, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA)
                            )
p

ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_r003_new.jpg"),
       width=10, height=10, dpi=300)




# Interpolate

s2 <- aggregate(consensus75 ~ init.placement + init.ccount, data = mydata)
wireframe(consensus75 ~ init.ccount + init.placement, data = s2,
shade = TRUE,
scales = list(arrows = FALSE),
pretty = TRUE)

s3 <- aggregate(ccount50 ~ init.placement + init.ccount, data = mydata, FUN = mean)
wireframe(ccount50 ~ init.ccount + init.placement, data = s3)

shade = TRUE,
scales = list(arrows = FALSE),
pretty = TRUE)


###########


mydata <- data[data$R == 0.03,]
mydata$moreThan9 <- mydata$init.ccount > 9

title <- "Only R = 0.03. More than 9 clusters?"
p <- ggplot(mydata, aes(init.ccount, consensus75))
p <- p + geom_smooth(alpha=0.5, method="lm", aes(color = as.factor(init.placement)))
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p <- p + facet_grid( ~ moreThan9, labeller = myLabeller)
p

ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_r003_more_than9.jpg"),
       width=10, height=12, dpi=300)

###############

data$clbr <- cut(data$init.ccount, breaks=seq(0,30,2))
data$plbr <- cut(data$init.placement, breaks=seq(0,0.5,0.05))

mydata <- summarySE(data[data$R == 0.03 & data$alpha == 0.01,], "consensus75",
                    c("simcount", "run", "init.placement", "init.ccount"))

mydataR3 <- summarySE(data[data$R == 0.3 & data$alpha == 0.01,], "consensus75",
                    c("simcount", "run", "init.placement", "init.ccount"))


mydata <- data[data$R == 0.03 & data$alpha == 0.5 & data$init.ccount == 12,]

mydata2 <- summarySE(mydata, "consensus75", c("init.placement"), na.rm=TRUE)


p <- ggplot(mydata, aes(init.placement, consensus75))
p <- p + stat_smooth(aes(outfit=fit<<-..y..))
p

p <- ggplot(mydata2, aes(init.placement, consensus75))
p <- p + geom_smooth()
#p <- p + geom_line()
p <- p + geom_jitter()
#p <- p + geom_errorbar(aes(ymin = consensus75 - se, ymax = consensus75 + se))
p



wireframe(consensus75 ~ init.placement * init.ccount, data[data$R == 0.03 & data$alpha == 0.5,],
          shade = TRUE,
          scales = list(arrows = FALSE),
          pretty = TRUE)
          #main=title,
          #zlim=c(0,0.25),
          #zlab=list(cex=1.3, label="DIST"), xlab=list(cex=1.3, label=xlabstr),
          #ylab=list(cex=1.3, label="Velocity\n multiplier", distance=3, at=10))

wireframe(consensus75 ~ plbr * clbr, data[data$R == 0.03 & data$alpha == 0.01,],
          shade = TRUE,
          scales = list(arrows = FALSE),
          pretty = TRUE)

#####################

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
p <- ggplot(mydata, aes(init.band.i, ccount25, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p


ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_exp2.jpg"),
       width=10, height=10, dpi=300)

mydata <- data
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(as.factor(init.band.i), ccount50, color = R))
p <- p + geom_boxplot(notch=TRUE)
#p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p

ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_exp2_boxplot.jpg"),
       width=10, height=10, dpi=300)


mydata <- data
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount50, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(alpha  ~ R, labeller = myLabeller)
p

ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_exp2_by_alpha.jpg"),
       width=10, height=10, dpi=300)



colnamesData <- colnames(data)


data.subset <- subset(data, select=colnamesData[seq(1,27)])
data.melted <- melt(data.subset, id=colnamesData[seq(1,14)])

data.melted.names <- colnames(data.melted)
data.melted.names[15] <- "consensus.share"
data.melted.names[16] <- "time"
colnames(data.melted) <- data.melted.names

data.melted$consensus.share <- substring(data.melted$consensus.share, 10)
data.melted$consensus.share <- as.factor(data.melted$consensus.share)


data.subset <- subset(data, select=colnamesData[c(seq(1,14),seq(28,40))])
data.melted2 <- melt(data.subset, id=colnamesData[seq(1,14)])

data.melted2.names <- colnames(data.melted2)
data.melted2.names[15] <- "consensus.share"
data.melted2.names[16] <- "ccount"
colnames(data.melted2) <- data.melted2.names

data.melted2 <- subset(data.melted2, select=data.melted2.names[c(2,3,15,16)])
data.melted2$consensus.share <- substring(data.melted2$consensus.share, 7)
data.melted2$consensus.share <- as.factor(data.melted2$consensus.share)


data.melted <- merge(data.melted, data.melted2, id=c("simcount","run","consensus.share"))

mydata <- data.melted[data.melted$consensus.share != 1 &
                      data.melted$consensus.share != 10 &
                      data.melted$consensus.share != 100 &
                      data.melted$consensus.share != 25 &
                      data.melted$consensus.share != 75 &
                      data.melted$R == 0.03
                      ,]

title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount, color = consensus.share))
#p <- p + geom_jitter(size=2)
p <- p + geom_smooth(method="lm", se = FALSE, size = 2)
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod #+ ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1)) + ylim(c(0,12))
p <- p + scale_color_hue(name="Share of\nconsensus")
p <- p + theme(legend.position = c(0.5, 0.9), legend.direction = "horizontal",
legend.background = element_rect(fill = "white", colour = "grey"),
legend.title = element_text(vjust=3, size=16,face="bold"),
legend.text = element_text(size=14),
legend.key.width = unit(1.5, "cm"),
legend.key =  element_rect(fill = "white", colour = NA)
)
p <- p + guides(col=guide_legend(ncol=4))
p

ggsave(filename = paste0(IMGPATH, "clusters_vs_progress_exp2_rainbow.jpg"), dpi=300)


title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount, color = consensus.share))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(method="lm", se = FALSE, size = 2)
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
#p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p


# Interpolate

s2 <- aggregate(consensus75 ~ init.band.i + init.ccount, data = mydata, FUN = mean)
wireframe(consensus75 ~ init.ccount + init.placement, data = s2,
shade = TRUE,
scales = list(arrows = FALSE),
pretty = TRUE)
