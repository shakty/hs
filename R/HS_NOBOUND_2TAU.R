# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR 
DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'

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
loadData <- function(DUMPDIR, DIR, TWO_THOUSANDS = 0) {

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


  if (TWO_THOUSANDS) {
    macro <- macro[macro$t == 2000,]
  }
  
  clu <- merge(params, macro, by=c("simname","simcount", "run"))
 

  clu$simname <- as.character(clu$simname)
  clu$simname <- substr(clu$simname, nchar(clu$simname)-1, nchar(clu$simname))
  clu$simname <- as.factor(clu$simname)
  clu$simcount <- as.factor(clu$simcount)
  #clu$t <- as.factor(clu$t)

  cl <- clu[clu$t == 2000,]

  return(clu)
}

theme_white <- function() {
  theme_update(panel.background = element_blank())
}



myLabeller <- function(var, value){
  value <- as.character(value)
  if (var == "R") {
    value[value== 0.03] <- "Small radius (R = 0.03)"
    value[value== 0.3] <- "Large Radius (R = 0.3)"
  } 
  return(value)
}

myThemeMod <- theme(axis.title.x = element_text(vjust=-1, size=24),
                    axis.title.y = element_text(vjust=-0.1, size=24),
                    plot.margin=unit(c(10,10,10,10),"mm"),                    
                    plot.title = element_text(vjust=3, size=24,face="bold"),
                    legend.background = element_rect(fill = "white", color="grey"),
                    legend.title = element_blank(),
                    legend.text = element_text(size=16),
                    legend.key.width = unit(1.5, "cm"),
                    legend.key = element_rect(fill = "white", colour = "white")
                    )

limits <- aes(ymax = count + se, ymin = count - se)
limitsFt <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)
limitsCI <- aes(ymax = count + ci, ymin = count - ci)

theme_set(theme_bw(base_size = 30))
theme_white()

XINTERCEPT <- 0.15

## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/NOBOUND/twotaus/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

## R ##
#######
cl <- loadData(DUMPDIR, 'nobound_R_tau/', 1)

clTau1 <- loadData(DUMPDIR, 'nobound_R_tau1/', 1)

cl <- rbind(cl,clTau1)
#cl <- rbind(cl[cl$tau == 50,], clTau1[clTau1$alpha == 0.5,])


# CL
summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("R", "tau"), na.rm=TRUE)
summaryCl.50 <- summarySE(cl[cl$t == 2000 & cl$tau == 50,], c("count"), c("R", "tau"), na.rm=TRUE)
summaryCl.1 <- summarySE(cl[cl$t == 2000 & cl$tau == 1,], c("count"), c("R", "tau"), na.rm=TRUE)

title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryCl[summaryCl$tau %in% c(1,10,50,100),] , aes(R, count, color = as.factor(tau)))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Number of Clusters')
p <- p  + myThemeMod
p

title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryCl[summaryCl$tau %in% c(1,10,50,100),] , aes(R, count))
p <- p + geom_bar(stat = "identity", position = "dodge", aes(fill=as.factor(tau)), width=0.01)
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Number of Clusters')
p <- p  + myThemeMod
p



c1 <- "#54aef3ff"
c2 <- "#18cad06b"
c3 <- "#d4392dff"

title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryCl.50, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=as.factor(tau), width=0.01))
p <- p + geom_errorbar(limits, position="dodge")
p <- p + geom_bar(data = summaryCl.1, stat = "identity", position="dodge", aes(fill=as.factor(tau), width=0.01))
p <- p + geom_errorbar(data = summaryCl.1, limits, position="dodge")
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Number of Clusters')
p <- p + scale_fill_discrete(name = "", labels = c(
                                          expression(paste(tau," = 1")),
                                          expression(paste(tau," = 50"))))
p <- p  + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "nobound_R_tau_cc.svg"),
       plot = p, width=10, height=5, dpi=300)


# FT
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("R", "tau"), na.rm=TRUE)
summaryFt.50 <- summarySE(cl[cl$t == 2000 & cl$tau == 50 & cl$R <= 0.3,], c("fromtruth.avg"), c("R", "tau"), na.rm=TRUE)
summaryFt.1 <- summarySE(cl[cl$t == 2000 & cl$tau == 1 & cl$R <= 0.3,], c("fromtruth.avg"), c("R", "tau"), na.rm=TRUE)

title <- 'Distance from truth by radius of influence'
p <- ggplot(summaryFt.50, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=as.factor(tau), width=0.01))
p <- p + geom_errorbar(limitsFt, position="dodge")
p <- p + geom_bar(data = summaryFt.1, stat = "identity", position="dodge", aes(fill=as.factor(tau), width=0.01))
p <- p + geom_errorbar(data = summaryFt.1, limitsFt, position="dodge")
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
#p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Distance from Truth')
#p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "nobound_R_tau_ft.svg"),
       plot = p, width=10, height=5, dpi=300)



title <- 'Cluster counts by radius of influence'
p <- ggplot(summaryFt[summaryFt$tau %in% c(1,10,50,100),] , aes(R, fromtruth.avg, color = as.factor(tau)))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(limitsFt)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
#p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Number of Clusters')
p <- p  + myThemeMod + scale_x_log10()
p

library("scales")

summaryFt$rescaleOffset <- summaryFt$fromtruth.avg + summaryFt$tau
scalerange <- range(summaryFt$fromtruth.avg)
gradientends <- scalerange + rep(c(0,100,200), each=2)
colorends <- c("green", "red", "white", "green", "white", "blue")


title <- 'Distance from truth by radius of influence'
p <- ggplot(summaryFt, aes(R, fromtruth.avg, group = as.factor(tau)))
p <- p + geom_bar(stat = "identity", aes(fill=rescaleOffset, width=0.01))
p <- p + geom_errorbar(limitsFt )
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
p <- p + xlab('Radius of Influence') + ylab('Avg. Distance from Truth')
#p <- p + scale_fill_gradientn(colours = colorends, values = rescale(gradientends))
p <- p + myThemeMod 
p



# Test

nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba$Name <- with(nba, reorder(Name, PTS))

nba.m <- melt(nba)
nba.s <- ddply(nba.m, .(variable), transform,
               rescale = scale(value))

## ALPHA ##
###########

cl <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1)

clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)

cl <- rbind(cl, clTau1)

cl$alphabrk <- cut(cl$alpha, seq(0,1,0.05))

# CL
summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("alphabrk", "R", "tau"), na.rm=TRUE)

# Reverse the order of a discrete-valued axis
# Get the levels of the factor
flevels <- levels(summaryCl$alphabrk)
# Reverse the order
flevels <- rev(flevels)


title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl[summaryCl$tau %in% c(1,10,50,100),],
            aes(alphabrk, count, color = as.factor(tau), group = as.factor(tau)))
p <- p + geom_point() + geom_line()
p <- p + facet_grid(~ R, labeller = myLabeller)
#p <- p + ylim(0,29)
#p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + scale_x_discrete(limits=flevels, labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + xlab('Strength of Social Influence') + ylab('Number of Clusters')
p <- p + myThemeMod + theme(strip.background = element_blank())
p

# CL
summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("alpha", "R"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
#p <- p + ylim(0,29)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab('Strength of Social Influence') + ylab('Number of Clusters')
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "nobound_alpha_tau1_cc.svg"),
       plot = p, width=10, height=5, dpi=300)

# FT

# does not work
myalphas <- c(0.01, seq(0.1,0.9,0.1),0.99)
# works
myalphas <- c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99)

summaryFt <- summarySE(cl[cl$t == 2000 &
                          cl$alpha %in% myalphas &
                          cl$tau %in% c(1,10,50,100)
                          ,], c("fromtruth.avg"), c("alpha", "R", "tau"), na.rm=TRUE)

#original
title <- 'Distance from truth vs Strength of social influence'
p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg,
                           color = as.factor(tau), group = as.factor(tau)))
p <- p + geom_point() + geom_line(alpha = 0.5) + geom_errorbar(limitsFt, width = 0.01)
p <- p + facet_grid(~ R, labeller = myLabeller)
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab("Strength of social influence") + ylab('Avg. Distance from Truth')
#p <- p + ylim(0,0.3)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
#p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod + theme(strip.background = element_blank())
p


title <- 'Distance from truth vs Strength of social influence'
p <- ggplot(summaryFt[summaryFt$tau %in% c(1,10,50,100),],
            aes(alphabrk, fromtruth.avg, color = as.factor(tau), group = as.factor(tau)))
p <- p + geom_point() + geom_line()
p <- p + facet_grid(~ R, labeller = myLabeller)
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
p <- p + xlab("Strength of social influence") + ylab('Avg. Distance from Truth')
#p <- p + ylim(0,0.3)
#p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
#p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "nobound_alpha_tau1_ft.svg"),
       plot = p, width=10, height=5, dpi=300)




# Save all Taus Cl
taus <- 1:100
for (t in taus) {
  summaryCl <- summarySE(cl[cl$tau == t,], c("count"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryCl, aes((1 - alpha), count))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
  p <- p + geom_errorbar(limits)
  p <- p + facet_grid(~ R, labeller = myLabeller)
  p <- p + ylim(0,100)
  p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  p <- p + xlab('Strength of Social Influence') + ylab('Number of Clusters')
  p <- p + myThemeMod + theme(strip.background = element_blank()) + ggtitle(title)
  #
  ggsave(filename=paste0(IMGPATH, "alpha/alpha_tau_cc_", sprintf("%04d", t), ".jpg"),
         plot = p)
}
# Movie
system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'alpha/alpha_tau_cc_%04d.jpg ',
              IMGPATH, 'alpha/movie_alpha_tau_cc.avi'))
# Save all Taus Ft
taus <- 1:100
for (t in taus) {
  summaryFt <- summarySE(cl[cl$tau == t,], c("fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
  p <- p + geom_errorbar(limitsFt)
  p <- p + facet_grid(~ R, labeller = myLabeller)                                        
  p <- p + xlab("Strength of social influence") + ylab('Avg. Distance from Truth')
  p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  p <- p + scale_fill_continuous(name="Distance\nfrom truth")
  p <- p + ylim(0, 6.3) + ggtitle(title)
  p <- p + myThemeMod + theme(strip.background = element_blank()) + ggtitle(title)
  #
  ggsave(filename=paste0(IMGPATH, "alpha/alpha_tau_", sprintf("%04d", t), ".jpg"),
         plot = p)
}
# Movie
system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'alpha/alpha_tau_%04d.jpg ',
              IMGPATH, 'alpha/movie_alpha_tau_ft.avi'))




# Alpha and different Taus
cl$alphabrk <- cut(cl$alpha, seq(0.01,0.99,0.01))

mycl <- cl[cl$t == 2000 &
           (cl$tau == 1 | cl$tau == 5 | cl$tau == 10 | cl$tau == 25 |
            cl$tau == 50 | cl$tau == 100),]

mycl <- cl[cl$t == 2000 & cl$tau == 1,]
summaryCl.tau1 <- summarySE(mycl, c("count"), c("alpha", "R", "tau"), na.rm=TRUE)

mycl <- cl[cl$t == 2000 & cl$tau == 50,]
summaryCl.tau50 <- summarySE(mycl, c("count"), c("alpha", "R", "tau"), na.rm=TRUE)
summaryCl.tau50$count2 <- summaryCl.tau50$count

mycl <- cl[cl$t == 2000 & (cl$tau == 50 | cl$tau == 1),]
summaryCl <- summarySE(mycl, c("count"), c("alpha", "R", "tau"), na.rm=TRUE)



colours <- c("#FC9272", "#FB6A4A", "#EF3B2C", "#9ECAE1", "#6BAED6", "#4292C6")


title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl.tau50, aes((1-alpha), count, group = as.factor(tau),
                                 fill = as.factor(tau)))
p <- p + geom_bar(stat = "identity", position = "dodge", width=0.01)
p <- p + geom_bar(data = summaryCl.tau1, stat = "identity",
                  position = "dodge", width=0.01)

#p <- p + geom_errorbar(limits, position = "dodge")
#p <- p + geom_errorbar(data = summaryCl.tau1, limits, position = "dodge")
p <- p + facet_grid(. ~ R, labeller = myLabeller)
p

#p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
#xlabText <- expression(paste('Strength of social influence (1-',alpha,')'))
#p <- p + xlab('Strength of Social Influence') + ylab('Number of Clusters')
#p <- p + myThemeMod + theme(strip.background = element_blank())
p <- p + scale_fill_manual(breaks = c(1,50),
                           palette = c(gg_color_hue2,
                             gg_color_hue3))
p


gg_color_hue <- function(n) {
  hues = seq(0, 100, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

gg_color_hue2 <- function(n) {
  hues = rep(260,n+1)
  luminances = seq(0,100, length=n+1)
  hcl(h=hues, l=65, c=luminances)[1:n]
}


gg_color_hue3 <- function(n) {
  hues = rep(0,n+1)
  luminances = seq(0,100, length=n+1)
  hcl(h=hues, l=65, c=luminances)[1:n]
}

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl.tau50, aes((1-alpha), count, fill = count))
p <- p + geom_bar(stat = "identity", position = "dodge", aes(group = as.factor(tau)))
p <- p + geom_errorbar(limits, position = "dodge")
p <- p + facet_grid(tau ~ R, labeller = myLabeller)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + xlab('Strength of Social Influence') + ylab('Number of Clusters')
p <- p + myThemeMod + theme(strip.background = element_blank())
p <- p + scale_fill_gradient(low = "red", high = "yellow", space = "Lab", na.value = "grey50", 
                             guide = "colourbar")
p

# TEST

set.seed(1234)
data <-
expand.grid(month = month.abb,
            building = c("Building A", "Building B", "Building C"),
            hc = c("Heating", "Cooling"))
data$value <- rnorm(nrow(data), 60, 10)

ggplot(data, aes(building,value,group=month)) + 
  geom_bar(stat = 'identity',
           position = 'dodge',
           aes(fill = interaction(building, hc)))

library("RColorBrewer")
colours <- c(brewer.pal(9,"Reds")[4:6], brewer.pal(9,"Blues")[4:6])

## NOISES ##
############

cl <- loadData(DUMPDIR, 'nobound_noises_tau/', 1)
  
clTau1 <- loadData(DUMPDIR, 'nobound_noises_tau1/', 1)

cl <- rbind(cl, clTau1[clTau1$alpha == 0.5,])

summaryCl <- summarySE(cl[cl$t == 2000 & cl$tau == 1,], c("count"), c("sigma", "epsilon", "R"), na.rm=TRUE)

title <- 'Cluster counts vs Angular noise and \nPosition noise'
p <- ggplot(summaryCl, aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Angular Noise') + ylab('Number of Clusters')
#p <- p + scale_x_continuous(labels = c("0", "0.02", "0.04", "0.06", "0.08", "0.1"),
#                            breaks = seq(0,0.1,0.02))
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
p <- p + myThemeMod + theme(strip.background = element_blank())
p


# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[22]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0"))
grob[["grobs"]][[23]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.1"))
grob[["grobs"]][[24]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.2"))
grob[["grobs"]][[25]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.3"))
grob[["grobs"]][[26]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.4"))
grob[["grobs"]][[27]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.5"))

svg(filename = paste0(IMGPATH, "nobound_noises_cc.svg"),
     width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()


# FT
summaryFt <- summarySE(cl[cl$t == 2000 & cl$tau == 1,], c("fromtruth.avg"), c("sigma", "epsilon", "R"), na.rm=TRUE)

title <- 'Distance from truth vs Angular noise and \nPosition noise'
p <- ggplot(summaryFt, aes(sigma, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
p <- p + xlab('Angular Noise') + ylab('Avg. Distance from Truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
#p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05", "0.075", "0.1"))
#p <- p + scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15))
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
p <- p + myThemeMod +  theme(strip.background = element_blank())
#p <- p + theme(axis.text.y = element_text(size=18))
p

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[22]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0"))
grob[["grobs"]][[23]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.1"))
grob[["grobs"]][[24]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.2"))
grob[["grobs"]][[25]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.3"))
grob[["grobs"]][[26]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.4"))
grob[["grobs"]][[27]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.5"))

svg(filename = paste0(IMGPATH, "nobound_noises_ft.svg"),
     width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()



# Save all Taus Cl
taus <- 1:100
for (t in taus) {
  summaryCl <- summarySE(cl[cl$tau == t,], c("count"), c("sigma", "epsilon", "R"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryCl, aes(sigma, count))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
  p <- p + geom_errorbar(limits)
  p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
  p <- p + xlab('Angular Noise') + ylab('Number of Clusters')
  p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
  p <- p + myThemeMod + theme(strip.background = element_blank())
  p <- p + ggtitle(title) + ylim(0,95)
  #
  ggsave(filename=paste0(IMGPATH, "noises/noises_tau_cc_", sprintf("%04d", t), ".jpg"),
         plot = p)
}

system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'noises/noises_tau_cc_%04d.jpg ',
              IMGPATH, 'noises/movie_noises_tau_cc.avi'))

# Save all Taus Ft
taus <- 1:100
for (t in taus) {
  summaryFt <- summarySE(cl[cl$tau == t,], c("fromtruth.avg"), c("sigma", "epsilon", "R"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryFt, aes(sigma, fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
  p <- p + geom_errorbar(limitsFt)
  p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
  p <- p + xlab('Angular Noise') + ylab('Avg. Distance from Truth')
  p <- p + scale_fill_continuous(name="Distance\nfrom truth")
  p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
  p <- p + myThemeMod +  theme(strip.background = element_blank())
  p <- p + ggtitle(title) + ylim(0,6)
  #
  ggsave(filename=paste0(IMGPATH, "noises/noises_tau_ft_", sprintf("%04d", t), ".jpg"),
         plot = p)
}

system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'noises/noises_tau_ft_%04d.jpg ',
              IMGPATH, 'noises/movie_noises_tau_ft.avi'))

## TAU ##
#########

# Using Alpha
cl <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1)
clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)
cl <- rbind(cl, clTau1)

# CL
summaryCl <- summarySE(cl[cl$t == 2000 & (cl$alpha == 0.5 | cl$alpha == 0.01 | cl$alpha == 0.99),], c("count"), c("tau", "R", "alpha"), na.rm=TRUE)

title <- 'Cluster counts by strength of the truth'
p <- ggplot(summaryCl, aes((100 - tau), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=1))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(alpha ~  R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + myThemeMod + theme(strip.background = element_blank())
p


# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[13]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.01"))
grob[["grobs"]][[14]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.5"))
grob[["grobs"]][[15]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.99"))

svg(filename = paste0(IMGPATH, "nobound_tau_cc.svg"),
    width=9)
grid.newpage()
grid.draw(grob)
dev.off()

summaryFt <- summarySE(cl[cl$t == 2000 & (cl$alpha == 0.5 | cl$alpha == 0.01 | cl$alpha == 0.99),], c("fromtruth.avg"), c("tau", "R", "alpha"), na.rm=TRUE)


# FT
title <- 'Distance from truth vs Strength of the truth\'s signal'
p <- ggplot(summaryFt, aes((100 - tau), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=1))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(alpha ~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod + theme(strip.background = element_blank())
p

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[13]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.01"))
grob[["grobs"]][[14]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.5"))
grob[["grobs"]][[15]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.99"))

svg(filename = paste0(IMGPATH, "nobound_tau_ft.svg"),
    width=9)
grid.newpage()
grid.draw(grob)
dev.off()


## VSCALING ##
##############

cl <- loadData(DUMPDIR, 'nobound_vscale_tau/', 1)
#cl <- loadData(DUMPDIR, 'final_vscaling/')

cl$vbr <- cut(cl$init.vscaling,  breaks=c(0,0.5,1, 1.5,2,10))
cl$vbr <- as.factor(cl$vbr)

summaryCl <- summarySE(cl[cl$t == 2000 & cl$tau == 1,], c("count"),
                       c("init.vscaling", "R", "alpha"), na.rm=TRUE)

title <- 'Cluster counts by initial velocity'
p <- ggplot(summaryCl, aes(as.factor(init.vscaling), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(alpha ~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Cluster counts')
p <- p + ggtitle(title) + myThemeMod + theme(legend.position = "none",
                                             strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "nobound_vscale_tau1_cc.svg"), plot = p, width = 9)

summaryFt <- summarySE(cl[cl$t == 2000 & cl$tau == 1,], c("fromtruth.avg"),
                       c("init.vscaling", "R", "alpha"), na.rm=TRUE)

title <- 'Distance from truth by initial velocity'
p <- ggplot(summaryFt, aes(as.factor(init.vscaling), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(alpha ~ R, labeller = myLabeller)
p <- p + xlab('Initial Velocity') + ylab('Distance from truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth") + theme(legend.position = "none") 
p <- p + ggtitle(title)
p


ggsave(filename = paste0(IMGPATH, "nobound_vscale_tau1_ft.svg"),
       plot = p, width=9)



# Save all Taus Cl
taus <- 1:100
for (t in taus) {
  summaryCl <- summarySE(cl[cl$tau == t,], c("count"), c("init.vscaling", "R", "alpha"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryCl, aes(as.factor(init.vscaling), count))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
  p <- p + geom_errorbar(limits)
  p <- p + facet_grid(alpha ~ R, labeller = myLabeller)
  p <- p + xlab('Velocity Multiplier') + ylab('Number of Clusters')
#  p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
  p <- p + myThemeMod + theme(strip.background = element_blank())
  p <- p + ggtitle(title) + ylim(0,100)
  #
  ggsave(filename=paste0(IMGPATH, "vscale/vscale_tau_cc_", sprintf("%04d", t), ".jpg"),
         plot = p, width = 9)
}
system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'vscale/vscale_tau_cc_%04d.jpg ',
              IMGPATH, 'vscale/movie_vscale_tau_cc.avi'))

# Save all Taus Ft
taus <- 1:100
for (t in taus) {
  summaryFt <- summarySE(cl[cl$tau == t,], c("fromtruth.avg"), c("init.vscaling", "R", "alpha"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryFt, aes(as.factor(init.vscaling), fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
  p <- p + geom_errorbar(limitsFt)
  p <- p + facet_grid(alpha ~ R, labeller = myLabeller)
  p <- p + xlab('Velocity Multiplier') + ylab('Avg. Distance from Truth')
  p <- p + scale_fill_continuous(name="Distance\nfrom truth")
#  p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
  p <- p + myThemeMod +  theme(strip.background = element_blank())
  p <- p + ggtitle(title) + ylim(0,60)
  #
  ggsave(filename=paste0(IMGPATH, "vscale/vscale_tau_ft_", sprintf("%04d", t), ".jpg"),
         plot = p, width = 9)
}
system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'vscale/vscale_tau_ft_%04d.jpg ',
              IMGPATH, 'vscale/movie_vscale_tau_ft.avi'))
