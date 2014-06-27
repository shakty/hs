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
    value[value == 0.03] <- "Small radius (R = 0.03)"
    value[value == 0.3] <- "Large Radius (R = 0.3)"
  } else if (var == "clbr") {
    value[value == "(0,1]"] <- "1\ncluster"
    value[value == "(1,5]"] <- "2-5\nclusters"
    value[value == "(5,10]"] <- "6-10\nclusters"
    value[value == "(10,20]"] <- "11-20\nclusters"
    value[value == "(20,30]"] <- "21-30\nclusters"
  }
  return(value)
}

nClustersLabeller <- function(var, value){
  value <- as.character(value)
   
  return(value)
}

myThemeMod <- theme(legend.position = "none",
                    axis.title.x = element_text(vjust=-1, size=24),
                    axis.title.y = element_text(vjust=-0.1, size=24),
                    plot.margin=unit(c(10,10,10,10),"mm"),                    
                    plot.title = element_text(vjust=3, size=24,face="bold")
                    )

limits <- aes(ymax = count + se, ymin = count - se)
limitsFt <- aes(ymax = fromtruth.avg + se, ymin = fromtruth.avg - se)
limitsCI <- aes(ymax = count + ci, ymin = count - ci)

theme_set(theme_bw(base_size = 30))
theme_white()

XINTERCEPT <- 0.15

## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/NOBOUND/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

## R ##
#######
cl <- loadData(DUMPDIR, 'nobound_R_tau/', 1)

clTau1 <- loadData(DUMPDIR, 'nobound_R_tau1/', 1)

cl <- rbind(cl, clTau1[clTau1$alpha == 0.5,])

# Only TAU = 1
cl <- clTau1

# CL
summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("R"), na.rm=TRUE)


p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
xlabText <- expression(paste('Radius of Influence ',R))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
p <- p  + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "nobound_R_tau1_cc.svg"),
       plot = p, width=10, height=5, dpi=300)


# FT
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
xlabText <- expression(paste('Radius of Influence ',R))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod 
p

ggsave(filename = paste0(IMGPATH, "nobound_R_tau1_ft.svg"),
       plot = p, width=10, height=5, dpi=300)



# Save all Taus Cl
taus <- 1:100
for (t in taus) {
  summaryCl <- summarySE(cl[cl$tau == t,], c("count"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryCl, aes(R, count))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
  p <- p + geom_errorbar(limits)
  p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
  p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
  p <- p + xlab('Radius of Influence') + ylab('Avg. Number of Clusters')
  p <- p  + myThemeMod + ggtitle(title) + ylim(0,100)
  p
  ggsave(filename=paste0(IMGPATH, "R/R_tau_cc_", sprintf("%04d", t), ".jpg"),
         plot = p)
}
system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'R/R_tau_cc_%04d.jpg ',
              IMGPATH, 'R/movie_R_tau_cc.avi'))
# Save all Taus Ft
taus <- 1:100
for (t in taus) {
  summaryFt <- summarySE(cl[cl$tau == t,], c("fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("Tau: ", t)
  p <- ggplot(summaryFt, aes(R, fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
  p <- p + geom_errorbar(limitsFt)
  p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
  p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
  p <- p + xlab('Radius of Influence') + ylab('Avg. Distance from Truth')
  p <- p + scale_fill_continuous(name="Distance\nfrom truth")
  p <- p + myThemeMod + ggtitle(title) + ylim(0,6.3)
  #
  ggsave(filename=paste0(IMGPATH, "R/R_tau_", sprintf("%04d", t), ".jpg"),
         plot = p)
}
system(paste0('ffmpeg -qscale 1 -r 2 -b 9600 -y -i ',
              IMGPATH, 'R/R_tau_%04d.jpg ',
              IMGPATH, 'R/movie_R_tau_ft.avi'))


## ALPHA ##
###########

cl <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1)
clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)
cl <- rbind(cl, clTau1)
# Only tau 1
cl <- clTau1

# CL
summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("alpha", "R"), na.rm=TRUE)

p <- ggplot(summaryCl, aes((1 - alpha), count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(~ R, labeller = myLabeller)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of Social Influence ',alpha))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
p <- p + myThemeMod + theme(strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "nobound_alpha_tau1_cc.svg"),
       plot = p, width=10, height=5, dpi=300)

# FT
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(~ R, labeller = myLabeller)
xlabText <- expression(paste('Strength of Social Influence ',alpha))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
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
# Only Tau 1
cl <- clTau1

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("sigma", "epsilon", "R"), na.rm=TRUE)

p <- ggplot(summaryCl, aes(sigma, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
xlabText <- expression(paste('Angular Noise ',sigma))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
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
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("sigma", "epsilon", "R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(sigma, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(epsilon ~ R, labeller = myLabeller)
xlabText <- expression(paste('Angular Noise ',sigma))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
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



#### SPEEDTEST
##############

library(texreg)
library(mediation)
library(grid)
library(lattice)
library(rgl)

DIR <- 'clusters_vs_progress_nobound/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)

data <- read.table('speedtest.csv', head = T, sep = ",")

# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
#data$init.placement <- as.factor(data$init.placement)

# If consensus is not reached it has value -1. Replace with NA
data[ , 15:28 ][ data[ , 15:28 ] == -1 ] <- NA
# 20.000
data[ , 15:28 ][ data[ , 15:28 ] == -1 ] <- 20000

# If no consensus was reached replace ccounts with NA
data[ , 29:40 ][ data[ , 29:40 ] == -1 ] <- NA

# Replaces everything, not good
#data[] <- lapply(data, function(x){replace(x, x == -1, 20000)})

data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)
# Add this if wanna see by cluster as well.
data$clbr <- cut(data$init.ccount, breaks=c(0,1,5, 10, 20, 30))


###########################
# Select ALPHA
###########################

dataAlpha <- data[data$alpha == 0.01,]
dataAlpha <- data[data$alpha == 0.5,]
dataAlpha <- data[data$alpha == 0.99,]
dataAlpha <- data

############
# Divide Tau
############

data1 <- dataAlpha[dataAlpha$tau == 1,]
data50 <- dataAlpha[dataAlpha$tau == 50,]

###########################
# Select Dist: ALL / vs 0.5
###########################
# TAU = 1
mydata <- data1[data1$init.placement == 1,]
mydata <- data1
# TAU = 50
mydata <- data50[data50$init.placement == 1,]
mydata <- data50
##################

### Exclude init.placement < 0.2
mydata <- mydata[mydata$init.placement >= 0.2,]


mydata <- mydata[mydata$init.placement %in% seq(0.2,0.9,0.1),]


p <- ggplot(mydata, aes(init.ccount, consensus75, color = as.factor(init.placement)))
#p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm")
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
#p <- p + myThemeMod
p <- p + facet_grid(.~ R, labeller = myLabeller)
p


ggsave(filename = paste0(IMGPATH, "SPEEDTEST/st_alpha001tau1.jpg"),
       width=10, height=20, dpi=300)



title <- "Number of clusters vs progress and time to reach a consensus"
p <- ggplot(mydata, aes(init.ccount, consensus75, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(alpha  ~ R, labeller = myLabeller)
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST/st_alpha001tau1_all.jpg"),
       width=10, height=12, dpi=300)




# Add this if wanna see by cluster as well.
mydata$clbr <- cut(mydata$init.ccount, breaks=c(0,1,5, 10, 20, 30))

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
  summaryNew <- summarySE(mydata[mydata$alpha != 0.99,], c(myVar), c('R', 'clbr'), na.rm=TRUE)
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


summaryData$R <- as.factor(summaryData$R)

summaryDataR003 <- summaryData[summaryData$R == 0.03,]
summaryDataR03 <- summaryData[summaryData$R == 0.3,]

title <- 'Temporal evolution of consensus building \n by number of initial clusters'
p <- ggplot(summaryDataR003, aes(x = share, weight=time, fill=R))
p <- p + geom_bar(aes(y=time), stat = "identity")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se), width = 5)
p <- p + geom_bar(data = summaryDataR03, aes(y=time), stat = "identity")
p <- p + geom_errorbar(data = summaryDataR03, aes(ymax = time + se, ymin = time - se), width=5)
p <- p + xlab('Share of Consensus') + ylab('Avg. Time Passed')
# Remove clbr if not needed
p <- p + facet_grid(.~clbr, labeller = myLabeller)
p <- p + scale_fill_discrete(name="Radius of\nInfluence")
p  + myThemeMod + theme(legend.background = element_rect(fill = "white",
                          colour = "grey"),
                        legend.position = c(1.1, 0.5),
                        plot.margin = unit(c(10,30,10,10),"mm"),
                        axis.text.x = element_text(size=15),
                        strip.background = element_blank())
p <- p + guides(col=guide_legend(ncol=2))

ggsave(filename = paste0(IMGPATH, "SPEEDTEST/st_alpha001tau1_temporal_evo.svg"),
       width=12, height=8, dpi=300)


mydata.summary <- summarySE(mydata, "consensus75", c("R", "alpha", "init.ccount", "init.placement"), na.rm = TRUE)

mydata.summary <- mydata.summary[mydata.summary$init.placement %in% seq(0.2,1,0.1),]

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary, aes(init.ccount, consensus75))
p <- p + geom_point(aes(color = as.factor(init.placement)), size=3)
p <- p + geom_line(aes(color = as.factor(init.placement)), alpha=0.5)
p <- p + xlab("Initial Number of Clusters") + ylab("Avg. Time to Consensus")
p <- p + facet_grid(alpha~R, labeller=myLabeller)
#p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\nDistance\nfrom Truth")
p <- p + myThemeMod + theme(legend.position = c(0.8, 0.3),
                            legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(vjust=3, size=16,face="bold"),
                            legend.text = element_text(size=14),
                            legend.key.width = unit(1.5, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA),
                            strip.background = element_blank()
                            )
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST/st_tau1_consensus_by_progress.svg"),
       width=10, height=7, dpi=300)

mydata003 <- mydata[mydata$R == 0.03,]
mydata003$moreThan9 <- mydata003$init.ccount > 9

title <- "Only R = 0.03. Smoothed functions."
p <- ggplot(mydata003, aes(init.ccount, consensus75))
p <- p + geom_smooth(alpha=0.5, aes(color = as.factor(init.placement)))
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p


title <- "Only R = 0.03. More than 9 clusters?"
p <- ggplot(mydata003, aes(init.ccount, consensus75))
p <- p + geom_smooth(alpha=0.5, method="lm", aes(color = as.factor(init.placement)))
p <- p + xlab("Initial number of clusters") + ylab("Time to Consensus")
p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p <- p + facet_grid( ~ moreThan9, labeller = myLabeller)
p



## NEW PLOTS


mydata.summary <- summarySE(mydata, "ccount25", c("init.ccount", "init.placement"), na.rm = TRUE)

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary, aes(init.ccount, ccount25))
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

mydata.summary <- summarySE(mydata, "ccount30", c("R", "init.ccount", "init.placement"), na.rm = TRUE)

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary[mydata.summary$R == 0.3,], aes(as.factor(init.placement), ccount30))
p <- p + geom_boxplot(aes(color = as.factor(init.placement)), size=4)
#p <- p + geom_line(aes(color = as.factor(init.placement)))
p <- p + xlab("Initial Distance from Truth") + ylab("Number of Clusters at 25")
#p <- p + facet_grid(alpha~.)
#p <- p + ggtitle(title)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p <- p +  theme(
                legend.background = element_rect(fill = "white", colour = "grey"),
                legend.title = element_text(vjust=3, size=16,face="bold"),
                legend.text = element_text(size=14),
                legend.key.width = unit(1.5, "cm"),
                legend.key =  element_rect(fill = "white", colour = NA)
                )
p

## Regressions

fit001 <- lm(consensus75 ~ ccount30, data = mydata[mydata$alpha == 0.01,])
fit05 <- lm(consensus75 ~ ccount30, data = mydata[mydata$alpha == 0.5,])
fit099 <- lm(consensus75 ~ ccount30, data = mydata[mydata$alpha == 0.99,])
summary(fit099)

plot(fit099)

texreg(list(fit001,fit05,fit099))


#########################
### Clusters vs Progress: RBANDS
#########################

DIR <- 'clusters_vs_progress_nobound_rbands/'

PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
## IMG DIR

IMGPATH <- paste0(DUMPDIR, "imgs/NOBOUND/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

data <- read.table('speedtest.csv', head = T, sep = ",")
# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
data$init.placement <- as.factor(data$init.placement)
data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)

# If consensus is not reached it has value -1. Replace with NA
data[ , 15:28 ][ data[ , 15:28 ] == -1 ] <- 20000

# If no consensus was reached replace ccounts with NA
data[ , 29:40 ][ data[ , 29:40 ] == -1 ] <- NA

# If consensus is not reached it has value -1. Replace with NA
#data[] <- lapply(data, function(x){replace(x, x == -1, NA)})

data$band <- as.factor(data$init.band.i)
data$bandbr <- cut(data$init.band.i, breaks=c(0,0.3, 0.5, 0.7, 0.9, 1))

# SELECT ALPHA
##############
dataAlpha <- data[data$alpha == 0.01,]
dataAlpha <- data[data$alpha == 0.5,]
dataAlpha <- data[data$alpha == 0.99,]
dataAlpha <- data


dataNA <- data[!complete.cases(data),]

dataNA[] <- lapply(dataNA, function(x){replace(x, x == -1, NA)})


# Converting all consensus shares and all clusters counts to two columns
mydata <- dataAlpha
colnamesData <- colnames(mydata)

data.subset <- subset(mydata, select=colnamesData[seq(1,27)])
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
                      data.melted$consensus.share != 75
                      ,]

mydata <- summarySE(mydata, "ccount", c("consensus.share", "init.band.i", "R"), na.rm = TRUE)

title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount, color = consensus.share))
p <- p + geom_point(size=5)
p <- p + geom_line(size=2, alpha=0.5)
p <- p + geom_errorbar(width=0.02, aes(ymin=ccount-ci,ymax=ccount+ci))
#p <- p + geom_smooth(alpha=0.5, color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + ggtitle(title) # + myThemeMod
p <- p + scale_y_discrete(breaks=seq(1,20,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p

title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount, color = consensus.share))
p <- p + geom_smooth(method="lm", se=FALSE, size=2)
p <- p + xlab("Initial Distance from Truth") + ylab("Avg. Number of Clusters")
p <- p + myThemeMod
p <- p + scale_y_discrete(breaks=seq(1,20,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p <- p + scale_color_hue(name="Share of\nconsensus")
p <- p +  theme(legend.position = c(0.3, 0.9), legend.direction = "horizontal",
                legend.background = element_rect(fill = "white", colour = "grey"),
                legend.title = element_text(vjust=3, size=16,face="bold"),
                legend.text = element_text(size=14),
                legend.key.width = unit(1.5, "cm"),
                strip.background = element_blank(),
                legend.key =  element_rect(fill = "white", colour = NA)
                )
p <- p + guides(col=guide_legend(ncol=4))
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST_RBANDS/clusters_by_progress_lm.jpg"),
       width=12, height=10, dpi=300)

title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount, color = consensus.share))
p <- p + geom_point(size=5)
p <- p + geom_line(size=2, alpha=0.5)
p <- p + geom_errorbar(width=0.02, aes(ymin=ccount-ci,ymax=ccount+ci))
#p <- p + geom_smooth(alpha=0.5, color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + ggtitle(title) # + myThemeMod
p <- p + scale_y_discrete(breaks=seq(1,20,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p


mydata = data[data$alpha == 0.5,]
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
  summaryNew <- summarySE(mydata[mydata$alpha != 0.99,], c(myVar), c('R', 'bandbr'), na.rm=TRUE)
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
                       
title <- 'Temporal evolution of consensus building \n by initial distance from truth'
p <- ggplot(summaryData, aes(x = share, weight=time, group=R, fill=R))
p <- p + geom_bar(aes(y=time), stat = "identity", position = "dodge")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se), position = "dodge")
p <- p + xlab('AvgShare of Consensus') + ylab('Avg. Time Passed')
# Comment bandbr if not needed
p <- p + facet_grid(.~bandbr)
p <- p + ggtitle(title)  + myThemeMod + theme(legend.position = c(1.1, 0.5),
                                              plot.margin = unit(c(10,30,10,10),"mm"),
                                              axis.text.x = element_text(size=15),
                                              strip.background = element_blank())
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST_RBANDS/temporal_evo_consensus_by_initial_distance.jpg"),
       width=10, height=8, dpi=300)


mydata <- data[data$alpha == 0.5,]
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount25, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,20,1))
p <- p + facet_grid(.  ~ R, labeller = myLabeller)
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST_RBANDS/clusters_vs_progress.svg"),
       width=10, height=10, dpi=300)

mydata <- data[data$alpha == 0.5,]
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(as.factor(init.band.i), ccount50, color = R))
p <- p + geom_boxplot(size=2)
#p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + scale_y_discrete(breaks=seq(1,15,1))
p <- p + facet_grid(R ~ ., labeller = myLabeller)
p <- p + theme(axis.text.x = element_text(angle = 90))
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST_RBANDS/clusters_vs_progress_boxplot.svg"),
       width=10, height=10, dpi=300)


mydata <- data
title <- "Number of clusters vs progress"
p <- ggplot(mydata, aes(init.band.i, ccount50, color = R))
p <- p + geom_jitter(size=2)
p <- p + geom_smooth(alpha=0.5, method="lm", color="black")
p <- p + xlab("Initial distance from truth") + ylab("Number of clusters")
p <- p + myThemeMod + ggtitle(title)
p <- p + facet_grid(alpha  ~ R, labeller = myLabeller)
p <- p + theme(axis.text.x = element_text(angle = 90))
p

ggsave(filename = paste0(IMGPATH, "SPEEDTEST_RBANDS/clusters_vs_progress_by_alpha.svg"),
       width=10, height=10, dpi=300)



## TAU 2 ##

cl <- loadData(DUMPDIR, 'scan_tau_again2/')

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("tau", "R", "boundaries"), na.rm=TRUE)

title <- 'Cluster counts by strength of the truth'
p <- ggplot(summaryCl[summaryCl$boundaries == 1,], aes(tau, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
p <- p + xlab('Truth strength in percentage') + ylab('Cluster counts')
p <- p + myThemeMod + theme(strip.background = element_blank())
p

# To be taken from TAU 20000
#ggsave(filename = paste0(IMGPATH, "scan_tau_upto_10.jpg"),
#       plot = p, width=10, height=5, dpi=300)

## R 2 ##

cl <- loadData(DUMPDIR, 'scan_R_again2/')

clTau1 <- loadData(DUMPDIR, 'nobound_R_tau1/', 1)
clTau1$boundaries <- 0

cl <- clTau1

cl <- rbind(cl, clTau1[clTau1$alpha == 0.5,])

# To be taken from TAU 20000
#ggsave(filename = paste0(IMGPATH, "scan_tau_upto_10.jpg"),
#       plot = p, width=10, height=5, dpi=300)

# CL
summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("R"), na.rm=TRUE)


p <- ggplot(summaryCl, aes(R, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
p <- p + geom_errorbar(limits)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 22, label = "Convergence Zone", size=8)
xlabText <- expression(paste('Radius of Influence ',R))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
p <- p  + myThemeMod
p

ggsave(filename = paste0(IMGPATH, "nobound_R_tau1_cc.svg"),
       plot = p, width=10, height=5, dpi=300)


# FT
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(R, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + geom_vline(xintercept = XINTERCEPT, colour="red", linetype = "longdash", size = 1)
p <- p + annotate("text", x = 0.55, y = 0.22, label = "Convergence Zone", size=8)
xlabText <- expression(paste('Radius of Influence ',R))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_fill_continuous(name="Distance\nfrom truth")
p <- p + myThemeMod 
p

ggsave(filename = paste0(IMGPATH, "nobound_R_tau1_ft.svg"),
       plot = p, width=10, height=5, dpi=300)

