# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")
source("/opt/MATLAB_WORKSPACE/hs/R/init.NOBOUND.R")
library(scales)

# DUMPDIR 
DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'
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

ggsave(filename = paste0(IMGPATH, "nobound_alpha_tau1_cc_TITLE_FIX.svg"),
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
data.old <- read.table('speedtest.csv', head = T, sep = ",")
#data <- data.old[data.old$init.placement < 0.5,]
DIR <- 'clusters_vs_progress_nobound_again/'
PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
data.new <- read.table('speedtest.csv', head = T, sep = ",")
data.old.subset <- data.old[data.old$tau == 1 & data.old$init.placement < 0.5,]
data <- rbind(data.old.subset, data.new)


## Mixin in
dd <- data
##

DIR <- 'clusters_vs_progress_nobound_biggap/'
PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
data <- read.table('speedtest.csv', head = T, sep = ",")

## Mixin in
dd.new <- data
data <- rbind(dd[dd$R == 0.03,], dd.new)
##

DIR <- 'clusters_vs_progress_nobound_biggap003/'
PATH <- paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
data <- read.table('speedtest.csv', head = T, sep = ",")


#data <- data.new

# Replace -1 in init.placement with 0
data$init.placement[data$init.placement == -1] <- 0
#data$init.placement <- as.factor(data$init.placement)
# If consensus is not reached it has value -1. Replace with NA
#data[ , 15:28 ][ data[ , 15:28 ] == -1 ] <- NA
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

mydata <- data1

mydata <- data1[data1$init.placement == 1,]

# TAU = 50
mydata <- data50[data50$init.placement == 1,]
mydata <- data50
##################

### Exclude init.placement < 0.2
mydata <- mydata[mydata$init.placement >= 0.2,]

### Adapta alpha to formulation in the paper
mydata$alpha2 <- mydata$alpha
alpha001 <- mydata$alpha == 0.01
alpha099 <- mydata$alpha == 0.99
mydata[alpha001,]$alpha2 <- 0.99
mydata[alpha099,]$alpha2 <- 0.01

limitPlacement <- sort(unique(mydata$init.placement))[2]

#mydata <- mydata[mydata$init.placement %in% seq(0.2,0.9,0.1),]



p <- ggplot(mydata[mydata$init.placement == 1.5892,], aes(init.ccount, consensus75, color = as.factor(init.placement)))
p <- p + geom_jitter(alpha = 0.5, color="grey")
p <- p + geom_smooth(alpha=0.5, method="lm", size = 2)
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
  # A
  summaryNew <- summarySE(mydata[mydata$alpha != 0.99,], c(myVar), c('R', 'clbr'), na.rm=TRUE)
  # B
  #summaryNew <- summarySE(mydata[mydata$alpha == 0.01,], c(myVar), c('R', 'init.ccount', 'init.placement'), na.rm=TRUE)
  # C
  #summaryNew <- summarySE(mydata[mydata$alpha != 0.99 & mydata$init.placement > limitPlacement,], c(myVar), c('R', 'init.ccount'), na.rm=TRUE)
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

summaryData$R <- factor(summaryData$R, levels = c("0.03","0.3"))
summaryDataR003 <- summaryData[summaryData$R == 0.03,]
summaryDataR03 <- summaryData[summaryData$R == 0.3,]


# To produce use A (see above)
title <- 'Temporal evolution of consensus building \n by number of initial clusters'
p <- ggplot(summaryDataR003, aes(x = share, weight=time, fill=R))
p <- p + geom_bar(aes(y=time), stat = "identity")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se), width = 5)
#p <- p + geom_bar(data = summaryDataR03, aes(y=time), stat = "identity")
#p <- p + geom_errorbar(data = summaryDataR03, aes(ymax = time + se, ymin = time - se), width=5)
p <- p + xlab('Share of Consensus') + ylab('Avg. Time Passed')
# Remove clbr if not needed
p <- p + facet_grid(.~clbr, labeller = myLabeller)
p <- p + scale_fill_discrete(name="Radius of\nInfluence")
#p <- p + scale_x_discrete(labels=(0,0.25,0.5,0.75,1))
p <- p  + myThemeMod + theme(
                        legend.background = element_rect(fill=alpha('white', 0.4)),
                        legend.position = c(0.2, 0.8),
                        plot.margin = unit(c(10,30,10,10),"mm"),
                        axis.text.x = element_text(size=15),
                        strip.background = element_blank())
p <- p + guides(col=guide_legend(ncol=2))
p


ggsave(filename = paste0(IMGPATH, "SPEEDTEST/st_alpha001tau1_temporal_evo_B.svg"),
       width=12, height=6, dpi=300)


# To produce use B (see above)
title <- 'Temporal evolution of consensus building \n by number of initial clusters'
p <- ggplot(summaryData, aes(x = share, weight=time, fill=R))
p <- p + geom_bar(aes(y=time), stat = "identity")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se), width = 5)
p <- p + xlab('Share of Consensus') + ylab('Avg. Time Passed')
# Remove clbr if not needed
p <- p + facet_wrap(~init.ccount)
p <- p + scale_fill_discrete(name="Radius of\nInfluence")
p  + myThemeMod + theme(legend.background = element_rect(fill = "white",
                          colour = "grey"),
                        legend.position = c(1.1, 0.5),
                        plot.margin = unit(c(10,30,10,10),"mm"),
                        axis.text.x = element_text(size=15),
                        strip.background = element_blank())

# To produce use B (see above)
# Reversing init.ccount and share
SHARE = 90
R003 <- summaryDataR003[summaryDataR003$share == SHARE & summaryDataR003$init.placement < 0.3,]
R03 <- summaryDataR03[summaryDataR03$share == SHARE & summaryDataR003$init.placement %in% c(0.7,0.8,0.9,1),]

R03 <- summaryDataR03[summaryDataR03$share == SHARE,]

p <- ggplot(R03, aes(x = init.ccount, weight=time))
p <- p + geom_bar(aes(y=time), stat = "identity", position="dodge")
p <- p + geom_errorbar(aes(ymax = time + se, ymin = time - se))
p <- p + geom_smooth(aes(y=time), size = 1.5, color = "red")
p <- p + xlab('Number of Initial Clusters') + ylab('Avg. Time Passed')
# Remove clbr if not needed
#p <- p + facet_wrap(~init.placement)
p <- p + scale_fill_discrete(name="Radius of\nInfluence")
p  + myThemeMod + theme(legend.background = element_rect(fill = "white",
                          colour = "grey"),
                        legend.position = c(1.1, 0.5),
                        plot.margin = unit(c(10,30,10,10),"mm"),
                        axis.text.x = element_text(size=15),
                        strip.background = element_blank())


ggsave(filename = paste0(IMGPATH, "SPEEDTEST/consensus_wave_R03.svg"),
       width=8, height=8, dpi=300)


# SMALL R 3 cases
p <- ggplot(R003[R003$init.placement %in% c(0.1,0.15,0.5),], aes(x = init.ccount, time, color=init.placement))
p <- p + geom_point(aes(ymax = time + se, ymin = time - se), alpha = 0.5)
p <- p + geom_smooth(aes(group=init.placement))
p <- p + xlab('Number of Initial Clusters') + ylab('Avg. Time Passed')
#p <- p + facet_wrap(~init.placement)
p <- p + scale_fill_discrete(name="Radius of\nInfluence")
p  + myThemeMod + theme(legend.background = element_rect(fill = "white",
                          colour = "grey"),
                        legend.position = c(1, 0.5),
                        plot.margin = unit(c(10,30,10,10),"mm"),
                        axis.text.x = element_text(size=15),
                        strip.background = element_blank())

R03b <- R03[R03$init.placement %in% c(0.7,0.9,1),]
R003b <- R003[R003$init.placement %in% c(0.1,0.15,0.2,1),]

p <- ggplot(R03b, aes(x = init.ccount, time, color=as.factor(init.placement)))
p <- p + geom_point(alpha = 0.5)
p <- p + geom_smooth(aes(group=init.placement), se = FALSE, size=2)
p <- p + xlab('Number of Initial Clusters') + ylab('Avg. Time Passed')
#p <- p + facet_wrap(~init.placement)
p <- p + scale_color_discrete(name="Initial Distance\nFrom Truth")
p  + myThemeMod + theme(legend.background = element_rect(fill=alpha('white', 0.4)),
                        legend.position = c(0.3, 0.8),
                        plot.margin = unit(c(10,30,10,10),"mm"),
                        legend.key.width = unit(1,"cm"),
                        legend.key.height = unit(1,"cm"),
                        axis.text.x = element_text(size=15),
                        strip.background = element_blank())


mydata.summary <- summarySE(mydata, "consensus75", c("R", "alpha2", "init.ccount", "init.placement"), na.rm = TRUE)
#mydata.summary <- mydata.summary[mydata.summary$init.placement %in% seq(0.2,1,0.1),]

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary, aes(init.ccount, consensus75))
p <- p + geom_jitter(aes(color = init.placement), size=2)
#p <- p + geom_line(aes(color = init.placement), alpha=0.5)
p <- p + xlab("Initial Number of Clusters") + ylab("Avg. Time to Consensus")
p <- p + geom_smooth(color="red", alpha = 0.2, method="lm")
p <- p + facet_grid(alpha2~R, labeller=myLabeller)
p <- p + scale_color_gradient(name="Initial Distance\nfrom Truth")
p <- p + myThemeMod + theme(legend.position = "right", #c(0.15, 0.80),
                            #legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(size=18,face="bold"),
                            legend.text = element_text(size=16),
                            legend.key.height = unit(0.7, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA),
                            strip.background = element_blank()
                            )
p


mydata.subset <- mydata.summary[mydata.summary$R == 0.03 & mydata.summary$init.placement %in% c(0.1,0.15,0.2,0.5) & mydata.summary$init.ccount < 30 & mydata.summary$alpha == 0.99,]

mydata.subset <- mydata.summary[mydata.summary$init.placement %in% c(0.8,0.85,0.9,0.95,1) & mydata.summary$init.ccount < 16,]

mydata.subset <- mydata.summary[mydata.summary$init.placement %in% c(0.1,0.15,0.2) & mydata.summary$init.ccount < 16,]

mydata.subset <- mydata.summary[mydata.summary$init.ccount < 16,]
extremes <- c(min(mydata$init.placement),max(mydata$init.placement))
mydata.subset011 <- mydata.summary[mydata.summary$init.ccount < 16 & mydata.summary$init.placement %in% extremes,]

mydata.subset011b <- mydata.summary[mydata.summary$init.placement %in% c(0.1,0.5,0.7, 0.9,1),]


# high social influence, small R, c < 16, d = 1
fit99 <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.01 & mydata$R == 0.03 & mydata$init.placement == 1,]); summary(fit99)

# Medium social influence, small R, c < 16, d = 1
fit5 <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.5 & mydata$R == 0.03 & mydata$init.placement == 1,]); summary(fit5)

# High social influence, big R, c < 16, d = 1
fit99b <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.01 & mydata$R == 0.3 & mydata$init.placement == 1,]); summary(fit99b)

# Medium social influence, big R, c < 16, d = 1
fit5b <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.5 & mydata$R == 0.3 & mydata$init.placement == 1,]); summary(fit5b)

# Low social influence, big R, c < 16, d = 1
fit01b <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.99 & mydata$R == 0.3 & mydata$init.placement == 1,]); summary(fit01b)


title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary[mydata.summary$R == 0.03,], aes(init.ccount, init.placement, fill = consensus75))
p <- p + facet_grid(R~.)
p <- p + geom_tile(colour = "white") + scale_fill_manual()
p

wireframe(consensus75 ~ init.ccount *init.placement, data = mydata.summary[mydata.summary$R == 0.3,])

#p <- p + geom_line(aes(color = init.placement), alpha=0.5)
p <- p + xlab("Initial Number of Clusters") + ylab("Avg. Time to Consensus")
p <- p + geom_smooth(data = mydata.subset011b[mydata.subset011b$R == 0.3,],
                     size = 1.1, alpha = 0.2,
                     se=FALSE,
                     #linetype = "dashed",
                     aes(color = as.factor(init.placement)))
p <- p + facet_grid(alpha2~R, labeller=myLabeller)
p <- p + scale_color_discrete(name="Initial Distance\nfrom Truth")
#p <- p + scale_x_continuous(breaks=c(1,5,10,15))
p <- p + myThemeMod + theme(legend.position = c(0.15, 0.80),
                            legend.background = element_rect(fill=alpha('white', 0.4)),
                            #legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(size=18,face="bold"),
                            legend.text = element_text(size=18),
                            legend.key.height = unit(0.7, "cm"),
                            legend.key.width = unit(1, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA),
                            strip.background = element_blank()
                            )
p



# test
myBreaks = c(1.4892, 1.5892, 1.6892, 1.7892, 1.8892, 1.9892, 2.0892, 2.1892, 2.2892)
title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary[mydata.summary$R == 0.03 & mydata.summary$init.placement %in% myBreaks, ], aes(init.ccount, consensus75, group=init.placement))
#p <- p + geom_point(size=2, alpha=0.2)
#p <- p + geom_line(alpha=0.5)
p <- p + geom_smooth(aes(color=as.factor(init.placement)), se=FALSE, method="lm")
p <- p + xlab("Initial Number of Clusters") + ylab("Avg. Time to Consensus")
p <- p + geom_smooth(data = mydata.subset011b[mydata.subset011b$R == 0.3,],
                     size = 1.1, alpha = 0.2,
                     se=FALSE,
                     #linetype = "dashed",
                     aes(color = as.factor(init.placement)))
p <- p + facet_grid(alpha2~R, labeller=myLabeller)
p <- p + scale_color_discrete(name="Initial Distance\nfrom Truth")
#p <- p + scale_x_continuous(breaks=c(1,5,10,15))
p

p <- p + myThemeMod + theme(legend.position = c(0.15, 0.80),
                            legend.background = element_rect(fill=alpha('white', 0.4)),
                            #legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(size=18,face="bold"),
                            legend.text = element_text(size=18),
                            legend.key.height = unit(0.7, "cm"),
                            legend.key.width = unit(1, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA),
                            strip.background = element_blank()
                            )
p

# ORIGINAL

# [mydata.subset011$R == 0.3,]

mydata.subset <- mydata.summary[mydata.summary$init.ccount < 16,]
extremes <- c(min(mydata$init.placement),max(mydata$init.placement))
mydata.subset011 <- mydata.summary[mydata.summary$init.ccount < 16 & mydata.summary$init.placement %in% extremes,]


mydata.subset <- mydata.summary[mydata.summary$init.ccount < 16 & mydata.summary$init.placement > limitPlacement,]
extremes <- c(min(mydata.subset$init.placement),max(mydata.subset$init.placement))
mydata.subset011 <- mydata.summary[mydata.summary$init.ccount < 16 & mydata.summary$init.placement %in% extremes,]

# Fake small R
mydata.subset <- rbind(mydata.subset, c(0.03, 0.01, NA, NA, NA, NA, NA, NA, NA))

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.subset, aes(init.ccount, consensus75))
p <- p + geom_jitter(size=2)
#p <- p + geom_line(aes(color = init.placement), alpha=0.5)
p <- p + xlab("Initial Number of Clusters") + ylab("Avg. Time to Consensus")
p <- p + geom_smooth(data = mydata.subset011,
                     size = 1.1, alpha = 0.2,
                     se=FALSE, method="lm",
                     aes(color = as.factor(init.placement)))
p <- p + facet_grid(alpha2~R, labeller=myLabeller)
p <- p + scale_color_discrete(name="Initial Distance\nfrom Truth")
p <- p + scale_x_continuous(breaks=c(1,5,10,15))
p <- p + myThemeMod + theme(legend.position = c(0.15, 0.80),
                            legend.background = element_rect(fill=alpha('white', 0.4)),
                            #legend.background = element_rect(fill = "white", colour = "grey"),
                            legend.title = element_text(size=18,face="bold"),
                            legend.text = element_text(size=18),
                            legend.key.height = unit(0.7, "cm"),
                            legend.key.width = unit(1, "cm"),
                            legend.key =  element_rect(fill = "white", colour = NA),
                            strip.background = element_blank()
                            )
p

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[13]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.01"))
grob[["grobs"]][[14]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.5"))
grob[["grobs"]][[15]][["children"]][[2]][["label"]] <- expression(paste(alpha," = 0.99"))


name <- "SPEEDTEST/non_overlapping_clusters.svg"
name <- "SPEEDTEST/non_overlapping_clusters_biggap.svg"

svg(filename = paste0(IMGPATH, name), width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()

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



## NEW PLOTS (old)


mydata.summary <- summarySE(mydata, "ccount50", c("init.ccount", "init.placement"), na.rm = TRUE)

title <- "The effect of progress and clustering on consensus"
p <- ggplot(mydata.summary, aes(init.ccount, ccount50))
p <- p + geom_point(aes(color = as.factor(init.placement)), size=4)
p <- p + geom_line(aes(color = as.factor(init.placement)))
p <- p + xlab("Initial number of clusters") + ylab("Avg. Number of Clusters at Consensus Share 50%")
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

## NEW PLOTS

mydata.summary <- summarySE(mydata, "ccount50", c("alpha", "R", "init.ccount", "init.placement", "clbr"), na.rm = TRUE)

title <- "Number of clusters by distance from truth and initial number of clusters"
p <- ggplot(mydata.summary[mydata.summary$alpha != 0.01 & mydata.summary$R == 0.03 & mydata.summary$init.placement > limitPlacement & mydata.summary$clbr != "(0,1]",], aes(as.factor(init.placement), ccount50, color=as.factor(init.placement)))
p <- p + geom_boxplot()
p <- p + xlab("Initial Distance from Truth") + ylab("Avg. Number of Clusters at Consensus Share 50%")
p <- p + facet_grid(clbr~., labeller = myLabeller)
p <- p + scale_color_hue(name="Initial\ndistance\nfrom Truth")
p <- p + scale_x_discrete(breaks = c(1.4892, 1.5892, 1.6892, 1.7892,
                            1.8892, 1.9892, 2.0892, 2.1892, 2.2892),
                          labels = c(1.48,1.58,1.68,1.78,1.88,1.98,2.08,2.18,2.28))
p <- p + ylim(1,5.75)
#p <- p + scale_x_discrete(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
p <- p + myThemeMod + theme(
                            strip.background = element_blank()
                            )
p


ggsave(filename = paste0(IMGPATH, "SPEEDTEST/progress_leads_to_clustering_boxplot_biggap_c50.svg"),
       width=10, height=10, dpi=300)

na50 <- mydata[is.na(mydata$ccount50) & mydata$alpha != 0.99,]


## Regressions

fit001 <- lm(consensus75 ~ ccount30, data = mydata[mydata$alpha == 0.01 & mydata$R == 0.03,]); summary(fit001)
fit05 <- lm(consensus75 ~ ccount30, data = mydata[mydata$alpha == 0.5 & mydata$R == 0.03,]); summary(fit05)

fit099 <- lm(consensus75 ~ ccount30, data = mydata[mydata$alpha == 0.01 & mydata$R == 0.3,]); summary(fit099)


# high social influence
fit001b <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.01 & mydata$R == 0.03,]); summary(fit001b)

# medium social influence
fit05b <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.5 & mydata$R == 0.03,]); summary(fit05b)

#bogus
fit099b <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.99 & mydata$R == 0.03,]); summary(fit099b)


fit05bRbig <- lm(consensus75 ~ init.ccount, data = mydata[mydata$alpha == 0.5 & mydata$R == 0.3,]); summary(fit05bRbig)

# init.dist = 1
fit05bRbigInit1 <- lm(consensus75 ~ init.ccount, data = mydata[mydata$init.placement == 0.6 & mydata$alpha == 0.5 & mydata$R == 0.3,]); summary(fit05bRbigInit1)

# init.dist = 1
fit05bRbigInit1 <- lm(consensus75 ~ init.ccount, data = mydata[mydata$init.placement == 0.2 & mydata$alpha == 0.5 & mydata$R == 0.3,]); summary(fit05bRbigInit1)



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

# Old
# cl <- loadData(DUMPDIR, 'scan_tau_again2/', 1)
# cl <- loadData(DUMPDIR, 'final_tau_20000_nobound/', 1)


clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)
clTau1$boundaries <- 0

clAgain3 <- loadData(DUMPDIR, 'scan_tau_again3/', 1)

cl <- rbind(clTau1, clAgain3)

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("tau", "R"), na.rm=TRUE)

p <- ggplot(summaryCl, aes(1-tau/100, count))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.011))
p <- p + geom_errorbar(limits)
p <- p + facet_grid(.~ R, labeller = myLabeller)
xlabText <- expression(paste('Strength of Attraction to Ground Truth 1/',tau))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
p <- p + scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
labels = c("0.1", "0.25", "0.5", "0.75", "1"))
p <- p + myThemeMod + theme(strip.background = element_blank(),
legend.position = "none"
)
p

ggsave(filename = paste0(IMGPATH, "scan_tau_nobound_cc.svg"),
plot = p, width=10, height=5, dpi=300)

## FT

summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("tau", "R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(1-tau/100, fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(.~ R, labeller = myLabeller)
xlabText <- expression(paste('Strength of Attraction to Ground Truth 1/',tau))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
labels = c("0.1", "0.25", "0.5", "0.75", "1"))
p <- p + myThemeMod + theme(strip.background = element_blank(),
legend.position = "none")
p

ggsave(filename = paste0(IMGPATH, "scan_tau_nobound_ft.svg"),
       plot = p, width=10, height=5, dpi=300)


## R 2 ##

cl2 <- loadData(DUMPDIR, 'scan_R_again2/', 1)
cl3 <- loadData(DUMPDIR, 'scan_R_again3/', 1)

clTau1 <- loadData(DUMPDIR, 'nobound_R_tau1/', 1)
clTau1$boundaries <- 0

cl <- clTau1

cl <- rbind(cl2, cl3, clTau1[clTau1$alpha == 0.5,])

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

## alpha 2 ##

cl <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1)
clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)


cl <- rbind(cl, clTau1)
# Only tau 1
cl <- clTau1


clTauAgain <- loadData(DUMPDIR, 'scan_alpha_again3/', 1)

clTau1$boundaries <- 0
cl <-  rbind(clTauAgain, clTau1)

cl <- clTauAgain

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


## NOISES ##
############

clTauAgain <- loadData(DUMPDIR, 'scan_noises_again2/', 1)
  
clTau1 <- loadData(DUMPDIR, 'nobound_noises_tau1/', 1)

clTau1$boundaries <- 0
cl <- rbind(clTauAgain, clTau1[clTau1$alpha == 0.5,])

                        
# cl <- clTauAgain
                

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




### Mediation and Moderation Analysis

DIR <- 'speedtest_final_R_alpha/'

DIR <- 'nobound_mediation_R_alpha/'
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
data[ , 15:27 ][ data[ , 15:27 ] == -1 ] <- 20000
# If no consensus was reached replace ccounts with NA
data[ , 28:40 ][ data[ , 28:40 ] == -1 ] <- NA
# Replaces everything, not good
#data[] <- lapply(data, function(x){replace(x, x == -1, 20000)})
data$R <- as.factor(data$R)
data$alpha <- as.factor(data$alpha)


mydata <- data[complete.cases(data[,c("consensus75","ccount50")]) & data$alpha == 0.5,]


# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(consensus75 ~ smallR, data = mydata)
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(ccount50 ~ smallR, data = mydata)
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR: R and alpha on consensus75

fit1c <- lm(consensus75 ~ smallR + ccount50, data = mydata)
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="smallR", mediator="ccount50")

summary(contcont)


# Test whether progress causes the mediator

fit1c <- lm(ccount50 ~ consensus75 + smallR, data = mydata)
summary(fit1c)


# All Experiment 1 data


cl <- loadData(DUMPDIR, 'nobound_R_tau/', 1)
clTau1 <- loadData(DUMPDIR, 'nobound_R_tau1/', 1)
cl <- rbind(cl, clTau1[clTau1$alpha == 0.5,])
clOld <- cl
cl <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1)
clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)
cl <- rbind(cl, clTau1)
# Only tau 1
cl <- clTau1
clOld <- rbind(cl, clOld)
cl <- loadData(DUMPDIR, 'nobound_noises_tau/', 1)  
clTau1 <- loadData(DUMPDIR, 'nobound_noises_tau1/', 1)
cl <- rbind(cl, clTau1[clTau1$alpha == 0.5,])
# Only Tau 1
cl <- clTau1
clOld <- rbind(cl, clOld)
cl <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1)
clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 1)
cl <- rbind(cl, clTau1)
clOld <- rbind(cl, clOld)
data <- clOld

data$convZone <- ifelse(data$R < 0.11,0,1)
data$taubrk <- cut(data$tau,c(0,1,20,50,100))
data$alpha2 <- 1 - data$alpha

p <- ggplot(data[data$t == 2000,], aes(count, fromtruth.avg, color=alpha2))
p <- p + geom_jitter(alpha=0.5)
p <- p + facet_grid(taubrk~convZone, labeller=myLabeller3)
p <- p + xlab('Avg. Number of Clusters') + ylab('Avg. Distance from Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.position = "right")
p

p <- ggplot(data[data$t == 2000,], aes(count, fromtruth.avg))
p <- p + geom_jitter(alpha=0.2)
p <- p + geom_smooth(method="lm", se = FALSE)
p <- p + xlab('Avg. Number of Clusters') + ylab('Avg. Distance from Truth')
p <- p + ylim(0,max(data$fromtruth.avg)+0.5)
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.position = "right")
p



ggsave(filename = paste0(IMGPATH, "nobound_clusters_vs_progress.jpg"),
       plot = p, dpi=300)


xlabText <- expression(paste("Estimated passage of ",italic("Fish fish"),"in 2001"))



# Mediation 2: R


clTau1 <- loadData(DUMPDIR, 'nobound_R_tau1/', 0)
data <- clTau1[clTau1$alpha == 0.5,]
data$id <- paste(data$simcount, data$run, sep='.')
data <- data[data$t %in% c(1000,2000),]
data <- subset(data, select=c("id","R", "t", "count", "fromtruth.avg"))
data1000 <- data[data$t == 1000,]
colnames(data1000) <- c("id1000","R1000","t","count1000","progress1000")
data1000$t <- NULL
data2000 <- data[data$t == 2000,]
colnames(data2000) <- c("id2000","R2000","t","count2000","progress2000")
data2000$t <- NULL
data2 <- cbind(data1000,data2000)
data2$ok <- ifelse(data2$id1000 != data2$id2000, 0,1)
data2 <- subset(data2, select=c("id2000","R2000","count1000","progress1000","count2000","progress2000"))
colnames(data2) <- c("id","R","count1000","progress1000","count2000","progress2000")
data2$smallR <- ifelse(data2$R <= 0.1, 1,0)

# Small R
# EFFECT 0.724
# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ R, data = data2[data2$smallR == 1,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ R, data = data2[data2$smallR == 1,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR

fit1c <- lm(progress2000 ~ R + count1000, data = data2[data2$smallR == 1,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="R", mediator="count1000")

summary(contcont)

# Big R
# Effect 0.71, but much smaller than Small R

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ R, data = data2[data2$smallR == 0,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ R, data = data2[data2$smallR == 0,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR: 

fit1c <- lm(progress2000 ~ R + count1000, data = data2[data2$smallR == 0,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="R", mediator="count1000")
summary(contcont)

# All R
# Effect 0.7
# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ R, data = data2)
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ R, data = data2[data2$smallR == 0,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR:

fit1c <- lm(progress2000 ~ R + count1000, data = data2[data2$smallR == 0,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="R", mediator="count1000")
summary(contcont)


# Mediation 2: alpha

data <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 0)
data$id <- paste(data$simcount, data$run, sep='.')
data <- data[data$t %in% c(1000,2000),]
data <- subset(data, select=c("id","alpha", "R", "t", "count", "fromtruth.avg"))
data1000 <- data[data$t == 1000,]
colnames(data1000) <- c("id1000","alpha1000","R1000","t","count1000","progress1000")
data1000$t <- NULL
data2000 <- data[data$t == 2000,]
colnames(data2000) <- c("id2000","alpha2000","R2000","t","count2000","progress2000")
data2000$t <- NULL
data2 <- cbind(data1000,data2000)
data2$ok <- ifelse(data2$id1000 != data2$id2000, 0,1)
mean(data2$ok)
data2 <- subset(data2, select=c("id2000","alpha2000","R2000","count1000","progress1000","count2000","progress2000"))
colnames(data2) <- c("id","alpha","R","count1000","progress1000","count2000","progress2000")
data2$smallR <- ifelse(data2$R <= 0.1, 1,0)

# Small R
# No mediation, the coefficient is not reduced in size.

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ alpha, data = data2[data2$smallR == 1,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ alpha, data = data2[data2$smallR == 1,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR:

fit1c <- lm(progress2000 ~ alpha + count1000, data = data2[data2$smallR == 1,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="alpha", mediator="count1000")
summary(contcont)

# Big R
# No mediation, no effect

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ alpha, data = data2[data2$smallR == 0,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ alpha, data = data2[data2$smallR == 0,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR:

fit1c <- lm(progress2000 ~ alpha + count1000, data = data2[data2$smallR == 0,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="alpha", mediator="count1000")
summary(contcont)

# All R
# No effect.
# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ alpha, data = data2)
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ alpha, data = data2)
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR: R and alpha on consensus75

fit1c <- lm(progress2000 ~ alpha + count1000, data = data2)
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="alpha", mediator="count1000")
summary(contcont)


# Mediation 2: noises

data <- loadData(DUMPDIR, 'nobound_noises_tau1/', 0)
data$id <- paste(data$simcount, data$run, sep='.')
data <- data[data$t %in% c(1000,2000),]
data <- subset(data, select=c("id","epsilon", "sigma", "R", "t", "count", "fromtruth.avg"))
data1000 <- data[data$t == 1000,]
colnames(data1000) <- c("id1000","epsilon1000","sigma1000","R1000","t","count1000","progress1000")
data1000$t <- NULL
data2000 <- data[data$t == 2000,]
colnames(data2000) <- c("id2000","epsilon2000","sigma2000","R2000","t","count2000","progress2000")
data2000$t <- NULL
data2 <- cbind(data1000,data2000)
data2$ok <- ifelse(data2$id1000 != data2$id2000, 0,1)
mean(data2$ok)
data2 <- subset(data2, select=c("id2000","epsilon2000","sigma2000", "R2000","count1000","progress1000","count2000","progress2000"))
colnames(data2) <- c("id","epsilon","sigma","R","count1000","progress1000","count2000","progress2000")
data2$smallR <- ifelse(data2$R <= 0.1, 1,0)


# EPSILON

# Small R

# Step1 Regress INDEPENDENT on DEPENDENT
# Not significant.
fit1a <- lm(progress2000 ~ epsilon, data = data2[data2$smallR == 1,])
summary(fit1a)

# Big R
# Minimal effect 0.09

# Step1 Regress INDEPENDENT on DEPENDENT
# Not significant.
fit1a <- lm(progress2000 ~ epsilon, data = data2[data2$smallR == 0,])
summary(fit1a)

# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ epsilon, data = data2[data2$smallR == 0,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR:

fit1c <- lm(progress2000 ~ epsilon + count1000, data = data2[data2$smallR == 0,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="epsilon", mediator="count1000")
summary(contcont)

# Big R

# Step1 Regress INDEPENDENT on DEPENDENT
# Not significant.
fit1a <- lm(progress2000 ~ epsilon, data = data2)
summary(fit1a)

# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ epsilon, data = data2)
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR:

fit1c <- lm(progress2000 ~ epsilon + count1000, data = data2)
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="epsilon", mediator="count1000")
summary(contcont)

# Small R

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ sigma, data = data2[data2$smallR == 1,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ sigma, data = data2[data2$smallR == 1,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR: 

fit1c <- lm(progress2000 ~ sigma + count1000, data = data2[data2$smallR == 1,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="sigma", mediator="count1000")
summary(contcont)

# Big R
# No Effect

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ sigma, data = data2[data2$smallR == 0,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ sigma, data = data2[data2$smallR == 0,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR: 

fit1c <- lm(progress2000 ~ sigma + count1000, data = data2[data2$smallR == 0,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="sigma", mediator="count1000")
summary(contcont)

# All R
# No Effect

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ sigma, data = data2)
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ sigma, data = data2)
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR: 

fit1c <- lm(progress2000 ~ sigma + count1000, data = data2)
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="sigma", mediator="count1000")
summary(contcont)

# Mediation 2: tau


# Need to load it in batches, otherwise I get memory error.
cl2000 <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1, 2000)
cl1000 <- loadData(DUMPDIR, 'nobound_alpha_tau/', 1, 1000)
cl <- rbind(cl1000[cl2000$alpha == 0.5,], cl2000[cl2000$alpha == 0.5,])
clTau1 <- loadData(DUMPDIR, 'nobound_alpha_tau1/', 0)
data <- rbind(cl, clTau1[clTau1$alpha == 0.5 & clTau1$t %in% c(1000,2000),])
data$id <- paste(data$simcount, data$run, sep='.')
data <- data[data$t %in% c(1000,2000),]
data <- subset(data, select=c("id","tau", "R", "t", "count", "fromtruth.avg"))
data1000 <- data[data$t == 1000,]
colnames(data1000) <- c("id1000","tau1000","R1000","t","count1000","progress1000")
data1000$t <- NULL
data2000 <- data[data$t == 2000,]
colnames(data2000) <- c("id2000","tau2000","R2000","t","count2000","progress2000")
data2000$t <- NULL
data2 <- cbind(data1000,data2000)
data2$ok <- ifelse(data2$id1000 != data2$id2000, 0,1)
mean(data2$ok)
data2 <- subset(data2, select=c("id2000","tau2000","R2000","count1000","progress1000","count2000","progress2000"))
colnames(data2) <- c("id","tau","R","count1000","progress1000","count2000","progress2000")
data2$smallR <- ifelse(data2$R <= 0.1, 1,0)

# SMALL R
# Little mediation 0.39

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ tau, data = data2[data2$smallR == 1,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ tau, data = data2[data2$smallR == 1,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR

fit1c <- lm(progress2000 ~ tau + count1000, data = data2[data2$smallR == 1,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="tau", mediator="count1000")
summary(contcont)

# LARGE R
# High mediation 0.77, but the effect is really small

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ tau, data = data2[data2$smallR == 0,])
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ tau, data = data2[data2$smallR == 0,])
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR

fit1c <- lm(progress2000 ~ tau + count1000, data = data2[data2$smallR == 0,])
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="tau", mediator="count1000")
summary(contcont)


# ANY R
# The effect even changes sign

# Step1 Regress INDEPENDENT on DEPENDENT

fit1a <- lm(progress2000 ~ tau, data = data2)
summary(fit1a)


# Step2 Regress MEDIATOR on INDEPENDENT 

fit1b <- lm(count1000 ~ tau, data = data2)
summary(fit1b)


# Step3 Regress DEPENDENT on INDEPENDENT + MEDIATOR

fit1c <- lm(progress2000 ~ tau + count1000, data = data2)
summary(fit1c)


# Estimation via quasi-Bayesian approximation
contcont <- mediate(fit1b, fit1c, sims=500, treat="tau", mediator="count1000")
summary(contcont)
