# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")
source("/opt/MATLAB_WORKSPACE/hs/R/init.NOBOUND.R")
library(scales)

# DUMPDIR 
DUMPDIR <- '/home/stefano/Documents/mypapers/swarm_science/data/'
## IMG DIR
IMGPATH <- paste0(DUMPDIR, "imgs/NOBOUND/aft/")
# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(IMGPATH))
}

limitsFtA <- aes(ymin = a.fromtruth.avg - se, ymax = a.fromtruth.avg + se)

## TAU 2 ##

# Weird Values
# cl <- loadData(DUMPDIR, 'scan_tau_again2/', 1)

# cl <- loadData(DUMPDIR, 'final_tau_20000_nobound/', 1)

clTau1 <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau1/', 1)
clTau1$boundaries <- 0

clAgain3 <- loadDataAgents(DUMPDIR, 'scan_tau_again3/', 1)

cl <- rbind(clTau1[clTau1$alpha == 0.5,], clAgain3)

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

summaryFt <- summarySE(cl[cl$t == 2000,], c("a.fromtruth.avg"), c("tau", "R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(1-tau/100, a.fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=a.fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFtA)
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

cl2 <- loadDataAgents(DUMPDIR, 'scan_R_again2/', 1)
cl3 <- loadDataAgents(DUMPDIR, 'scan_R_again3/', 1)
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_R_tau1/', 1)
clTau1$boundaries <- 0
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
summaryFt <- summarySE(cl[cl$t == 2000,], c("a.fromtruth.avg"), c("R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(R, a.fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=a.fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFtA)
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

#cl <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau/', 1)

clTau1 <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau1/', 1)
clTau1$boundaries <- 0
# Only tau 1
#cl <- clTau1

clTauAgain <- loadDataAgents(DUMPDIR, 'scan_alpha_again3/', 1)
cl <-  rbind(clTauAgain, clTau1)

# cl <- clTauAgain

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
summaryFt <- summarySE(cl[cl$t == 2000,], c("a.fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes((1 - alpha), a.fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=a.fromtruth.avg, width=0.01))
p <- p + geom_errorbar(limitsFtA)
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

clTauAgain <- loadDataAgents(DUMPDIR, 'scan_noises_again2/', 1)
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_noises_tau1/', 1)
clTau1$boundaries <- 0
cl <- rbind(clTauAgain, clTau1[clTau1$alpha == 0.5,])
                

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
summaryFt <- summarySE(cl[cl$t == 2000,], c("a.fromtruth.avg"), c("sigma", "epsilon", "R"), na.rm=TRUE)

p <- ggplot(summaryFt, aes(sigma, a.fromtruth.avg))
p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=a.fromtruth.avg))
p <- p + geom_errorbar(limitsFtA)
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


# All Experiment 1 data (What is this???? Seems incorrect)


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

p <- ggplot(data[data$t == 2000,], aes(count, a.fromtruth.avg, color=alpha2))
p <- p + geom_jitter(alpha=0.5)
p <- p + facet_grid(taubrk~convZone, labeller=myLabeller3)
p <- p + xlab('Avg. Number of Clusters') + ylab('Avg. Distance from Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.position = "right")
p

p <- ggplot(data[data$t == 2000,], aes(count, a.fromtruth.avg))
p <- p + geom_jitter(alpha=0.2)
p <- p + geom_smooth(method="lm", se = FALSE)
p <- p + xlab('Avg. Number of Clusters') + ylab('Avg. Distance from Truth')
p <- p + ylim(0,max(data$a.fromtruth.avg)+0.5)
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
data <- subset(data, select=c("id","R", "t", "count", "a.fromtruth.avg"))
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
data <- subset(data, select=c("id","alpha", "R", "t", "count", "a.fromtruth.avg"))
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
data <- subset(data, select=c("id","epsilon", "sigma", "R", "t", "count", "a.fromtruth.avg"))
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
data <- subset(data, select=c("id","tau", "R", "t", "count", "a.fromtruth.avg"))
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



##########

# INTERACTION WITH TAU

## ALPHA ##
###########

cl <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau/', 1)
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau1/', 1)
cl <- rbind(cl, clTau1)
cl$alphabrk <- cut(cl$alpha, seq(0,1,0.05))
cl$taubrk <- cut(cl$tau, c(0,1,5,20,100))

myalphas <- c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99)

# CL

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("alpha", "R", "taubrk"), na.rm=TRUE)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes((1 - alpha), count, color=taubrk, group = taubrk))
#p <- p + geom_bar(stat = "identity", position="dodge", aes(group=taubrk, fill=count, width=0.01))
p <- p + geom_point()
p <- p + geom_line(alpha=0.5)
p <- p + geom_errorbar(limits)
p <- p + facet_grid(. ~ R, labeller = myLabeller)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of Social Influence ',alpha))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
p <- p + scale_color_hue(labels=c("1", "2-5", "6-20", "21-100"), name="Tau Levels")
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = c(0.5,0.5))
p

ggsave(filename = paste0(IMGPATH, "nobound_interaction_tau_alpha_cc.svg"),
       plot = p, width=10, height=5, dpi=300)


# FT

summaryFt <- summarySE(cl[cl$t == 2000 &
                          cl$alpha %in% myalphas &
                          cl$tau %in% c(1,10,50,100)
                          ,], c("fromtruth.avg"), c("alpha", "R", "tau"), na.rm=TRUE)

summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("alpha", "R", "taubrk"), na.rm=TRUE)

p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg, color=taubrk, group = taubrk))
#p <- p + geom_bar(stat = "identity", position="dodge", aes(group=taubrk, fill=count, width=0.01))
p <- p + geom_point()
p <- p + geom_line(alpha=0.5)
p <- p + geom_errorbar(limitsFt)
p <- p + facet_grid(. ~ R, labeller = myLabeller)
p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
xlabText <- expression(paste('Strength of Social Influence ',alpha))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_color_hue(labels=c("1", "2-5", "6-20", "21-100"), name="Tau Levels")
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = c(0.5,0.5)
                            )
p

ggsave(filename = paste0(IMGPATH, "nobound_interaction_tau_alpha_ft.svg"),
       plot = p, width=10, height=5, dpi=300)



## NOISES ##
############

cl <- loadDataAgents(DUMPDIR, 'nobound_noises_tau/', 1)  
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_noises_tau1/', 1)
cl <- rbind(cl, clTau1[clTau1$alpha == 0.5,])
cl$taubrk <- cut(cl$tau, c(0,1,5,20,100))

summaryCl <- summarySE(cl[cl$t == 2000,], c("count"), c("sigma", "epsilon", "R", "taubrk"), na.rm=TRUE)


# Only for epsilon = 1
title <- 'Cluster counts vs Angular noise and \nPosition noise'
p <- ggplot(summaryCl[summaryCl$epsilon == 0.1,], aes(sigma, count, group = taubrk, color = taubrk))
#p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_point(size=3)
p <- p + geom_line(alpha=0.5)
p <- p + geom_errorbar(limits, width=0.005)
p <- p + facet_grid(. ~ R, labeller = myLabeller)
xlabText <- expression(paste('Angular Noise ',sigma))
p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
p <- p + scale_color_hue(labels=c("1", "2-5", "6-20", "21-100"), name="Tau Levels")
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.text = element_text(size=18),
                            legend.position = c(0.7,0.8)
                            )
#p <- p + guides(col=guide_legend(ncol=4))
p


ggsave(filename = paste0(IMGPATH, "nobound_interaction_tau_alpha_cc_epsilon01.svg"),
       plot = p, width=10, height=5, dpi=300)

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[22]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0"))
grob[["grobs"]][[23]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.1"))
grob[["grobs"]][[24]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.2"))
grob[["grobs"]][[25]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.3"))
grob[["grobs"]][[26]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.4"))
grob[["grobs"]][[27]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.5"))

svg(filename = paste0(IMGPATH, "nobound_interaction_tau_noises_cc.svg"),
     width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()


# FT
summaryFt <- summarySE(cl[cl$t == 2000,], c("fromtruth.avg"), c("sigma", "epsilon", "R", "taubrk"), na.rm=TRUE)


p <- ggplot(summaryFt[summaryFt$epsilon == 0.1,], aes(sigma, fromtruth.avg, group = taubrk, color = taubrk))
#p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count))
p <- p + geom_point(size=3)
p <- p + geom_line(alpha=0.5)
p <- p + geom_errorbar(limitsFt, width=0.005)
p <- p + facet_grid(. ~ R, labeller = myLabeller)
xlabText <- expression(paste('Angular Noise ',sigma))
p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
p <- p + scale_x_continuous(labels = c("0", "0.025", "0.05","0.075", "0.1"))
p <- p + scale_y_continuous(breaks=c(0,2,4), labels = c("0", "2", "4"))
p <- p + scale_color_hue(labels=c("1", "2-5", "6-20", "21-100"), name="Tau Levels")
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = c(0.5,0.5)
                            )
#p <- p + guides(col=guide_legend(ncol=4))
p



ggsave(filename = paste0(IMGPATH, "nobound_interaction_tau_noises_ft_epsilon01.svg"),
       plot = p, width=10, height=5, dpi=300)

# Unfortunately, have to use this weird way of setting the labels, because facet labeller
# has a problem with the expression method.
grob <- ggplotGrob(p)

grob[["grobs"]][[22]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0"))
grob[["grobs"]][[23]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.1"))
grob[["grobs"]][[24]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.2"))
grob[["grobs"]][[25]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.3"))
grob[["grobs"]][[26]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.4"))
grob[["grobs"]][[27]][["children"]][[2]][["label"]] <- expression(paste(epsilon," = 0.5"))

svg(filename = paste0(IMGPATH, "nobound_interaction_tau_noises_ft.svg"),
     width=10, height=10)
grid.newpage()
grid.draw(grob)
dev.off()




###

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



# All Experiment 1 data

# Tau
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau1/', 1)
clTau1$boundaries <- 0
clAgain3 <- loadDataAgents(DUMPDIR, 'scan_tau_again3/', 1)
#cl <- rbind(clTau1[clTau1$alpha == 0.5,], clAgain3)
cl <- rbind(clTau1, clAgain3)
clOld <- cl
# R
cl2 <- loadDataAgents(DUMPDIR, 'scan_R_again2/', 1)
cl3 <- loadDataAgents(DUMPDIR, 'scan_R_again3/', 1)
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_R_tau1/', 1)
clTau1$boundaries <- 0
#cl <- rbind(cl2, cl3, clTau1[clTau1$alpha == 0.5,])
cl <- rbind(cl2, cl3, clTau1)
clOld <- rbind(cl, clOld)
# Alpha
clTauAgain <- loadDataAgents(DUMPDIR, 'scan_alpha_again3/', 1)
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_alpha_tau1/', 1)
clTau1$boundaries <- 0
cl <- rbind(clTauAgain, clTau1)
clOld <- rbind(cl, clOld)
# Noises
clTauAgain <- loadDataAgents(DUMPDIR, 'scan_noises_again2/', 1)
clTau1 <- loadDataAgents(DUMPDIR, 'nobound_noises_tau1/', 1)
clTau1$boundaries <- 0
#cl <- rbind(clTauAgain, clTau1[clTau1$alpha == 0.5,])
cl <- rbind(clTauAgain)
clOld <- rbind(cl, clOld)
# All Experiment 1 data
data <- clOld

data$convZone <- ifelse(data$R < 0.11,0,1)
data$taubrk <- cut(data$tau,c(0,1,20,50,100))
data$alpha2 <- 1 - data$alpha

p <- ggplot(data[data$t == 2000,], aes(count, a.fromtruth.avg, color=alpha2))
p <- p + geom_jitter(alpha=0.5)
p <- p + facet_grid(taubrk~convZone, labeller=myLabeller3)
p <- p + xlab('Avg. Number of Clusters') + ylab('Avg. Distance from Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.position = "right")
p

p <- ggplot(data[data$t == 2000,], aes(count, a.fromtruth.avg))
p <- p + geom_jitter(alpha=0.8)
p <- p + geom_smooth(method="lm", se = FALSE, size = 2)
p <- p + xlab('Avg. Number of Clusters') + ylab('Avg. Distance from Truth')
p <- p + ylim(0,max(data$fromtruth.avg)+0.5)
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.position = "right")
p



ggsave(filename = paste0(IMGPATH, "nobound_clusters_vs_progress_aft.jpg"),
       plot = p, dpi=300)



#### 20000


cl <- loadDataAgents(DUMPDIR, 'nobound_alpha_20000/')




summaryFt <- summarySE(cl[cl$t == 20000,], c("a.fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)


p



# Save all Taus Cl
taus <- seq(1000,20000,1000)
for (t in taus) {
  summaryFt <- summarySE(cl[cl$t == t,], c("a.fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("T: ", t)
  p <- ggplot(summaryFt, aes((1 - alpha), a.fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=a.fromtruth.avg, width=0.01))
  p <- p + geom_errorbar(limitsFtA)
  p <- p + facet_grid(~ R, labeller = myLabeller)
  xlabText <- expression(paste('Strength of Social Influence ', alpha))
  p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
  p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  p <- p + scale_fill_continuous(name="Distance\nfrom truth")
  p <- p + myThemeMod + theme(strip.background = element_blank())
  p <- p + ggtitle(title)
  #
  ggsave(filename=paste0(IMGPATH, "alpha20000/ft/alpha_t_", sprintf("%04d", t), ".jpg"),
         plot = p)
}

# Save all Taus Cl
taus <- seq(1000,20000,1000)
for (t in taus) {
  summaryFt <- summarySE(cl[cl$t == t,], c("fromtruth.avg"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("T: ", t)
  p <- ggplot(summaryFt, aes((1 - alpha), fromtruth.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=fromtruth.avg, width=0.01))
  p <- p + geom_errorbar(limitsFt)
  p <- p + facet_grid(~ R, labeller = myLabeller)
  xlabText <- expression(paste('Strength of Social Influence ', alpha))
  p <- p + xlab(xlabText) + ylab('Avg. Distance from Truth')
  p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  p <- p + scale_fill_continuous(name="Distance\nfrom truth")
  p <- p + myThemeMod + theme(strip.background = element_blank())
  p <- p + ggtitle(title)
  #
  ggsave(filename=paste0(IMGPATH, "alpha20000/ft_clusters/alpha_t_", sprintf("%04d", t), ".jpg"),
         plot = p)
}

# Save all Taus Cl
taus <- seq(1000,20000,1000)
for (t in taus) {
  summaryCl <- summarySE(cl[cl$t == t,], c("count"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("T: ", t)
  p <- ggplot(summaryCl, aes((1 - alpha), count))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=count, width=0.01))
  p <- p + geom_errorbar(limits)
  p <- p + facet_grid(~ R, labeller = myLabeller)
  p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  xlabText <- expression(paste('Strength of Social Influence ',alpha))
  p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
  p <- p + myThemeMod + theme(strip.background = element_blank())
  p <- p + ggtitle(title)
  #
  ggsave(filename=paste0(IMGPATH, "alpha20000/cl/alpha_t_", sprintf("%04d", t), ".jpg"),
         plot = p)
}

# Save all Taus Cl
taus <- seq(1000,20000,1000)
for (t in taus) {
  summaryCl <- summarySE(cl[cl$t == t,], c("size.avg"), c("alpha", "R"), na.rm=TRUE)
  #
  title <- paste0("T: ", t)
  p <- ggplot(summaryCl, aes((1 - alpha), size.avg))
  p <- p + geom_bar(stat = "identity", position="dodge", aes(fill=size.avg, width=0.01))
  p <- p + geom_errorbar(aes(ymin = size.avg - se, ymax = size.avg + se))
  p <- p + facet_grid(~ R, labeller = myLabeller)
  p <- p + scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  xlabText <- expression(paste('Strength of Social Influence ',alpha))
  p <- p + xlab(xlabText) + ylab('Avg. Number of Clusters')
  p <- p + myThemeMod + theme(strip.background = element_blank())
  p <- p + ggtitle(title)
  #
  ggsave(filename=paste0(IMGPATH, "alpha20000/sizeavg/alpha_t_", sprintf("%04d", t), ".jpg"),
         plot = p)
}



#### 20000 - 500


cl <- loadDataAgents(DUMPDIR, 'nobound_alpha_20000_500/')



summaryCl <- summarySE(cl[cl$alpha %in% c(0.01, 0.5, 0.99),], c("count"), c("alpha", "t"), na.rm=TRUE)
summaryCl$alpha <- as.factor(summaryCl$alpha)

title <- 'Cluster counts vs Strength of social influence'
p <- ggplot(summaryCl, aes(t, count, color=alpha, group = alpha))
#p <- p + geom_bar(stat = "identity", position="dodge", aes(group=taubrk, fill=count, width=0.01))
p <- p + geom_point()
p <- p + geom_line(alpha=0.5)
p <- p + geom_errorbar(limits)
p <- p + xlab("Time") + ylab('Avg. Number of Clusters')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = "right")
p



summaryFt <- summarySE(cl[cl$alpha %in% c(0.01, 0.5, 0.99),], c("a.fromtruth.avg"), c("alpha", "t"), na.rm=TRUE)
summaryFt$alpha <- as.factor(summaryFt$alpha)

title <- 'Avg Dist. From Truth vs Strength of social influence'
p <- ggplot(summaryFt, aes(t, a.fromtruth.avg, color=alpha, group = alpha))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(limitsFtA)
p <- p + xlab("Time") + ylab('Avg. Distance From Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = "right")
p

summarySizeAvg <- summarySE(cl[cl$alpha %in% c(0.01, 0.5, 0.99),], c("size.avg"), c("alpha", "t"), na.rm=TRUE)
summarySizeAvg$alpha <- as.factor(summarySizeAvg$alpha)

title <- 'Avg Dist. From Truth vs Strength of social influence'
p <- ggplot(summarySizeAvg, aes(t, size.avg, color=alpha, group = alpha))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(aes(ymin=size.avg - se, ymax = size.avg + se))
p <- p + xlab("Time") + ylab('Avg. Distance From Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = "right")
p

summarySizeSd <- summarySE(cl[cl$alpha %in% c(0.01, 0.5, 0.99),], c("size.sd"), c("alpha", "t"), na.rm=TRUE)
summarySizeSd$alpha <- as.factor(summarySizeSd$alpha)

title <- 'Avg Dist. From Truth vs Strength of social influence'
p <- ggplot(summarySizeSd, aes(t, size.sd, color=alpha, group = alpha))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(aes(ymin=size.sd - se, ymax = size.sd + se))
p <- p + xlab("Time") + ylab('Avg. Distance From Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = "right")
p

summaryPDist.avg <- summarySE(cl[cl$alpha %in% c(0.01, 0.5, 0.99),], c("a.pdist.mean"), c("alpha", "t"), na.rm=TRUE)
summaryPDist.avg$alpha <- as.factor(summaryPDist.avg$alpha)

title <- 'Avg Dist. From Truth vs Strength of social influence'
p <- ggplot(summaryPDist.avg, aes(t, a.pdist.mean, color=alpha, group = alpha))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(aes(ymin=a.pdist.mean - se, ymax = a.pdist.mean + se))
p <- p + xlab("Time") + ylab('Avg. Distance From Truth')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = "right")
p

summaryPDist.sd <- summarySE(cl[cl$alpha %in% c(0.01, 0.5, 0.99),], c("a.pdist.sd"), c("alpha", "t"), na.rm=TRUE)
summaryPDist.sd$alpha <- as.factor(summaryPDist.sd$alpha)


p <- ggplot(summaryPDist.sd, aes(t, a.pdist.sd, color=alpha, group = alpha))
p <- p + geom_point()
p <- p + geom_line()
p <- p + geom_errorbar(aes(ymin=a.pdist.sd - se, ymax = a.pdist.sd + se))
p <- p + xlab("Time") + ylab('Avg. Std. Pairwise Distance')
p <- p + myThemeMod + theme(strip.background = element_blank(),
                            legend.title = element_text(vjust=3,
                              size=18, face="bold"),
                            legend.position = "right")
p
