# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

# DUMPDIR NEW
DUMPDIR = "/mnt/tmp/dump/NEW/"

# NO TRUTH
#DIR = "attrZero_nv_Kseed_rndseq_tm_Rleft/"

# TM
DIR = "attrLinear_nv_rndseq_tm_Rleft/"
#DIR = "attrK_nv_Kseed_rndseq_tm_Rleft/"
#DIR = "attrHard_nv_rndseed_rndseq_tm_Rleft/"
#DIR = "attrMillean_nv_rndseq_tm_Rleft/"
#DIR = "attrFunnel_nv_Kseed_rndseq_tm_Rleft/"
#DIR = "attrGentle_nv_Kseed_rndseq_tm_Rleft/" # some simulations are not saved in the results
#DIR = "attrExpo_nv_Kseed_rndseq_tm_Rleft/"

# TC
#DIR = "attrExpo_nv_seedFixed_rndseq_tc_Rleft/"
#DIR = "attrHard_nv_Kseed_rndseq_tc_Rleft/"
#DIR = "attrK_nv_Kseed_rndseq_tc_Rleft/"
#DIR = "attrMillean_nv_Kseed_rndseq_tc_Rleft/"
#DIR = "attrFunnel_nv_Kseed_rndseq_tc_Rleft/"
#DIR = "attrGentle_nv_Kseed_rndseq_tc_Rleft/" 
#DIR = "attrLinear_nv_Kseed_rndseq_tc_Rleft/"

# next
#DIR = "attrLinear_nv_rndseq_tm_Rleft/"

# DUMPDIR VELOCITY
#DUMPDIR = "/mnt/tmp/dump/VELOCITY/"

#DIR = "attrLinear_nv_Kseed_rndseq_tm_Rleft_n10/"
#DIR = "attrLinear_nv_Kseed_rndseq_tm_Rleft_n50/"
#DIR = "attrLinear_nv_Kseed_rndseq_tm_Rleft_n150/"
#DIR = "attrLinear_nv_Kseed_rndseq_tm_Rleft_n200/"

DIR = "attrLinear_nv_rndseq_tm_Rleft/attrLinear_nv_rndseq_tm_Rleft_s5/"

DUMPDIR = "/mnt/tmp/dump/RND_SEED/"
DIR = "attrZero_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0/attrZero_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0_s1"

INTERACTIVE = FALSE
PATH = paste0(DUMPDIR,DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");

# Create IMG dir if not existing
if (!file.exists(IMGPATH)) {
  dir.create(file.path(PATH, "/img/"))
} 

######################################
## MACRO AVG_SPLIT Aggregated: Loading
######################################

params <- read.table('params.csv', head=TRUE, sep=",")
#params$simname <- as.factor(params$simname)
#params$simcount <- as.factor(params$simcount)
#params$run <- as.factor(params$run)

macro <- read.table('clusters_macro.csv', head=TRUE, sep=",")
macro$simname <- as.character(macro$simname)
#macro$simname <- substr(macro$simname, nchar(macro$simname)-1, nchar(macro$simname))
#macro$simname <- as.factor(macro$simname)
#macro$simcount <- as.factor(macro$simcount)
#macro$t <- as.factor(macro$t)

cl <- macro
#cl <- subset(macro, t %% 100 == 0)
#cl$t <- as.factor(cl$t)

clOrig <- cl
cl$SIM <- as.numeric(cl$simcount)


INTERACTIVE = FALSE


## Testing weirdness in movements
A <- macro[macro$move.avg ==  max(macro$move.avg),]
simcountsHIGH <- A$simcount

A <- macro[macro$move.avg > 0.2,]
simcountsHIGH <- A$simcount



macro$HIGHMOVE <- macro$simcount %in% simcountsHIGH

title = "Mean and std. agents movements"
p.move <- ggplot(macro, aes(t, y = move.avg, colour=HIGHMOVE))
p.move <- p.move + geom_jitter()
p.move

title = "Mean and std. agents movements"
p.move <- ggplot(macro, aes(t, y = move.avg, colour=HIGHMOVE))
p.move <- p.move + geom_line()
p.move

title = "Mean and std. agents movements"
p.move <- ggplot(macro[macro$simcount == 15,], aes(t, y = move.avg))
p.move <- p.move + geom_line()
p.move



title = "Mean and std. agents movements"
p.move <- ggplot(macro, aes(t, y = move.avg, colour=HIGHMOVE))
p.move <- p.move + geom_smooth(size=2)
p.move


clu <- merge(params, A, by=c("simname","simcount","run"))

plot(1:nrow(clu), clu$R)

# it seems that most of the weird run have R = 0
# but some has not. let's check those

clu.Rnot0 <- clu[clu$R != 0,]

plot(1:nrow(clu.Rnot0), clu.Rnot0$alpha)

# those are the runs for alpha = 1

# so all the cases in which agents are not influenced by other agents


# other tests



#MOVEMENTS + SD
title = "Evolution of movements (+Std.Dev.) by sigma"
p.move <- ggplot(cl, aes(t))
p.move <- p.move + geom_line(aes(y = move.avg, colour=simname, group=simname))
p.move <- p.move + facet_grid(simname~.,margins=F)
p.move <- p.move + ggtitle(title) + xlab("Rounds") + ylab("Move")
if (INTERACTIVE) {
  p.move
}

# MOVEMENTS + SE
limits <- aes(ymax = cl$move.avg + cl$move.se, ymin=cl$move.avg - cl$move.se)
title = "Evolution of movements (+Std.Err.) by sigma"
p.move.se <- ggplot(cl, aes(t,group=simname))
p.move.se <- p.move.se + geom_point(aes(y = move.avg, colour=simname, group=simname))
#p.move.se <- p.move.se + geom_line(aes(y = move.avg, colour=simname, group=simname))
#p.move.se <- p.move.se + geom_errorbar(limits, width=0.2)
p.move.se <- p.move.se + ggtitle(title) + xlab("Rounds") + ylab("Move")
if (INTERACTIVE) {
  p.move.se
}
