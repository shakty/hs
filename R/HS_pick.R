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

INTERACTIVE = FALSE
PATH = paste0(DUMPDIR,DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");


params <- read.table('params_all.csv', head=TRUE, sep=",")
#params$simname <- as.factor(params$simname)
#params$simcount <- as.factor(params$simcount)
#params$run <- as.factor(params$run)

macro <- read.table('clusters_macro_all.csv', head=TRUE, sep=",")
#macro$simname <- as.factor(macro$simname)
#macro$simcount <- as.factor(macro$simcount)
#macro$run <- as.factor(macro$run)
# Merging macro and params: clu
clu <- merge(params, macro, by=c("simname","simcount","run"))
# Factorising
#for (n in names(clu[1:25])) {
#  clu[, n] <- as.factor(clu[, n])      
#}


AA <- clu[clu$sigma == 0.1 & clu$R == 0.01 & clu$t == 100 & clu$alpha == 0.01,]$simcount ; AA
