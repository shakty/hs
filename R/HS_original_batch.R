# HS analysis
source("/opt/MATLAB_WORKSPACE/hs/R/init.R")

DIR = "R-alpha-noA-noB-2013-3-6-20-8/"
DIR = "A-alpha-R-2013-3-10-10-25/"

DIR = "tau-sigma-by-alpha-R-noA-noB-truth_middle-2013-3-9-18-42/"

DIR = "A-B-alpha-R-tau-2013-3-10-0-29/"

DIR = "sigma_tau-2013-3-6-10-36/"



DIR = "vscaling-2013-3-10-22-3/"
DIR = "vscaling-tau-easy-to-converge-2013-3-11-9-31/"

DIR = "R-alpha-01-2013-3-10-21-4/"

DIR = "R-alpha-99-2013-3-10-21-37/"

DIR ="bottom-R-alpha-2013-3-11-10-33/"

DIR = "upper-R-alpha-2013-3-11-17-45/"

DIR = "limited_sigma_R/"

DIR = "limited_sigma_R_matlab13/"

DIR = "limited_sigma_R_sequantial_random_avv1/"

DIR = "limited_sigma_R_seq_rnd_avv1/"
DIR = "limited_sigma_R_seq_rnd_avv0/"

DIR = "limited_sigma_R_sim_avv0_B/"
DIR = "limited_sigma_R_sim_avv1/"

DIR = "limited_sigma_R_seq_det_avv0/"
DIR = "limited_sigma_R_seq_det_avv1/"

DIR = "limited_sigma_R_sim_avv0_SEED/"

DIR = 'alpha1_tau_vinit_av1/'
DIR = 'alpha1_tau_vinit_av0/'

DIR = "cluster_zone_sigma_R_alpha_av1/"
DIR = "cluster_zone_sigma_R_alpha/"

DIR = "truth_corner_alpha_Rright_sigma_av0/"
DIR = "truth_corner_alpha_Rright_sigma_av1/"

DIR = "truth_corner_alpha_R_av0/"

DIR = "truth_corner_alpha_Rmiddle_sigma_av1/"

# too big
#DIR = "truth_corner_alpha_R_sigma_av1_ALL/"

DIR = "truth_corner_alpha_R_av1/"

DIR = "noisev_truth_corner_alpha_Rleft_sigma_av1/"
DIR = "noisev_truth_middle_alpha_Rleft_sigma_av1/"
DIR = "noisev_truth_exactcorner_alpha_Rleft_sigma_av1/"

DIR = "simul_tmiddle_noisev_av1_rleft/"
DIR = "simul_tcorner_noisev_av1_rleft/"
DIR = "simul_texactcorner_noisev_av1_rleft/"

DIR = "attr0_av1_nv_seqrnd_thm/"

DIR = "tec_np_seqrnd_av1_Rleft/"

DUMPDIR = "/opt/MATLAB_WORKSPACE/hs/dump/"
PATH = paste0(DUMPDIR,DIR)
setwd(PATH)
IMGPATH <- paste0(PATH, "img/");
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

# START


v1 <- "R"
v2 <- "alpha"
v3 <- "sigma"
#data <- clu
data <- clu[clu$t == 21,]
paramsData <- params
heatmapFacets(v1,v2,v3, data)


OLDPATH = IMGPATH
for (S in unique(params$sigma)) {
  curDir <-  paste0("sigma_",S,"/")
  dir.create(file.path(OLDPATH, curDir), showWarnings = FALSE)
  IMGPATH <- paste0(OLDPATH, curDir)
  heatmap2by2Detail(v1,v2, data = clu[clu$sigma == S,], paramsData = params[params$sigma == S,])
}
IMGPATH = OLDPATH
