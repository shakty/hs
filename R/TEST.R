# TEST
library(plot3D)

M <- mesh(seq(0, 6*pi, length.out = 80),
seq(pi/3, pi, length.out = 80))

u <- M$x ; v <- M$y
x <- u/2 * sin(v) * cos(u)
y <- u/2 * sin(v) * sin(u)
z <- u/2 * cos(v)

surf3D(x, y, z, colvar = z, colkey = FALSE, box = FALSE)


M <- mesh(seq(0, 2*pi, length.out = 80),
seq(0, 2*pi, length.out = 80))

u <- M$x ; v <- M$y
x <- sin(u)
y <- sin(v)

z <- sin(u + v)

surf3D(x, y, z, colvar = z, border = "black", colkey = FALSE)



# DUMPDIR 
DUMPDIR = "/mnt/tmp/dump/NAVNP/"

# Linear
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon/"
DIR = "attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v/"

INTERACTIVE = FALSE
PATH = paste0(DUMPDIR, DIR, "aggr/")
setwd(PATH)
IMGPATH <- paste0(PATH, "img/")

z <- cluall[cluall$sigma == 0.01, "fromtruth.avg"]

M <- mesh(sort(unique(cluall$R)), sort(unique(cluall$alpha)))
u <- M$x ; v <- M$y



surf3D(as.matrix(M$x), as.matrix(M$y), as.matrix(M$z))

# car
scatter3d(fromtruth.avg ~ R + alpha, data=cluall[cluall$t == 2000,])
