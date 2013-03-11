library(plyr)
library(ggplot2)
library(gridExtra)
library(reshape)
library(MASS)
library(car)


## AXIS SCALES
plotScaleDis <- scale_y_discrete(breaks=seq(0, 1, 0.05))
plotScaleSize <- scale_y_discrete(breaks=seq(0, 100, 5))
plotScaleCount <- scale_y_discrete(breaks=seq(0, 100, 5))

plotXscale <- scale_x_discrete(breaks=seq(0, 30, 5))
reducedXScale <- scale_x_discrete(breaks=seq(5, 25, 10))

reducedYScaleCount <- scale_y_discrete(breaks=seq(0, 100, 10))
reducedYScaleDis <- scale_y_discrete(breaks=seq(0, 1, 0.2))
reducedYScaleSize <- scale_y_discrete(breaks=seq(0, 100, 10))

## X and Y LABS
yLabCount <- ylab('number of clusters')
yLabSize <- ylab('cluster size')
yLabDis <- ylab('distance from truth')

hs.setwd <- function(DIR, session){
  DATADIR = sprintf("%s%s/csv/", DIR, session)
  setwd(DATADIR)
  getwd()
}

hs.source <- function(sourcefile) {
 FULLSOURCE = sprintf("%s/%s", scriptdir, sourcefile)
 source(FULLSOURCE)
}



printParams <- function(varlist) {
  out <- '';
  for (n in names(params[c(-1:-9,-20:-length(params))])) {
    if (!(n %in% varlist)) {
      
      if (n == "init.vscaling") {
        N <- "V"
      }
      else {
        N <- n
      }

      if (sd(params[, n]) == 0) {
        V = params[1, n]
      }
      else {
        minN <- min(params[, n])
        maxN <- max(params[, n])
        if (class(params[,n]) != "integer") {
          minN <- sprintf('%2.3f',minN)
          maxN <- sprintf('%2.3f',maxN)
        }
        V = paste0("[", minN, "-", maxN, "]")
      }
     
      out <- sprintf('%s %s: %s, ', out, N, V)
    }
  }
  
  return(substr(out, 2, nchar(out)-2))
}

hs.makeggtitle <- function(text, vars) {
  paramString <- printParams(vars)
  return(ggtitle(paste0(text,"\n\n",paramString)))
}

saveOrPlot <- function(save, p, title, path) {
  if (save) {
    filename <- gsub(" ", "_", title)
    filename <- gsub("/n", "", filename)
    filename = tolower(paste0(filename,".jpg"))
    path <- paste0(path, filename)
    ggsave(filename=path, plot=p)
  }
  else {
    p
  }
}


# TIME EVOLUTION

plotTS2by2 <- function (v1, v2, save = TRUE) {
  
#count
  title <- paste0("Cluster counts in time by ", v1)
  p <- ggplot(clu, aes_string(x="t", y="count", group=v1, colour=v1))
  p <- p + geom_smooth() + yLabCount +
    plotScaleCount + plotXscale +
      hs.makeggtitle(title, c(v1))

  saveOrPlot(save, p, paste0("ts_", title), IMGPATH)
  
#size
  title <- paste0("Cluster sizes in time by ", v1)
  p <- ggplot(clu, aes_string(x="t", y="size.avg", group=v1, colour=v1))
  p <- p + geom_smooth() + yLabSize +
    plotScaleSize + plotXscale +
      hs.makeggtitle(title, c(v1))

  saveOrPlot(save, p, paste0("ts_", title), IMGPATH)

#from truth
  title <- paste0("Convergenge in time by ", v1)
  p <- ggplot(clu, aes_string(x="t", y="fromtruth.avg", group=v1, colour=v1))
  p <- p + geom_smooth() + yLabDis + plotXscale + # no plotScaleDis
      hs.makeggtitle(title, c(v1))

  saveOrPlot(save, p, paste0("ts_", title), IMGPATH)
}

boxplot2by2 <- function (v1, v2, save = TRUE) {

# BOXPLOTS

#count
  title <- paste0("Distribution of cluster counts by ", v1)
  p <- ggplot(clu, aes_string(x=v1, y="count", group=v1, colour=v1))
  p <- p + geom_boxplot() + yLabCount +
    plotScaleCount + 
      hs.makeggtitle(title, c(v1))

  saveOrPlot(save, p, paste0("boxplot_", title), IMGPATH)
  
#size
  title <- paste0("Distribution of cluster sizes by ", v1)
  p <- ggplot(clu, aes_string(x=v1, y="size.avg", group=v1, colour=v1))
  p <- p + geom_boxplot() + yLabSize +
    plotScaleSize + 
      hs.makeggtitle(title, c(v1))

  saveOrPlot(save, p, paste0("boxplot_", title), IMGPATH)

#from truth
  title <- paste0("Distribution of convergence levels by ", v1)
  p <- ggplot(clu, aes_string(x=v1, y="fromtruth.avg", group=v1, colour=v1))
  p <- p + geom_boxplot() + yLabDis + hs.makeggtitle(title, c(v1)) # no plotScaleDis

  saveOrPlot(save, p, paste0("boxplot_", title), IMGPATH)
    
}

#HEATMAP

heatmap2by2 <- function(v1, v2, save = TRUE) {
  
#count.avg  
  title <- paste0("Cluster counts by combinations of ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=count), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(clu$count)), guide="legend", breaks= seq(0,100,5), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))
  
  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)
  
#size.avg  
  title <- paste0("Cluster sizes by combinations of ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x=v1,y=v2))
  p <- p + geom_tile(aes(fill=size.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(clu$size.avg)), guide="legend", breaks= seq(0,100,5), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  
#fromtruth.avg  
  title <- paste0("Convergence levels by combinations of ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=fromtruth.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(0,max(clu$fromtruth.avg)), guide="legend", breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  # scale_fill_brewer() + hs.makeggtitle(title, c(v1, v2))
}


heatmap2by2Detail <- function(v1, v2, save = TRUE, data = clu) {

  
  
#count.avg  
  title <- paste0("Cluster counts by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=count), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$count)), breaks= seq(0,100,1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))
  
  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  
#size.avg  
  title <- paste0("Cluster sizes by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1,y=v2))
  p <- p + geom_tile(aes(fill=size.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$size.avg)), breaks= seq(0,100,1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

#size.sd
  title <- paste0("Std. dev. cluster sizes by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1,y=v2))
  p <- p + geom_tile(aes(fill=size.sd), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$size.sd)), breaks= seq(0,100,1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  
#fromtruth.avg  
  title <- paste0("Convergence levels by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=fromtruth.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(0,max(data$fromtruth.avg)), breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)


#fromtruth.avg  
  title <- paste0("Std. dev. of convergence levels by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=fromtruth.sd), colour = "white")
  p <- p + scale_fill_continuous(limits=c(0,max(data$fromtruth.sd)), breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  # scale_fill_brewer() + hs.makeggtitle(title, c(v1, v2))
}



## FACETS

facets2by2 <- function (v1, v2, save = TRUE) {

  facetFormula <- as.formula(sprintf('%s~%s', v2, v1))
  
#count ts
  title <- paste0("Cluster counts in time by ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x="t", y="count", group=v1, colour=v1))
  p <- p + geom_smooth()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedYScaleCount + reducedXScale + yLabCount
  p <- p  + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)

#count boxplot
  title <- paste0("Distributions of cluster counts by ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x=v1, y="count", group=v1, colour=v1))
  p <- p + geom_boxplot()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedYScaleCount + reducedXScale + yLabCount
  p <- p  + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)

#size ts
  title <- paste0("Cluster sizes in time by ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x="t", y="size.avg", group=v1, colour=v1))
  p <- p + geom_smooth()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedYScaleSize + reducedXScale + yLabSize
  p <- p + hs.makeggtitle(title, c(v1, v2))
  
  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)
  
#size boxplot
  title <- paste0("Distributions of cluster sizes by  ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x=v1, y="size.avg", group=v1, colour=v1))
  p <- p + geom_boxplot()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedYScaleSize + reducedXScale + yLabSize
  p <- p  + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)

#fromtruth ts
  title <- paste0("Convergence levels in time by ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x="t", y="fromtruth.avg", group=v1, colour=v1))
  p <- p + geom_smooth()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedXScale + yLabDis # no reducedYScaleDis
  p <- p  + hs.makeggtitle(title, c(v1, v2))
  
  saveOrPlot(save, p, paste0("facets_", title), IMGPATH)
  
#fromtruth boxplot
  title <- paste0("Distributions of convergence levels by  ", v1, " and ", v2)
  p <- ggplot(clu, aes_string(x=v1, y="fromtruth.avg", group=v1, colour=v1))
  p <- p + geom_boxplot()
  p <- p + facet_grid(facetFormula, margins = T)
  p <- p + reducedXScale + yLabDis # no reducedYScaleDis
  p <- p  + hs.makeggtitle(title, c(v1, v2))

  saveOrPlot(save, p, paste0("facets_", title), IMGPATH) 
}

allPlots <- function(v1, v2, save = TRUE) {
   plotTS2by2(v1,v2,save)
   boxplot2by2(v1,v2,save)
   facets2by2(v1,v2,save)
   heatmap2by2(v1,v2,save)
}
