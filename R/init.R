library(plyr)
library(ggplot2)
#library(gridExtra)
library(reshape)
library(MASS)
library(car)
library(plm)

## AXIS SCALES
plotScaleDis <- scale_y_discrete(breaks=seq(0, 1, 0.05))
plotScaleSize <- scale_y_discrete(breaks=seq(0, 100, 5))
plotScaleCount <- scale_y_discrete(breaks=seq(0, 100, 5))

plotXscale <- scale_x_discrete(breaks=seq(0, 30, 5))
reducedXScale <- scale_x_discrete(breaks=seq(5, 25, 10))

RXScale <- scale_x_discrete(breaks=seq(0,100,0.05))

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


lagg <- function(series, nlag=1) {
  return(c(rep(NA,nlag),series[1:obs-nlag]))
}

printParams <- function(varlist, paramsData = params) {
  out <- '';
  for (n in names(paramsData[c(-1:-9,-20:-length(paramsData))])) {
    if (!(n %in% varlist)) {
      
      if (n == "init.vscaling") {
        N <- "V"
      }
      else {
        N <- n
      }

      if (sd(paramsData[, n]) == 0) {
        V = paramsData[1, n]
      }
      else {
        minN <- min(paramsData[, n])
        maxN <- max(paramsData[, n])
        if (class(paramsData[,n]) != "integer") {
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

hs.makeggtitle <- function(text, vars, paramsData = params) {
  paramString <- printParams(vars, paramsData)
  # can't be longer than 2 linesc
  limit = 30
  if (nchar(paramString) > limit) {
    str <- paste0(substr(paramString, 1, limit), "\n", substring(paramString, limit+1))
  }
  str <- paste0(text,"\n\n", paramString)
  return(ggtitle(str))
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


heatmap2by2Detail <- function(v1, v2, save = TRUE, data = clu, paramsData = params) {

  
#count.avg  
  title <- paste0("Cluster counts by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=count), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$count)), breaks= seq(0,100,5), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2), paramsData) + RXScale
  
  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  
#size.avg  
  title <- paste0("Cluster sizes by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1,y=v2))
  p <- p + geom_tile(aes(fill=size.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$size.avg)), breaks= seq(0,100,10), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2), paramsData) + RXScale

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

#size.sd
  title <- paste0("Std. dev. cluster sizes by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1,y=v2))
  p <- p + geom_tile(aes(fill=size.sd), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$size.sd)), breaks= seq(0,100,10), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2), paramsData) + RXScale

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)

  
#fromtruth.avg  
  title <- paste0("Convergence levels by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=fromtruth.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(0,max(data$fromtruth.avg)), breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2), paramsData) + RXScale

  saveOrPlot(save, p, paste0("heat_", title), IMGPATH)


#fromtruth.avg  
  title <- paste0("Std. dev. of convergence levels by combinations of ", v1, " and ", v2)
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=fromtruth.sd), colour = "white")
  p <- p + scale_fill_continuous(limits=c(0,max(data$fromtruth.sd)), breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2), paramsData) + RXScale

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


## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                           c( N    = length2(xx[,col], na.rm=na.rm),
                              mean = mean   (xx[,col], na.rm=na.rm),
                              sd   = sd     (xx[,col], na.rm=na.rm)
                              )
                          },
                    measurevar,
                    na.rm
             )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

heatmapFacets<- function(v1,v2,v3,data = clu, paramsData = params, imgpath = IMGPATH, save = TRUE) {

  pt <- heatmapFacets_fromtruth(v1,v2,v3)
  saveOrPlot(save, pt$p, paste0("hf_", pt$t), imgpath)

  pt <- heatmapFacets_count(v1,v2,v3)
  saveOrPlot(save, pt$p, paste0("hf_", pt$t), imgpath)

  pt <- heatmapFacets_size(v1,v2,v3)
  saveOrPlot(save, pt$p, paste0("hf_", pt$t), imgpath)

}

heatmapFacets_fromtruth<- function(v1,v2,v3,data = clu, paramsData = params, imgpath = IMGPATH, save = TRUE) {

  title <- paste0("Convergence levels by combinations of ", v1, ", ", v2, ", and ", v3)
  facetFormula <- as.formula(sprintf('%s~.', v3))
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=fromtruth.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(0,max(data$fromtruth.avg)), breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2, v3), paramsData) + RXScale + scale_y_discrete(breaks= seq(0,100,0.1)) 
  p <- p + facet_grid(facetFormula)

  return(list(p=p,t=title))
}

heatmapFacets_count<- function(v1,v2,v3,data = clu, paramsData = params, imgpath = IMGPATH, save = TRUE) {

  title <- paste0("Cluster counts by combinations of ", v1, ", ", v2, ", and ", v3)
  facetFormula <- as.formula(sprintf('%s~.', v3))
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=count), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$count)), breaks= seq(0,1,0.1), high='lightblue',low='red')
  p <- p + hs.makeggtitle(title, c(v1, v2, v3), paramsData) + RXScale + scale_y_discrete(breaks= seq(0,100,0.1)) 
  p <- p + facet_grid(facetFormula)

  return(list(p=p,t=title))
}

heatmapFacets_size<- function(v1,v2,v3,data = clu, paramsData = params, imgpath = IMGPATH, save = TRUE) {

  title <- paste0("Cluster sizes  by combinations of ", v1, ", ", v2, ", and ", v3)
  facetFormula <- as.formula(sprintf('%s~.', v3))
  p <- ggplot(data, aes_string(x=v1, y=v2))
  p <- p + geom_tile(aes(fill=size.avg), colour = "white")
  p <- p + scale_fill_continuous(limits=c(1,max(data$size.avg)), breaks= seq(0,1,0.1), low='lightblue',high='red')
  p <- p + hs.makeggtitle(title, c(v1, v2, v3), paramsData) + RXScale + scale_y_discrete(breaks= seq(0,100,0.1)) 
  p <- p + facet_grid(facetFormula)

  return(list(p=p,t=title))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
                
    pushViewport(viewport(layout = grid.layout(nrow(layout)+1, ncol(layout))))

    grid.text("title of this panel", vp = viewport(layout.pos.row = 1, ))
 
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row + 1,
                                      layout.pos.col = matchidx$col))
    }
  }

}


tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")) {
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function() {
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}
