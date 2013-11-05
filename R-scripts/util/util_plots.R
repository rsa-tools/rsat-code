## ##############################################################
## Drawing and plotting utilities
##
## Loading:
##  source (file.path(dir.util, "util_plots.R"))

source(file.path(dir.util, 'util_descr_stats.R'))

################################################################
#
# Color palettes
#
################################################################


################################################################
## Define color palettes
green.to.red <- function() {
  red.levels <- c(rep(0,127), 0:127)/128
  green.levels <- c(127:0, rep(0,127))/128
  blue.levels <- 0
  green.to.red <- rgb(red.levels, green.levels, blue.levels)
  return(green.to.red)
}

blue.to.yellow <- function() {
  red.levels <- c(rep(0,127), 0:127)/128
  green.levels <-c(rep(0,127), 0:127)/128
  blue.levels <-  c(127:0, rep(0,127))/128
  blue.to.yellow <- rgb(red.levels, green.levels, blue.levels)
  return(blue.to.yellow)
}

blue.to.red <- function() {
  red.levels <-  c(seq(from=0,to=256,by=4),rep(256,64))/256
  green.levels <-  c(seq(from=0,to=256,by=4),seq(from=256,to=0,by=-4))/256
  blue.levels <-  c(seq(from=256,to=0,by=-4), rep(0,64))/256
  blue.to.red <- rgb(red.levels, green.levels, blue.levels)
  return(blue.to.red)
}

fire <- function() {
  red.levels <-  rep(1, 256)
  green.levels <-  c(seq(from=0,to=255,by=4/3),rep(256,64))/256
  blue.levels <- c(rep(0,192), seq(from=4,to=256,by=4))/256
  fire.palette <- rgb(red.levels, green.levels, blue.levels)
  return(fire.palette)
}

################################################################
## Green -> White -> Red palette
green.white.red <- function() {
  green.levels <- c(rep(127,127), 126:0)/128
  red.levels <- c(0:127,rep(127,127))/128
  blue.levels <- c(0:127, 126:0)/128
  palette <- rgb(red.levels, green.levels, blue.levels)
  return(palette)
}
################################################################
## Color palette mimicking the Blue->Yellow->Red palette from Berry et
## al. (2010)
palette.BYR <- function() {
  green.levels <- c(0:127, 126:0)/128
  red.levels <- c(0:127,rep(127,127))/128
  blue.levels <- c(127:0,rep(0,127))/128
  palette <- rgb(red.levels, green.levels, blue.levels)
  return(palette)
}


## associate a color to each group
get.palette <- function (x,
			 group.labels,
                         custom.palette = NULL
                         ) {
  if (is.null(custom.palette)) {
    custom.palette = palette
  }
  groups <-  levels(as.factor(group.labels))
  group.palette <- custom.palette[1:length(groups)]
  names(group.palette) <- groups
  return (group.palette)
}

################################################################
## Associate a color to each object on the basis of its group
get.colors <- function (x,group.labels,custom.palette=NULL) {
  n <- nrow(x)
  if (is.null(custom.palette)) {
    custom.palette = palette
  }
  group.palette <-  get.palette(x,group.labels,custom.palette)

  colors.all <- rep('#888888',length.out=n)
  names(colors.all) <- row.names(x)
  groups <-  levels(as.factor(group.labels))
  for (g in groups) {
    colors.all[group.labels==g] <- group.palette[g]
  }
  return (colors.all)
}

## Associate a color to each object on the basis of its group
get.colors.training <- function (x,group.labels) {
  colors.all <- get.colors(x,group.labels)
  return (colors.all[!is.na(group.labels)])
}


## ##############################################################
## Plot frequency polygons of several classes (one curve per class)
hist.by.class <- function (x,
                           group.labels,
                           breaks = pretty(x,40),
                           relative=F,
                           leg.pos=NA,
                           ... ### additional params are passed to plot and lines
                           ) {

  group.labels[is.na(group.labels)] <- "NA"
  
  x <- as.vector(x)
  x.min <- min(x,na.rm=T)
  x.max <- max(x,na.rm=T)
  groups <- unique(group.labels)


  h <- hist(x,breaks=breaks,plot=F)
  y <- h$counts
  if (relative) {
    y <- y / sum(y)
  }
  plot(h$mids,y,
       ylab="Counts",
       font.lab=2,
       col=1,
       lwd=2,
       type="n",
       ...
       )
  for (c in 1:length(groups)) {
    h <- hist( x[group.labels==groups[c]], breaks=breaks, plot=F)
    y <- h$counts
    if (relative) {
      y <- y / sum(y)
    }
    lines(h$mids,y,type="l",col=c+1,lwd=2,...)
    y.max <- max(y)
  }
  if (is.na(leg.pos)) {
    leg.pos = c((x.max+x.min)/2,y.max*0.9)
  }
  legend(leg.pos[1], leg.pos[2], legend=c(groups),col=(1:length(groups))+1,lwd=3)
}

################################################################
#
# Profile drawing
#
plot.profiles <- function(x,			# data frame
                          selection=NULL,		# a vector with selected row names from data
                          plot.type="l", # Can be changed to "b" for small selections
                          pch.per.profile=NULL, ## A vector specifying the pch for each profile. Must have the sema length as the number of profiles
                          main="profiles",	# main title for the plot
                          plot.legend = T,      # plot the legend
			  legend.labels = NULL,
                          ylim=range(x, na.rm=T),
			  xlab='',
			  ylab='',
			  lwd=1,
                          col.profiles="gray",		# a color (used for all the profiles) or a vector of colors (one per gene)
                          col.mean.profile="darkblue", # color of the mean profile 
                          col.median.profile="darkgreen", # color of the median profile 
                          plot.mean.profile=F, ## Plot the median profiles
                          plot.median.profile=T, ## Plot the median profiles
                          plot.sd.profile=T, ## plot curves at 1 std dev from the median
                          plot.quartile.profile=F, ## plot profiles of first (Q1) and third (Q3) quartiles
                          plot.conf.profile=F, ## Plot confidence interval around the median 
                          alpha=0.05, ## alpha risk for the confidence interval
                          X.axis=T, ## Plot the columns on the X axis
                          xlab.by = 1,
                          las=2,
                          ... ## Additional parameters are passed to plot()
                          ) {

  verbose("Plot profiles", 3)
  verbose(dim(x),3)
  verbose("Selection", 3)  
  verbose(length(selection), 3)
  verbose(selection, 4)  
  par(las=las)
  

  cols.for.axis=seq(from=1, to=ncol(x), by=xlab.by)
  
  if (is.null(legend.labels)) {
    legend.labels <- row.names(x)
  }
  legend.col <- vector()
  legend.lty <- vector()
  legend.lwd <- vector()
  legend.pch <- vector()
  
  x <- as.data.frame(x)

  if (is.null(selection)) {
    profiles <- x
  } else {
    if (length(selection) == 0) {
      return()
    } 
    profiles 	<- x[selection,]
  }
  if (!is.null(dim(profiles))) {
    p <- dim(profiles)[2]
    n <- dim(profiles)[1]
  } else {
    p <- 1
    n <- length(profiles)
  }
#   if (is.na(ylim)) {
#     ylim <- c(min(x,na.rm=T),max(x,na.rm=T))
#   }
  if (plot.legend) {
    xlim <- c(1,p*1.25)
  } else {
    xlim <- c(1,p)
  }

  palette.rep 	<- rep(col.profiles,length.out=n)
  plot(1:p,profiles[1,],
       main=main,
       xlim	= xlim,
       ylim	= ylim,
       type	= "n",
       xlab	= xlab,
       ylab	= ylab,
       lwd	= lwd,
       panel.first	= c(	grid(col=1),
         abline(h=0,col=1),
         abline(v=0,col=1)
         ),
       font.axis=2,
       font.lab=2,
       xaxt="n",
       ...
       )
  for (row in 1:n) {
    profile.pch <- 1
    if (!is.null(pch.per.profile)) {
      profile.pch <- pch.per.profile[row]
    }
    lines(1:p,profiles[row,],
          type	= plot.type,
	  lwd=lwd,
          col = palette.rep[row],
          pch = profile.pch)
  }

  if (X.axis) {
    axis(side=1, at=1:ncol(x), labels=NA,las=las) ## Draw tick bars for each column
    col.spacing <- round(ncol(x)/10)
    axis(side=1, at=cols.for.axis, labels=names(x)[cols.for.axis],las=las)
#    axis(side=1, at=1:ncol(x), labels=names(x),las=las)
  }
  
  ## Estimate central tendency and dispersion of the profiles
  profile.stats <- data.frame(
                              mean=apply(x, 2, "mean", na.rm=T), 
                              median=apply(x, 2, "median", na.rm=T),
                              sd=apply(x, 2, "sd", na.rm=T),
                              iqr=apply(x, 2, "IQR", na.rm=T),
                              Q1=apply(x, 2, "quantile", probs=0.25, na.rm=T),
                              Q3=apply(x, 2, "quantile", probs=0.75, na.rm=T)
                              )
  profile.stats$stderr <- profile.stats$sd/sqrt(nrow(x))
  profile.stats$conf <- profile.stats$sd*(-qnorm(alpha/2))/sqrt(nrow(x))
  
  ## Confidence interval about the mean profile
  profile.stats[,paste("confint", alpha, sep="")] <- profile.stats$stderr*qnorm(1-alpha/2)
#  setwd(dir.results); export.object(profile.stats, paste('profile_stats', group, sep="_"),export.formats="table")

  ## Plot mean profile
  legend.labels.first <- vector()
  if (plot.mean.profile) {
    lines(1:p, profile.stats$mean,col=col.mean.profile,lwd=2)
    legend.labels.first <- c(legend.labels.first, "mean profile")
    legend.col <- c(legend.col, "darkblue")
    legend.lty <- c(legend.lty, "solid")
    legend.lwd <- c(legend.lwd, 2)
    legend.pch <- c(legend.pch, NA)
  }

  ## Plot median profile
  if (plot.median.profile) {
    lines(1:p, profile.stats$median,col=col.median.profile,lwd=2)
    legend.labels.first <- c(legend.labels.first, "median profile")
    legend.col <- c(legend.col, "darkgreen")
    legend.lty <- c(legend.lty, "solid")
    legend.lwd <- c(legend.lwd, 2)
    legend.pch <- c(legend.pch, NA)
  }

  ## Plot sd estimate profile
  if (plot.sd.profile) {
    lines(profile.stats$median+profile.stats$sd, lwd=2, col="darkviolet",lty="dashed")
    lines(profile.stats$median-profile.stats$sd, lwd=2, col="darkviolet",lty="dashed")
    legend.labels.first <- c(legend.labels.first, "1 sd interval")
    legend.col <- c(legend.col, "darkviolet")
    legend.lty <- c(legend.lty, "dashed")
    legend.lwd <- c(legend.lwd, 2)
    legend.pch <- c(legend.pch, NA)
  }

  ## Plot quartile profile
  if (plot.quartile.profile) {
    lines(profile.stats$Q1, lwd=2, col="darkviolet",lty="dashed")
    lines(profile.stats$Q3, lwd=2, col="darkviolet",lty="dashed")
    legend.labels.first <- c(legend.labels.first, "1st and 3rd quartiles")
    legend.col <- c(legend.col, "darkviolet")
    legend.lty <- c(legend.lty, "dotted")
    legend.lwd <- c(legend.lwd, 2)
    legend.pch <- c(legend.pch, NA)
  }

  ## Plot profile of confidence interval
  if (plot.conf.profile) {
    lines(profile.stats$median+profile.stats$conf, lwd=2, col="darkred",lty="dotted")
    lines(profile.stats$median-profile.stats$conf, lwd=2, col="darkred",lty="dotted")
    legend.labels.first <- c(legend.labels.first, paste(1-alpha, "% confidence interval"))
    legend.col <- c(legend.col, "darkred")
    legend.lty <- c(legend.lty, "dotted")
    legend.lwd <- c(legend.lwd, 2)
    legend.pch <- c(legend.pch, NA)
  }

  ## add a legend to the plot
  if (plot.legend) {
    legend.labels <- c(legend.labels.first, legend.labels)
    legend.col <- c(legend.col, palette.rep)
    legend.lty <- c(legend.lty, rep("solid", length.out=n))
    legend.lwd <- c(legend.lwd, rep(lwd, length.out=n))
    if ((!is.null(pch.per.profile)) || (plot.type == "b")) {
      legend.pch <- c(legend.pch, pch.per.profile)
    } else {
      legend.pch <- NULL
    }
    if (length(legend.labels) > 25) {
      legend.labels <- c(legend.labels[1:25], "...")
      legend.col <- c(legend.col[1:25], "red")
      legend.lty <- c(legend.lty[1:25], "dotted")
      legend.lwd <- c(legend.lwd[1:25], 2)
      legend.pch <- c(legend.pch[1:25], 2)
    }
    legend("topright",legend.labels,col=legend.col,lty=legend.lty, lwd=legend.lwd,pch=legend.pch, bg='white', bty='o')
  }
  return(profile.stats)
}



################################################################
##
## Plot frequency distributions from a data frame.
## One polygon frequency per column
##
plot.frequency.distrib <- function(x,
                                   main = "frequency distribution",
                                   colnames  = names(x), # restrict the plot to selected columns
                                   col = rainbow(ncol(x)), # vector of colors
                                   class.number = 50,
                                   class.interval = NA,
                                   distrib.min = min(x, na.rm=T),
                                   distrib.max = max(x, na.rm=T),
                                   xlab = "values",
                                   ylab = "frequencies",
                                   xlim   = range(x,na.rm=T),
                                   ylim   = NA,
                                   plot.legend = T,  # plot the legend
                                   single.plot = TRUE, # when FALSE, one separate plot is drawed per variable
                                   draw.mean  = FALSE, #
                                   draw.median = FALSE,
                                   draw.mode  = FALSE,
                                   cumul  = FALSE
                                   ) {

  x <- na.omit(data.frame(x))
  col.number <- length(colnames)
  if (is.na(col)) {
    colors <- rep(rainbow(255),length.out=col.number)
  } else {
    colors <- rep(col,length.out=col.number)
  }
  
  ## calculate class frequencies
  breaks <- pretty(c(distrib.min,distrib.max),n=class.number)
  class.number <- length(breaks)-1
  curves <- matrix(nrow=class.number,ncol=dim(x)[2])
  
  for (i in 1:col.number) {
##     distrib <- class.grouping(as.vector(x[,i]),
##                               breaks <- breaks,
##                               class.interval=class.interval,
##                               display.plot=FALSE,
##                               return="descr")
    h <- hist(as.vector(x[,i]),breaks=breaks,plot=F)
    mids <- h$mids
    freq <- h$counts/sum(h$counts)
    
    if (cumul) {
      curves[,i] <- cusum(freq)
    } else {
      curves[,i] <- freq
    }
  }

  ## calculate limits of X and Y axis 
  if (is.na(ylim)) {
    if (cumul) {
      ylim <- c(0,1)
    } else {
      ylim <- c(0,max(curves,na.rm=T))
    }
  }

  ## plot the class frequencies
  if (single.plot) {
    ## a single plot with all the frequency polygons
    plot(mids, curves[,1],
         type = "l",
         xlim = xlim,
         ylim = ylim,
         xlab = xlab,
         ylab = ylab,
         main = main,
         col = colors[1],
         lwd = 2,
         cex = 2,
         font.axis = 2,
         font.main = 4,
         font.lab = 2,
         panel.first = c(grid(col=1),
           abline(h=0,col=1),
           abline(v=0,col=1))
         )
    for (i in 2:col.number) {
      lines(mids, curves[,i],
            type = "l",
            col = colors[i],
            lwd = 2)
    }
    ## draw the legend
    if (plot.legend) {
      legend.x <- xlim[1]# + 0.75*(xlim[2] - xlim[1])
      legend.y <- ylim[2]
      legend(legend.x,legend.y, colnames, col=colors,
             lty=1,lwd=2,cex=0.8,bg='white',bty='o')
    }
  } else {
    ## one separate plot per frequency polygon
    par(mfrow=(n2mfrow(col.number)))
    for (i in 1:col.number) {
      plot(distrib$mid, curves[,i],
           type = "l",
           xlim = xlim,
           ylim = ylim,
           xlab = xlab,
           ylab = ylab,
           main = colnames[i],
           col = colors[i],
           lwd = 2,
           cex = 2,
           font.axis = 2,
           font.main = 4,
           font.lab = 2,
           panel.first = c(grid(col=1),
             abline(h=0,col=1),
             abline(v=0,col=1))
           )
    }
    par(mfrow=c(1,1))
  }
  if (draw.mode) {
    mode <- distrib[which.max(distrib$freq),"mid"]
    print(paste("mode = ",mode))
    abline(v=mode,lty=1,col=colors[i])
  }
  if (draw.mean)  {
    print(paste("mean = ",mean(vect,na.rm=T)))
    abline(v=mean(vect,na.rm=T),col=colors[i])
  }
  if (draw.median) {
    print(paste("median = ",median(vect,na.rm=T)))
    abline(v=median(vect,na.rm=T),col=colors[i])
  }
}


## ##############################################################
## Plot the data in the plane of the two principal components
## with a specific color for each class
plot.prcomp.groups <- function(x, ### A data frame with the multivariate data
			       group.labels, ## A vector with one label per row of the data frame
			       group.palette = NA, ## vector with one color per group
			       plot.NA = TRUE, ## plot objects without label
			       xleg = NA, ## X position of the legend box
			       yleg = NA, ## Y position of the legend box
			       group.symbols=NA, ## Specific symbols for plotting each point
			       label.ref = T, ## calculate rotation with labelled objects only
			       scale.axes = F, ## rescale the axes before rotation
			       return.pc = F, ## Return the result of prcomp
			       ...
			       ) {
    ## Create a temporary vector for group.labels 
    group.labels.tmp <- group.labels

    ## Define a group palette if not provide
    if (is.na(group.palette)) {
        group.palette <- get.palette(x,group.labels.tmp)
    }

    ## group names
    groups <- names(table(group.labels.tmp))

    ## Temporarily convert NA to labels
    if (plot.NA) {
	group.labels.tmp[is.na(group.labels.tmp)] <- "undef"
	group.palette["undef"] <- "#AAAAAA"
	groups <- c("undef",groups)
    }

    object.colors <- group.palette[group.labels.tmp]
    names(object.colors) <- names(group.labels.tmp)

    ## calculate principal components
    if (label.ref) {
	## calculate the PC on the labeled objects only
	verbose("Calculating principal components on the labeled objects only, using correlation matrix", 2) 
	pc <- prcomp(x[!is.na(group.labels),])
        verbose("Rotating and scaling the whole data set", 2) 
        if (scale.axes) {
	    x.scaled <- scale(x,center=T,scale=T)
            pc$x <- as.matrix(x.scaled) %*% pc$rotation
	} else {
            pc$x <- as.matrix(x) %*% pc$rotation	    
	}
    } else {
        pc <- prcomp(x)
    }

    pc.comp1 <- pc$x[,"PC1"]
    pc.comp2 <- pc$x[,"PC2"]

    if (is.na(xleg)) { xleg <- min(pc.comp1) }
    if (is.na(yleg)) { yleg <- max(pc.comp2) }

    
    ## group.symbols
    if (is.na(group.symbols)) {
	group.symbols=19
    }

    ## plot the genes on the pla made of 2 first components
    plot(pc.comp1,
	 pc.comp2,
	 type	= "p",
	 pch	= group.symbols,
	 col	= object.colors,
	 xlab	= "First component", 
	 ylab	= "Second component",
	 panel.first	= c(grid(col=1),
			    abline(h=0,col=1),
			    abline(v=0,col=1)
			    ),
	 ...
	 )

    ## Redraw the objects with a real class to avoid them being
    ## hidden under the multitude of NA
    if (plot.NA) {
	lines(pc$x[!is.na(group.labels),"PC1"],
	      pc$x[!is.na(group.labels),"PC2"],
	      type	= "p",
	      pch	= 19,
	      col	= object.colors[!is.na(group.labels)],
	      )
    }
    legend(xleg,yleg,groups,col=group.palette[groups],pch=19,cex=1)
    rm(group.labels.tmp)
    if (return.pc) {
	return(pc)
    }
}


## ##############################################################
## Starting from a data frame (patttern amtrix), calculate coefficient
## of correlation between each pair of objects (row) and apply
## multidimensional sclaing with SVD to obtain a 2D plot.
plot.svd <- function(x, ## a data frame
		     col='gray', ## a vector of colors
		     pch=19,
                     xlab = "SVD1",
                     ylab = "SVD2",
		     ... ## Additional parameters are passed to the plot function
		     ) {

    c <- cor(t(x), use="pairwise")
    s <- svd(as.matrix(c))
        
    plot(s$u[,1],s$u[,2],
	 col = col,
	 type = "p",
	 xlab = xlab,
	 ylab = ylab,
	 pch=pch,
	 panel.first = c(	grid(col=1),
           abline(h=0,col=1),
           abline(v=0,col=1)),
	 ...
        )
}

################################################################
#
# perform SVD and plot data according to the 2 first dimensions
#
# TO BE DEBUGGED
plot.svd.2 <- function(x, ### a data frame
		     group.labels=NULL, ## A vector with one label per row of the data frame
                     object.colors=NULL, ### a vector of colors (one per object if the data frame)
                     sim=NA, ### a (dis)similarity matrix
		     ...
                     ) {
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  if (is.null(object.colors)) {
      if (is.null(group.labels)) {
	  object.colors <- "#ff0088"
      } else {
          object.colors <- get.colors(x,group.labels)
      }
  }
#  object.colors <- group.palette[group.labels]
  
  ## Calculate similarity matrix if required
  if (is.na(sim)) {
    sim <-  cor(t(x), use="pairwise.complete.obs")
  }

  ## Singular value decomposition
  s <-  svd(as.matrix(sim))

  ## Plot the result
  par(mai=c(0.8,0.8,0.2,0.2))
  plot(s$u[,1],s$u[,2],
       col = group.palette,
       type = "p",
       xlab = "SVD1",
       ylab = "SVD2",
       panel.first = c(	grid(col=1),
         abline(h=0,col=1),
         abline(v=0,col=1)),
       ...

       )
  par(mai=c(1,1,1,1))
  
  return(s)
}


################################################################
## Draw a dotplot from a matrix
dotplot.image <- function (x, ## a matrix with values to use for dotplot
                           col=NULL, ## Default colors: black and white
                           ...
                           ) {

  ## Define monochrome palette
  if (is.null(col)) {
    x.max <- max(x)
    col <- rgb(x.max:0,x.max:0,x.max:0,maxColorValue=x.max)
  }
                                        #  default.mai <- par("mai")
#  par(mai=c(1.4,1.2,0.5,0.2))
  image(as.matrix(x),col=col,axes=F,...)
  box()
  grid(nrow(x),ncol(x),col='#000000')
  nc1 <- ncol(x) - 1
  nr1 <- nrow(x) - 1
  axis(1,labels=row.names(x),at=(0:nr1)/nr1,las=2)
  axis(2,labels=names(x),at=(0:nc1)/nc1,las=2)
#  par(mai=default.mai)
}


## ##############################################################
## Plot a matrix with the values as labels, and a label text size
## proportional to the values
plot.table <- function (x, # a table object, such as those obtained with the function table()
                        ... ## additional parameters are passed to plot()
                        ) {
  
  d <- data.frame(x)
  print(d)
#  print(as.numeric(as.vector(d[,1])))
  d1 <-as.numeric(as.vector(d[,1])) 
  d2 <-as.numeric(as.vector(d[,2])) 
  
  plot(d1,d2,
       type="n", ...)
  
  f <- d[d$Freq > 0, "Freq"]
  c <- 2*(f/max(f))^0.25
  
#  d1 <- as.numeric(as.vector(d[d$Freq > 0,1]))
#  d2 <- as.numeric(as.vector(d[d$Freq > 0,2]))
  text(d1[d$Freq > 0], d2[d$Freq > 0], labels=f, cex=c)
}


## ##############################################################
## plot a Hinton diagram
##
## The Hinton diagram provides an intuitive representation of a weight
## matrix.
##
## Each cell of the matrix is represented by a circle whose area is
## proportional to the cell value.
##
## Distinct colors are used for positive and negative
## values.
##
plot.hinton.diagram <- function (x, ## A matrix or data frame witht he weigts to be plotted
                                 line.colors=c("black","grey"), ## line colors for circles associated to positive and negative value, resp. 

                                 fill.colors=c("black","white"), ## fill colors for circles associated to positive and negative value, resp. 
                                 scale=1, ## Factor for scaling the radius of the circles
                                 axes=c(1,2), ## Axes where the legends have to be plotted.
                                              ## Same numbering as the function axis(): 1=below, 2=left, 3=above and 4=right.
                                 shape="rectangle", ## supported: circle | square | rectangle | bars
                                 return.result = F, ## If true, a detailed result is returned as a data frame
                                 xlab=NA,ylab=NA, ## By default, no labels on the axes
                                 ... ## Additional parameters are passed to the plot function symbols()
                                 ) {

  ## Reverse the order ot the columns in order to o from top to bootom rather than from bottom to top
  x <- x[nrow(x):1,]
  
  ## Compute the number of columnx and X, Y positions of the circles
  nr <- nrow(x)
  nc <- ncol(x)
  x.centers <- ((1:nc)-0.5)/nc
  y.centers <- ((1:nr)-0.5)/nr

  ## Cast data frame of values to a vector
  value <- as.vector(as.matrix(t(x)))
  

  ## Define sign-specific colors
  x.sign <- sign(value)+1
  x.sign[x.sign==2] <- 1
  circle.line.colors <- line.colors[2-x.sign]
  circle.fill.colors <- fill.colors[2-x.sign]

  ## Compute a vector of X and Y coordinates for the circles
  centers <- expand.grid(x.centers,y.centers)
  names(centers) <- c("x","y")
  
  ## Rectangle sizes
  height <- scale*sqrt(abs(value))/(max(sqrt(abs(value)))*nr)
  width <- scale*sqrt(abs(value))/(max(sqrt(abs(value)))*nc)
  rect <- as.matrix(data.frame(width, height))
  min.side <- apply(rect,1,min)
  squares <- as.matrix(data.frame(min.side, min.side))
  
  ## Bar sizes
  if (shape == "bars") {
    bar.height <- scale*abs(value)/(max(abs(value))*nr)
    bar.width <- scale/nc
    rect <- as.matrix(data.frame(bar.width, bar.height))
  }
  
  ## Inch dimensions (necessary for circle)
  din <- data.frame(width*par("din")[1],height*par("din")[2])
  names(din) <- c("width","height")
  side.inches <- apply(din,1,min) ## For squares, the function symbol() uses the X dimensions

  ## Draw the circles
  if (shape == "circle") {
    symbols(x=centers[,1], y=centers[,2], fg=circle.line.colors, bg=circle.fill.colors,new=F,xlab=xlab, ylab=xlab,xaxt="n",yaxt="n", inches=max(side.inches/2), circles=side.inches/2, ...)
  } else if (shape == "square") {
    symbols(x=centers[,1], y=centers[,2], fg=circle.line.colors, bg=circle.fill.colors,new=F,xlab=xlab, ylab=ylab,xaxt="n",yaxt="n", inches=max(side.inches), squares=side.inches, ...)
    ## The square function is tricky, I simply use a rectangle with square coordinates
#    symbols(x=centers[,1], y=centers[,2], fg=circle.line.colors, bg=circle.fill.colors, new=F,xlab=xlab, ylab=ylab,xaxt="n",yaxt="n", inches=F, rectangles=squares,...)
  } else if ((shape == "rectangle") || (shape == "bars")){
    symbols(x=centers[,1], y=centers[,2], fg=circle.line.colors, bg=circle.fill.colors, new=F,xlab=xlab, ylab=ylab,xaxt="n",yaxt="n", inches=F, rectangles=rect,...)
  } else {
    stop('Invalid shape for plot.hinton.diagram(): only "circle" and "square" are supported.')
  }
  
  ## Draw the axes
  if (1 %in% axes) { axis(1,labels=names(x),at=x.centers,las=2) }
  if (2 %in% axes) { axis(2,labels=row.names(x),at=y.centers,las=2) }
  if (3 %in% axes) { axis(3,labels=names(x),at=x.centers,las=2) }
  if (4 %in% axes) { axis(4,labels=row.names(x),at=y.centers,las=2) }

  ## Compute a data frame with the results if required
  if (return.result) {
    result <- data.frame(value=value,
                         centers=centers,
                         rect=rect,
                         din=din,                       
                         side.inches=side.inches,
                         line.col=circle.line.colors,
                         fill.col=circle.fill.colors
                         )
    return(result)
  } 
}
