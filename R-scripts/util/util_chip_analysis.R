################################################################
##
## diverse functions for analysing gene expression data
##
## loading:
## 	source(file.path(dir.util, "util_chip_analysis.R"))
##
## demo
##	chip.util.demo()
library(stats)

source(file.path(dir.util,'palette.R'))
source(file.path(dir.util,'util_plots.R'))

################################################################
## Take as input a set of expression profiles + clusters defined on
## this set, and plot the profiles per cluster
plot.cluster.profiles <- function(x, ## A  N x P data frame
                                  clusters, ## Cluster assignation, a vector with N enties (N is hte number of rows of the profile table)
                                  plot.pca = T, ## Plot clusters on a PCA projection
                                  main.PCA.var='Var per PC', ## Main title for the PCA variance histogram
                                  main.PCA.plot='Cluster mapping', ## Main title for the PCA variance histogram
                                  cluster.col=NULL, ## a vector specifying the k colors associated to the respective clusters
                                  ... ## Additional parameters are passed to the plot function
                                  ) {
  
  
  

  ## Number of clusters
  cluster.ids <- sort(unique(clusters))
  k <- length(cluster.ids)
  
  ## Specify cluster colors
  if (is.null(cluster.col)) {
    cluster.col <- rainbow(n=k)
    names(cluster.col) <- cluster.ids
  }
  
  ## Defien the number of plots
  if (plot.pca) {
    plots <- k+2 ### one plot for the 2D visualization + one per cluster
  } else {
    plots <- k
  }
  plot.rows <- max(1,floor(sqrt(plots)))
  plot.cols <- max(1,ceiling(plots/plot.rows))
  par(mfrow=c(plot.rows,plot.cols))
  
  ## ############## k-mean clustering
  if (plot.pca) {
    ## ## prior transformation in principal compenent space 
    pcomp <- prcomp(x)    
    plot(pcomp, main=main.PCA.var)
    plot(pcomp$x[,"PC1"],pcomp$x[,"PC2"],
         col = cluster.col[clusters],
         main=main.PCA.plot,
         xlab = "PC1", 
         ylab = "PC2",
         panel.first = c(grid(col='#BBBBBB',lty='solid'),
           abline(h=0,col=1),
           abline(v=0,col=1)),
         ...
         )
  }
  
  ## ############## Profile plotting
  for (cl.nb in 1:k) {
    current.cluster <- cluster.ids[cl.nb]
    current.profiles <- x[clusters==current.cluster,]
    
    plot.profiles(current.profiles,
                  plot.legend=F,
                  main=paste("cl",cl.nb,"(n=",nrow(current.profiles),")"),
                  col.median.profile=cluster.col[cl.nb],
                  ... 
                  )

  }
  par(mfrow=c(1,1))
}



################################################################
#
# Perform k-means clustering 
# and plot resulting profiles by cluster
#
kmean.profiles	<- function(x, ## A data frame
                            k = 19, ## Number of clusters
                            iter = 200, ## Number of iterations
                            pca = F, ## PCA transformation before analyzing
			    pch=19,
                            cluster.col=rainbow(n=k), ## a vector specifying the k colors associated to the respective clusters
			    ... ## Additional parameters are passed to the plot function
			    ) {

  plots <- k+1 ### one plot for the 2D visualization + one per cluster
  plot.rows <- max(1,floor(sqrt(plots)))
  plot.cols <- max(1,ceiling(plots/plot.rows))
  par(mfrow=c(plot.rows,plot.cols))
  par(mai=c(0.3,0.5,0,0.05))

  nrow  <- dim(x)[1]
  ncol  <- dim(x)[2]

  ## ############## k-mean clustering
  if (pca) {
    ## ## prior transformation in principal compenent space 
    pcomp	 <- prcomp(x)
    km <- kmeans(pcomp$x,k,iter)
    
    plot(pcomp$x[,"PC1"],pcomp$x[,"PC2"],
         col = cluster.col[km$cluster],
         xlab = "PC1", 
         ylab = "PC2",
         pch=pch,
         panel.first = c(grid(col=1),
           abline(h=0,col=1),
           abline(v=0,col=1)),
         ...
         )
  } else {
    km <- kmeans(x,k,iter)
    ## select first and last columns for assigning coordinates
    plot(x[,1],
         x[,ncol],
         xlab = names(x)[1],
         ylab = names(x)[ncol],
         col = cluster.col[km$cluster],
         pch=pch,
         panel.first = c(	grid(col=1),
           abline(h=0,col=1),
           abline(v=0,col=1)),
         ...
         )
  }

  ## ############## Profile plotting
  for (cl.nb in 1:k) {
    cl <- x[km$cluster==cl.nb,]
    v <- 1:(dim(cl)[2])
    plot(v,
         cl[1,],
         ylim = c(min(x),max(x)),
         type = "l",
         xlab = "column",
         ylab = "value",
         col = cluster.col[cl.nb],
         panel.first = c(	grid(col=1),
           abline(h=0,col=1),
           abline(v=0,col=1)
           ),
         font.axis=2,
         font.lab=2,
         ...
         )
    for (row in 2:dim(cl)[1]) {
      lines(v,cl[row,],
            type = "l",
            col = cluster.col[cl.nb])
    }
    legend(min(v),max(x),paste("cl",cl.nb,"(",length(cl[,1]),"el)"), bty = "n",cex=1.5)
  }	
  par(mai=c(1,1,1,1))
  par(mfrow=c(1,1))

  return(km)
}




################################################################
## Plot profiles of median and IQR values per column
plot.median.IQR.profiles <- function(data.profiles, 
				     time.points=NULL, ## Time points. If not provided, the time interval is assumed constant
				     data.type="", ## A string used as prefix for the plot title
				     data.colors='#0000FF', ## Colors associated to each point

				     ## impose limits on the Y axis
				     median.ymin=NULL,
				     median.ymax=NULL,
				     IQR.ymin=0,
				     IQR.ymax=NULL,
				     ...) {

  ## interval
  if (is.null(time.points)) {
      time.points=1:ncol(data.profiles)
  }

  ################################################################
  ## Median profile

  ## Calculate median profiles
  (data.col.median <- apply(data.profiles, 2, median,na.rm=T))
  
  ## Calculate Y limits if not provided
  if (is.null(median.ymin)) {
      median.ymin = min(data.col.median)
  }
  if (is.null(median.ymax)) {
      median.ymax = max(data.col.median)
  }

  ## Plot the profile
  par(mfrow=c(2,1))
  plot(time.points.rep.shift,
       data.col.median,
       main=paste(data.type, "median profiles"),
       col=data.colors,
       ylab="median per column",
       xlab="time (min)",
       ylim=c(median.ymin, median.ymax),
       panel.first=grid(col='#000000'),
       ...
       )
  
  ################################################################
  ## IQR profile

  ## Calculate IQR profiles
  (data.col.IQR <- apply(data.profiles, 2, IQR, na.rm=T))

  ## Calculate Y limits if not provided
  if (is.null(IQR.ymax)) {
      IQR.ymax = max(data.col.IQR)
  }

  ## plot the profile
  plot(time.points.rep.shift,
       data.col.IQR,
       main=paste(data.type, "IQR profiles"),
       col=data.colors,
       ylab="IQR per column",
       xlab="time (min)",
       ylim=c(IQR.ymin, IQR.ymax),
       panel.first=grid(col='#000000'),
       ...
       )
  par(mfrow=c(1,1))
}

## ##############################################################
## Permute elements of a data frame for permutation tests
permute.frame <- function(x, ## matrix or data frame
                          by="table" ## Permutation mode. Supported: table|row|col
                          ) {

  ## Cast the input data into a matrix (in case it would be a data.frame)
  x.perm <- as.matrix(x)
  
  if (by == "table") {
    ## Permute all cells of the original data table
    x.perm <- matrix(sample(as.vector(x.perm)),nrow=nrow(x),ncol=ncol(x))

  } else if (by == "row") {
    ## Row-wise permutation
    for (i in 1:nrow(x)) {
      x.perm[i,] <- sample(as.vector(as.matrix(x[i,])))
    }

  } else if (by == "col") {
    ## Column-wise permutation
    for (i in 1:ncol(x)) {
      x.perm[,i] <- sample(as.vector(as.matrix(x[,i])))
    }
    
  } else {
    stop("Invalid value for the parameter 'by'. Supported: table|row|col.")
  }

  return(x.perm)
}


################################################################
#
# Demo
#

chip.util.demo <- function(nrow = 250,
			   ncol = 7,
			   k = 11,
			   iter = 200,
                           cluster.col=rainbow(n=k)
			   ) {
    palette.rep <- rep(palette,length.out=k)
    x <- data.frame(matrix(rnorm(nrow*ncol),ncol=ncol))
    km <- kmean.profiles(x,k,iter,return="kmeans")
    plot.svd(x,col=cluster.col[km$cluster])
}


