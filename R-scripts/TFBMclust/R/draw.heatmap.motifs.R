####################################################
## Draw a heatmap usig a particular allgomeration
## method.
draw.heatmap.motifs <- function(dist.table, 
                                method = "average", 
                                clusters.list, 
                                metric = "Ncor", 
                                tree.pos = "column",
                                color.palette = "YlOrRd",
                                color.classes = 9
                                ){

#   ################################################################
#   ## Blue -> White -> Red palette
#   blue.white.red <- function() {
#     blue.levels <- c(rep(127,127), 126:0)/128
#     red.levels <- c(0:127,rep(127,127))/128
#     green.levels <- c(0:127, 126:0)/128
#     palette <- rgb(red.levels, green.levels, blue.levels)
#     return(palette)
#   }
#

library("RColorBrewer")
palette <- colorRampPalette(rev(brewer.pal(color.classes, color.palette)), space="Lab")

  metric.definition <- NULL
  ## If required, convert similarities to distances
  ## Similarity sores bounded to 1
  if ((metric == "Ncor")
      || (metric=="NcorS")
      || (metric=="cor")
      || (metric=="logocor")
      || (metric=="Nlogocor")
      || (metric=="Icor")
      || (metric=="NIcor")
      || (metric=="logoDP")
      || (metric=="mean_zscore")
  ) {
    ## cor       Pearson correlation (computed on residue occurrences in aligned columns)
    ## Ncor   		Relative width-normalized Pearson correlation
    ## logocor 			correlation computed on sequence logos
    ## Nlogocor 			Relative width-normalized logocor
    ## Icor 			Pearson correlation computed on Information content
    ## NIcor 			Relative width-normalized Icor
    metric.definition <- "correlation"

  } else if ((metric == "logoDP")
             || (metric == "cov")) {
    ## logoDP 			dot product of sequence logos
    ## cov 			covariance between residues in aligned columns

    stop("logoDP and cov metrics are not supported yet")

  } else if ((metric == "dEucl")
             || (metric == "NdEucl")
             || (metric == "NsEucl")
             || (metric == "SSD")
             || (metric == "SW")
             || (metric == "NSW")
             || (metric == "rank_mean")
  ) {
    ## dEucl 			Euclidian distance between residue occurrences in aligned columns
    ## NdEucl 			Relative width-normalized dEucl
    ## NsEucl 			similarity derived from Relative width-normalized Euclidian distance
    ## SSD 			Sum of square deviations
    ## SW 			Sandelin-Wasserman
    ## NSW 			Relative width-normalized Sandelin-Wasserman

    metric.definition <- "distance"

  } else if (metric == "match_rank") {
    ## match_rank rank of current match among all sorted matches
    stop("match_rank metric is not supported yet")

  } else {
    stop(paste(metric, "is not a valid metric", sep="\t"))
  }


  ## The font size is adapted relative to the number of input motifs
  ## If there are less than 25, set the font size to 25
  nb.motifs <- length(dim(dist.table)[1])
  font.size <- nb.motifs
  if(nb.motifs < 25){
    font.size <- 25
  }

  dist.table <- as.matrix(dist.table)

  ## Set the suported colors
  nb.clusters <- length(clusters.list)
  color <- rainbow(nb.clusters)
  cluster.to.color <- list(nb.clusters)
  color.counter <- 0


  ## Fill a list where each element correspond to a cluster name
  ## and the value its corresponding color
  sapply(names(clusters.list), function(x){
    color.counter <<- color.counter + 1
    cluster.to.color[[x]] <<- color[color.counter]
  })

  ## Create a vector with the the corresponding color of each motif
  color.order <- vector()
  sapply(colnames(dist.table), function(y){
    match.exp <- paste("^", y, "$", sep = "")
    current.cluster <- names(unlist(sapply(clusters.list, function(x){
                                                            grep(match.exp, x, value = TRUE, perl = TRUE)
                                                          })))

    color.order <<- append(color.order, cluster.to.color[[current.cluster]])
  })


  ## Color gradient for the heatmap
#  grad <- colorRampPalette(c("blue", "black"))(n = 299)
  if(metric.definition == "distance"){
    grad <- palette
  }else if(metric.definition == "correlation"){
    grad <- palette
#    grad <- colorRampPalette(c("red", "white"))(n = 256)
  }


## Calculate the bottom border
par(oma=c(2,0.5,0.5,2), family="mono")

  # Get the aligned consensuses, which will be used as the Row names
  col.names <-sapply(colnames(dist.table), function(x){
    as.vector(get.name(x))
  })
  col.names <- as.vector(col.names)


  ## Draw the heatmap
  heatmap.2(dist.table,

          ## Display the tree symetrically
          trace = "none",
          symm=TRUE,

          ## Set the colors of columns, rows and cells
          ColSideColors = color.order,
          RowSideColors = color.order,
          col = grad,

          ## The order of the values is set according these dendrograms
          Rowv = as.dendrogram(tree),
          Colv =as.dendrogram(tree),

          ## To show only the dendrogram in row
          dendrogram = tree.pos,

          ## Set the col and row labels
          labRow = col.names,
          labCol = col.names,

          ## Set the font size
          cexRow = 11/font.size + 0.08,
          cexCol = 11/font.size + 0.08,

          ## Set the key with the values
          key = TRUE,
          keysize = 1.5,
          key.xlab = "Distance",
          key.ylab = "",
          key.title = "Color key",
          density.info = "none"
  )
}