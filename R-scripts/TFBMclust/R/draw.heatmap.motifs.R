####################################################
## Draw a heatmap usig a particular allgomeration
## method.
draw.heatmap.motifs <- function(dist.table, method = "average", clusters.list, alignment.list, score = "Ncor"){

  ################################################################
  ## Blue -> White -> Red palette
  blue.white.red <- function() {
    blue.levels <- c(rep(127,127), 126:0)/128
    red.levels <- c(0:127,rep(127,127))/128
    green.levels <- c(0:127, 126:0)/128
    palette <- rgb(red.levels, green.levels, blue.levels)
    return(palette)
  }

  metric.definition <- NULL
  ## If required, convert similarities to distances
  ## Similarity sores bounded to 1
  if ((score == "Ncor")
      || (score=="cor")
      || (score=="logocor")
      || (score=="Nlocogor")
      || (score=="Icor")
      || (score=="NIcor")
  ) {
    ## cor       Pearson correlation (computed on residue occurrences in aligned columns)
    ## Ncor   		Relative width-normalized Pearson correlation
    ## logocor 			correlation computed on sequence logos
    ## Nlogocor 			Relative width-normalized logocor
    ## Icor 			Pearson correlation computed on Information content
    ## NIcor 			Relative width-normalized Icor
    metric.definition <- "correlation"

  } else if ((score == "logoDP")
             || (score == "cov")) {
    ## logoDP 			dot product of sequence logos
    ## cov 			covariance between residues in aligned columns

    stop("logoDP and cov scores are not supported yet")

  } else if ((score == "dEucl")
             || (score == "NdEucl")
             || (score == "NdEucl")
             || (score == "NsEucl")
             || (score == "SSD")
             || (score == "SW")
             || (score == "NSW")
  ) {
    ## dEucl 			Euclidian distance between residue occurrences in aligned columns
    ## NdEucl 			Relative width-normalized dEucl
    ## NsEucl 			similarity derived from Relative width-normalized Euclidian distance
    ## SSD 			Sum of square deviations
    ## SW 			Sandelin-Wasserman
    ## NSW 			Relative width-normalized Sandelin-Wasserman

    metric.definition <- "distance"

  } else if (score == "match_rank") {
    ## match_rank rank of current match among all sorted matches
    stop("match_rank score is not supported yet")

  } else {
    stop(paste(score, "is not a valid score", sep="\t"))
  }


  ## The font size is adapted relative to the number of input motifs
  ## If there are less than 25, set the font size to 25
  nb.motifs <- length(alignment.list)
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
  if(metric.definition == "correlation"){
    grad <- blue.white.red()
  }else if(metric.definition == "correlation"){
    grad <- colorRampPalette(c("red", "white"))(n = 256)
  }


  ## Calculate the bottom border
  rigth <- round(170/font.size + (font.size/2 * 0.001), digits = 2)
  bottom <- round(190/font.size + (font.size/2 * 0.001), digits = 2)
#    rigth <- round(136, digits = 2)
#    bottom <- round(136, digits = 2)

  par(oma=c(bottom,0.5,0.5,rigth), family="mono")

  # Get the aligned consensuses, which will be used as the Row names
  consensus <-sapply(colnames(dist.table), function(x){
    as.vector(alignment.list[[x]][["consensus_d"]])
  })
  consensus <- as.vector(consensus)

  # Get the aligned consensuses, which will be used as the Row names
  columns.heatmap <- sapply(colnames(dist.table), function(x){
    if(x == alignment.list[[x]][["name"]]){
      as.vector(x)
    } else{
      paste(x, alignment.list[[x]][["name"]], sep = "  ")
    }
  })
  columns.heatmap <- as.vector(columns.heatmap)


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
          dendrogram = "column",

          ## Set the col and row labels
          labRow = consensus,
          labCol = columns.heatmap,

          ## Set the margins
          cexRow = 20/font.size + 0.1,
          cexCol = 20/font.size + 0.02,

          ## Set the key with the values
          key = TRUE,
          keysize = 1.5,
          key.xlab = "Distance",
          key.ylab = "",
          key.title = "Color key",
          density.info = "none"
  )
}



