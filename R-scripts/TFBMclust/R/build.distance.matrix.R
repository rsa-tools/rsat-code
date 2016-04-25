###########################################################
## Build a distance matrix from the distance metric list
build.distance.matrix <- function(metric="Ncor"){

  dist.table <- NULL
  distances.objects <- list()

  ## Extract metric values
  metric.values <- global.compare.matrices.table[,metric]

  ## If required, convert similarities to distances
  ## Similarity sores bounded to 1
  if ((metric == "Ncor")
      || (metric=="NcorS")
      || (metric=="cor")
      || (metric=="logocor")
      || (metric=="Nlogocor")
      || (metric=="Icor")
      || (metric=="NIcor")
  ) {
    ## cor       Pearson correlation (computed on residue occurrences in aligned columns)
    ## Ncor 			Relative width-normalized Pearson correlation
    ## logocor 			correlation computed on sequence logos
    ## Nlogocor 			Relative width-normalized logocor
    ## Icor 			Pearson correlation computed on Information content
    ## NIcor 			Relative width-normalized Icor
    metric.dist <- 1 - metric.values

  } else if (metric == "cov") {

    ## cov 			covariance between residues in aligned columns

    stop("cov metric is not supported yet")

  } else if ((metric == "dEucl")
             || (metric == "NdEucl")
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
    ## logoDP   		dot product of sequence logos

    metric.dist <- metric.values

  } else if ( (metric == "mean_zscore")
              || (metric == "logoDP")){

    metric.dist <- max(metric.values) - metric.values

  }else {
    stop(paste(metric, "is not a valid metric", sep="\t"))
  }


  ## Add a column with metric column to the compare matrices table, will
  ## be used to generate a cross-table
  global.compare.matrices.table$metric.dist <- metric.dist

  ## Build the distance table from the column metric
  dist.table <- xtabs(metric.dist ~ id1+id2, global.compare.matrices.table)

  ## Ensure that symmetrical distances are defined
  dist.table.sym <- pmax(dist.table,t(dist.table))

  ## Cast the distance table into an object of class "dist"
  dist.matrix <- as.dist(dist.table.sym)

  return(list(table = dist.table.sym, matrix = dist.matrix))
}

