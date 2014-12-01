##################################################################
## Build the tree by hierarchical clustering using Rcpp library
## if it is indicated, export the tree in Newick format
hclust.motifs <- function(dist.matrix, hclust.method = "average"){

  ## Require Rclusterpp if it is required
  if(!require("Rclusterpp")){
    install.packages("Rclusterpp")
  }

  ## Load library
  suppressPackageStartupMessages(library(Rclusterpp, warn.conflicts=FALSE))

  return(Rclusterpp.hclust(dist.matrix, method = hclust.method))
}