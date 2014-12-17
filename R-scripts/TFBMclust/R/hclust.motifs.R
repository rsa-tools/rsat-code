##################################################################
## Build the tree by hierarchical clustering using Rcpp library
## if it is indicated, export the tree in Newick format
hclust.motifs <- function(dist.matrix, hclust.method = "average"){
  return(Rclusterpp.hclust(dist.matrix, method = hclust.method))
}