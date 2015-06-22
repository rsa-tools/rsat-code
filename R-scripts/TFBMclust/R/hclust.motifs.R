##################################################################
## Build the tree by hierarchical clustering using Rcpp library
## if it is indicated, export the tree in Newick format
hclust.motifs <- function(dist.matrix, hclust.method = "average"){

#   if(!require("Rclusterpp")){
#     install.packages("Rclusterpp")
#   }

#    return(Rclusterpp.hclust(dist.matrix, method = hclust.method))
   return(hclust(dist.matrix, method = hclust.method))
}