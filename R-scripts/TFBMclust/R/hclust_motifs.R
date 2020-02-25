###################################################################
## Compute the hierarchical clustering on a given distance matrix
hclust.motifs <- function(dist.matrix, hclust.method = "average"){

   return(hclust(dist.matrix, method = hclust.method))
}