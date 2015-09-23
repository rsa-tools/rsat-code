####################################################################
## Following the order on the tree object, associate one
## number (1-4) corresponding to one color (black,red,blue,green)
## All the members of one cluster will be colored with the same
## color in the tree.
color.code.clusters <- function(clusters.ids, tree){

  tree.order <- get.id(tree$order)

  ## Following the order on the tree object, identify to which
  ## cluster corresponds each branch
  clusters.names.matrix <- sapply(clusters.ids, function(X){
    tree.order %in% X
  })
  clusters.names.matrix <- t(clusters.names.matrix)

  cluster.names.order <- apply(clusters.names.matrix, 2, function(X){
    names(which(X > 0))
  })

  ## Following the order on the tree object, associate one
  ## number (1-4) corresponding to one color (black,red,blue,green)
  ## All the members of one cluster will be colored with the same
  ## color in the tree.
  cluster.names.order.unique <- unique(cluster.names.order)
  color.code <- NULL
  color.counter <- 0
  cluster.counter <- 0
  sapply(1:length(cluster.names.order.unique), function(x){

    color.counter <<- color.counter + 1
    labels.rep <- length(which(cluster.names.order == cluster.names.order.unique[x]))
    color.code <<- append(color.code, rep(color.counter, times = labels.rep))

    ## Reset the color counter to include only black, red, green and blue colors
    ## in the dendrogram's branches
    if(color.counter == 4){
      color.counter <<- 0
    }
  })
  return(color.code)
}