################################################################
## Collect the list of leaves (motif) associated to each
## internal node (motif cluster) of the tree.
##
## This method works for any hclust result tree, irrespective of the
## fact that it refers to motifs or anything else.
##
## Usage: generate a vector with the list of leaves (string-formatted)
## per internal node
##   treenodes <- leaves.per.node(tree)
#
## The input tree must be an hclust result
##
## Then to get a vector with the leaves associated to a given internal
## node (e.g. node 4):
##   leaves.for.node.4 <- as.numeric(unlist(strsplit(tree.nodes[4], split=" ")))
##
leaves.per.node <- function (tree) {
  merge.table <- tree$merge
  leave.lists <- list()

  for (i in 1:nrow(merge.table)) {
    branch1 <- merge.table[i, 1]
    branch2 <- merge.table[i, 2]

    ## Depending on whether the left branch points to a leave or an
    ## internal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch1 < 0) {
      nodes1 <- -branch1 ## branch one only contains one leave
    } else {
      nodes1 <- leave.lists[branch1]
    }

    ## Depending on whether the right branch points to a leave or an
    ## internal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch2 < 0) {
      nodes2 <- -branch2 ## branch two only contains one leave
    } else {
      nodes2 <- leave.lists[branch2]
    }

    leave.lists[i] <- paste(nodes1, nodes2)
  }

  ## Transform the list to export it
  if(length(leave.lists) > 1){
    leave.lists <- sapply(leave.lists, function(x){ as.numeric(unlist(strsplit(x, " "))) })
  } else{
    leave.lists <- list(as.numeric(unlist(strsplit(leave.lists[[1]], " "))))
  }
  return (leave.lists)
}
