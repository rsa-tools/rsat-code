###################################################################
###################################################################
## Library with R functions required in R-scripts/cluster_motifs.R


#######################################################
#######################################################
## Etxract the tree from an hclust object 
createLeafNode <- function(hclust, i) {
  list(label = hclust$labels[[i]],
       order = hclust$order[[i]])
}

########################################################
########################################################
## Convert an hclust tree into a JSON format tree
hclustToTree <- function(hclust) {
  if (length(hclust$merge) == 0)
    return(NULL)
  merges <- list()
  for (index in 1:nrow(hclust$merge)) {
    left <- hclust$merge[index, 1]
    right <- hclust$merge[index, 2]
    if (left < 0)
      left <- createLeafNode(hclust, -left)
    else
      left <- merges[[left]]
    if (right < 0)
      right <- createLeafNode(hclust, -right)
    else
      right <- merges[[right]]
    if (left$order > right$order) {
      tmp <- left
      left <- right
      right <- tmp
    }
    merges[[index]] <- list(
                            children = list(
                              left,
                              right
                              ),
                            order = left$order
                            )
  }
  return(merges[nrow(hclust$merge)])
}


################################################################
################################################################
## Collect the list of leaves (motif) associated to each
## internal node (motif cluster) of the tree.
##
## This method works for any hclust result tree, irrespective of the
## fact that it refers to motifs or anything else.
##
## Usage: generate a vector with the list of leaves (string-formatted)
## per internal node
##   treenodes <- leaves.per.node(tree) ## The input tre must be an hclust result
##
## Then to get a vector with the leaves associated to a given internal
## node (e.g. node 4):
##   leaves.for.node.4 <- as.numeric(unlist(strsplit(tree.nodes[4], split=" ")))
##
leaves.per.node <- function (tree 
                             ) {
  merge.table <- tree$merge
  leave.lists <- list()

  for (i in 1:nrow(merge.table)) {
    branch1 <- merge.table[i, 1]
    branch2 <- merge.table[i, 2]

    ## Depending on whether the left branch points to a leave or an
    ## iternal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch1 < 0) {
      nodes1 <- -branch1 ## branch one only contains one leave
    } else {
      nodes1 <- leave.lists[branch1]
    }

    ## Dependingon whether the right branch points to a leave or an
    ## iternal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch2 < 0) {
      nodes2 <- -branch2 ## branch one only contains one leave
    } else {
      nodes2 <- leave.lists[branch2]
    }
   
    leave.lists[i] <- paste(nodes1, nodes2)
  }
  leave.lists <- sapply(leave.lists, function(x){ as.numeric(unlist(strsplit(x, " "))) })
  return (leave.lists)
}


################################################################
################################################################
## Identify the "central" leaf (motif) of a subtre (cluster),
## i.e. the motif with the smallest average distance to all
## other motifs of this cluster.
central.motifs.ids <- function(ids1, ids2){
  
  ## Get the comparison numbers
  compa.numbers <- sapply(ids2, function(X){
    sapply(ids1, function(Y){
      get.compa.nb(X, Y)
    })
  })
  compa.numbers <- as.vector(compa.numbers)

  ## Get the ids of the less distant nodes
  compa.info <- compare.matrices.table[compa.numbers,][which(compare.matrices.table[compa.numbers,score] == max(compare.matrices.table[compa.numbers,score])),c("id1","id2")]
  return(as.vector(unlist(compa.info)))
}


########################################################
########################################################
## Get the number of comparison between the input IDs
## in the compare-matrices results table
get.compa.nb <- function(ID1,ID2){
  
  return(which( (compare.matrices.table[,"id2"] == ID2 & compare.matrices.table[,"id1"] == ID1) | (compare.matrices.table[,"id1"] == ID2 & compare.matrices.table[,"id2"] == ID1) ))
}


########################################################
########################################################
## Given the number of leaf, get the id of the motif
get.id <- function(num){
  
  return(as.vector(sapply(num, function(X){
    description.table[X,"id"]
  })))
}


########################################################
########################################################
## Given a vector with IDs, return a list with the 
## information (consensus, number, id, strand, spacer)
## of the inverted alignment
inverted.aligment <- function(ids){

  temporal.list <- list()
  
  ## Get the higher length among the consensuses
  consensuses <- sapply(motifs.info[ids], function(X){
    nchar(X[["consensus"]])
  })
  max.cons.length <- max(consensuses)

  ## Create the new consensuses with the inverted orientation
  temporal.list <- sapply(ids, function(X){
    
    inverted.aligment.list <- list()
    inverted.consensus <- NULL
    spacer.length <- NULL
    spacer <- NULL

    ## Create the spacer
    spacer.length <- as.numeric(motifs.info[[X]][["spacer"]])
    spacer <- paste(rep("-", times = spacer.length), collapse = "")

    ## Get the consensus in the inverted orientation
    if(motifs.info[[X]][["strand"]] == "D"){
      inverted.consensus <- as.vector(description.table[[as.numeric(motifs.info[[X]][["number"]]), "rc_consensus"]])
      inverted.aligment.list[[X]][["strand"]] <- "R"
    }else{
      inverted.consensus <- as.vector(description.table[[as.numeric(motifs.info[[X]][["number"]]), "consensus"]])
      inverted.aligment.list[[X]][["strand"]] <- "D"
    }

    ## Fill the spaces with "-"
    new.consensus <- paste(inverted.consensus, spacer, sep = "")
    new.consensus <- paste(paste(rep("-", times = max.cons.length - nchar(new.consensus)), collapse = ""), new.consensus, sep = "")
    new.consensus <- paste(unlist(strsplit(new.consensus, "[A-Za-z]+", perl = TRUE))[1], inverted.consensus, sep = "")

    ## Fill the new list
    inverted.aligment.list[[X]][["consensus"]] <- inverted.consensus
    inverted.aligment.list[[X]][["number"]] <- as.numeric(motifs.info[[X]][["number"]])
    inverted.aligment.list[[X]][["spacer"]] <- length(unlist(strsplit(new.consensus, "-")))-1
    

    return(inverted.aligment.list)
           
  })

  names(temporal.list) <- sapply(strsplit(names(temporal.list), "\\."), function(X){return(X[1])})
  return(temporal.list)
}

