### Install rjson package
# install.packages("RJSONIO")
library("RJSONIO")
library("ctc")
library("reshape")
# library("dynamicTreeCut")
library("dendroextras")

## Define parameters
score <- "Ncor";
hclust.method <- "average"

#################################################################################
## Etxract the tree from an hclust object and convert it into json object.
createLeafNode <- function(hclust, i) {
  list(label = hclust$labels[[i]],
       order = hclust$order[[i]])
}

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


dir.results <- "/home/jaimecm/Documents/TAGC/Clustering_test/Prueba_Jacques/results"
file.prefix <- file.path(dir.results, "Testing_10_03_2014")

#dir.results <- "/Users/jvanheld/test/motif_clustering/results/peakmo_clustering"
#file.prefix <- file.path(dir.results, "peakmo_example")

setwd(dir.results)

infile <- paste(sep="", file.prefix, "_pairwise_compa.tab")
description.file <-  paste(sep="", file.prefix, "_pairwise_compa_matrix_descriptions.tab")

##################################
## Read matrix comparison table
compare.matrices.table <- read.csv(infile, sep = "\t", comment.char = ";")
names(compare.matrices.table)[1] <- sub("^X.", "", names(compare.matrices.table)[1])


################################################################
## Read description table 
description.table <- read.csv(description.file, sep = "\t", comment.char = ";")
names(description.table)[1] <- sub("^X.", "", names(description.table)[1])
matrix.labels <-  as.vector(description.table$label)
names(matrix.labels) <- as.vector(description.table$id)

## Extract score values
score.values <- compare.matrices.table[,score] 
score.dist <- 1 - score.values

## Add a column with score column to the compare matrices table, will
## be used to generate a cross-table
compare.matrices.table$score.dist <- score.dist

################################################################
## Build a distance matrix from the distance score list
dist.table <- t(xtabs(score.dist ~ name1+name2, compare.matrices.table) )
## Ensure that symmetrical distances are defined
for (i in 1:nrow(dist.table)) {
  for (j in i:ncol(dist.table)) {
    if (i==j) {next}
    dist.max <- max(dist.table[i,j], dist.table[j,i])
    dist.table[i,j] <- dist.max
    dist.table[j,i] <- dist.max
  }
}
print(dist.table)
#write.table(dist.table, file = distance.table, quote = FALSE, row.names = TRUE, col.names=NA, sep = "\t")

## Convert distance table into a distance matrix, required by hclust()
dist.matrix <- as.dist(dist.table)


##############################################
### Runs and plot the hierarchical cluster
tree <- hclust(dist.matrix, method = hclust.method)
tree$labels <- as.vector(description.table$label)
write(hc2Newick(tree, flat = TRUE), file="test_output_hclust_average.newick")


###########################
### Traversing the tree ###
###########################

##########################################
## Given a merge level return all the
## motifs clustered in that level
get.all.leaves.below <- function(num){

  X <- tree$merge[num,1]
  Y <- tree$merge[num,2]

  if(X < 0 & Y < 0){
    return(abs(c(X,Y)))
  }
  
  else if( (X < 0 & Y > 0) | (X > 0 & Y < 0) ){
    return(abs(c(min(c(X,Y)), get.all.leaves.below(max(c(X,Y))))))
  }

  else if (X > 0 & Y > 0){
    return(abs(c(get.all.leaves.below(X), get.all.leaves.below(Y))))
  }
}


########################################
## Given a cluster number returns all
## nodes within that cluster
#get.all.leaves.in.cluster <- function(cluster){
  
#}


#################################################################################
## Get the if of the motif with the leatest distance between the motif and the
## leaves in the cluster
minor.distance.id <- function(id,motif.numbers){
  
  id.2 <- get.id(motif.numbers)
  
  n.2 <- as.vector(sapply(id.2, function(X){
    get.compa.nb(id,X)
  }))

  compa.info <- compare.matrices.table[n.2,][which(compare.matrices.table[n.2,"Ncor"] == max(compare.matrices.table[n.2,"Ncor"])),c("id1","id2")]
  id.highest <- as.vector(compa.info[,which(compa.info != id)])
  return(id.highest)
}


########################################################
## Get the number of comparison between the input IDs
## in the compare-matrices results table
get.compa.nb <- function(ID1,ID2){
  return(which( (compare.matrices.table[,"id2"] == ID2 & compare.matrices.table[,"id1"] == ID1) | (compare.matrices.table[,"id1"] == ID2 & compare.matrices.table[,"id2"] == ID1) ))
}

#################################################
## Given the IDs of the leaves in two clusters
## return the pairs of closest IDs 
get.ids.minor.distance.two.clusters <- function(cluster1, cluster2){

}


########################################################
## Given the number of leaf, get the id of the motif
get.id <- function(num){
  return(as.vector(sapply(num, function(X){
    description.table[X,"id"]
  })))
}


#tree$labels <- as.vector(description.table$consensus)
tree$labels <- paste(as.vector(description.table$consensus), 1:length(description.table$consensus))

#############################################################
## Bottom-up traversal of the tree to orientate the logos
merge.level <- 1
merge.information.list <- list()
motif.consensus <- list()
motif.number <- list()

#for (merge.level in 1:nrow(tree$merge)) {
for (merge.level in 1:6) { 
  child1 <- tree$merge[merge.level,1]
  child2 <- tree$merge[merge.level,2]
  
  merge.id <- paste("merge_level_", merge.level, sep = "")
  merge.information.list[[merge.id]] <- list()
  
  
###############################################
  ## Simplest case: merging between two leaves ##
###############################################
  if ((child1 < 0) && (child2 < 0)) {
    
    ## Identify the two motifs
    n1 <- min(-child1,-child2) ## row number of the first motif in the description table
    n2 <- max(-child1,-child2) ## row number of the second motif in the description table
    
    id1 <- get.id(n1) ## Id of the first motif
    id2 <- get.id(n2) ## Id of the second motif
    
    ## Comparison number in the compare-matrices table
    compa.nb <- get.compa.nb(id1,id2)
    
    ## Choose the relative orientation of the two motifs
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])
    consensus1 <- as.vector(description.table[n1,"consensus"]) ## Consensus of the first motif
    if (strand == "R") {
      consensus2 <- as.vector(description.table[n2,"rc_consensus"]) ## Consensus of the second motif
    } else {
      consensus2 <- as.vector(description.table[n2,"consensus"]) ## Consensus of the second motif
    }
    
    ## Add the offset to the logos
    offset <- as.vector(compare.matrices.table[compa.nb, "offset"])
    spacer <- paste(collapse="",rep(x="-",times=abs(offset)))
    
    if (offset < 0) {
      consensus1 <- paste(sep="", spacer, consensus1)
    } else {
      consensus2 <- paste(sep="", spacer, consensus2)
    }
    
    ## Reset the consensus with the aligned and re-oriented consensus
    tree$labels[n1] <- paste(consensus1, n1)      
    tree$labels[n2] <- paste(consensus2, n2)
    
    ## Store the information in a list of lists
    merge.information.list[[merge.id]][[paste("motif", n1, sep = "_")]][["id"]] <- id1
    merge.information.list[[merge.id]][[paste("motif", n1, sep = "_")]][["consensus"]] <- consensus1
    
    merge.information.list[[merge.id]][[paste("motif", n2, sep = "_")]][["id"]] <- id2
    merge.information.list[[merge.id]][[paste("motif", n1, sep = "_")]][["consensus"]] <- consensus2
    
    merge.information.list[[merge.id]][[paste("motifs", n1, n2, sep = "_")]][["relative_strand"]] <- strand
    
    
    motif.consensus[[id1]] <- consensus1
    motif.consensus[[id2]] <- consensus2
    
    motif.number[[id1]] <- n1
    motif.number[[id2]] <- n2
  }
  
  
  ############################################
  ## Merging aligned logos with another one ##
  ############################################
  if (((child1 < 0) && (child2 > 0))
      || ((child1 > 0) && (child2 < 0)) ) {
    
    ## Identify the two motifs (n1 -> leaf) (N2 -> cluster)
    n1 <- abs(min(child1, child2)) 
    N2 <- max(child1, child2) 
    
    ## Id of the leaf
    id1 <- get.id(n1)
    
    
    ## Get all motifs in the cluster
    N2.motifs <- get.all.leaves.below(N2)
    
    ## Get the motif with the lowest distance between the leaf and the leaves in the cluster
    id2 <- minor.distance.id(id1,N2.motifs)
    compa.nb <- get.compa.nb(id1, id2)
    
    ## Get the orientation of the previous cluster
    prev.cluster.strand <- merge.information.list[[paste("merge_level_", N2, sep = "")]][[3]][["relative_strand"]]
    
    ## Choose the relative orientation of the motif
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])
    if (prev.cluster.strand == "R") {
      consensus1 <- as.vector(description.table[n1,"rc_consensus"]) 
    } else {
      consensus1 <- as.vector(description.table[n1,"consensus"])
    }
    
    motif.consensus[[id1]] <- consensus1
    motif.number[[id1]] <- n1
    tree$labels[n1] <- paste(consensus1, n1)   
    
    ## Add the offset to the logos
    if(prev.cluster.strand == strand){
      offset <- as.vector(compare.matrices.table[compa.nb, "offset"])
    }
    else{
      offset <- abs(as.vector(compare.matrices.table[compa.nb, "offset"]))
    }
    cluster.spacer <- length(unlist(strsplit(motif.consensus[[id2]], "-")))-1
    spacer <- paste(collapse="",rep(x="-",times = abs(offset) - cluster.spacer))
    
    if (offset > 0) {
      
      consensus1 <- paste(sep="", spacer, consensus1)
      tree$labels[n1] <- paste(consensus1, n1)
    } else {
      
      ## Get the ids of the leaves in the cluster
      ids.cluster <- sapply(get.all.leaves.below(N2), get.id)
      for (id in ids.cluster){
        motif.consensus[[id]] <- paste(paste(spacer, motif.consensus[[id]], sep=""), motif.number[[id]]) 
        tree$labels[motif.number[[id]]] <- motif.consensus[[id]] 
      }
    }

    merge.information.list[[merge.id]][[paste("motif", n1, sep = "_")]][["id"]] <- id1
    merge.information.list[[merge.id]][[paste("motif", n1, sep = "_")]][["consensus"]] <- consensus1
    
    merge.information.list[[merge.id]][[paste("cluster", N2, sep = "_")]][["id"]] <- N2
    
    merge.information.list[[merge.id]][[paste("motifs", paste(get.all.leaves.below(N2), collapse = "_"), sep = "_")]][["relative_strand"]] <- strand
    
  } 
} 

###############################################################################
## Plot the tree horzontally
dev.new(width=10, height=7)
par(mar=c(3,2,1,20),family="mono")
plot(as.dendrogram(tree), horiz=TRUE)

