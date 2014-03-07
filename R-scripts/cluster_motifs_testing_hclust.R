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


## dir.results <- "/home/jaimecm/Documents/TAGC/Clustering_test/Prueba_Jacques/results"
## file.prefix <- file.path(dir.results, "Testing_10_02_2014")
dir.results <- "/Users/jvanheld/test/motif_clustering/results/peakmo_clustering"
file.prefix <- file.path(dir.results, "peakmo_example")
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


######################################
### Testinf dynamicTreeCut library ###
######################################
#a <- cutreeDynamic(tree, minClusterSize = 2, method = "tree")
#cutreeDynamicTree(tree, maxTreeHeight = 1, deepSplit = TRUE, minModuleSize = 2)

#tree$labels <- as.vector(description.table$consensus)
tree$labels <- paste(as.vector(description.table$consensus), 1:length(description.table$consensus))

## Bottom-up traversal of the tree to orientate the logos
merge.level <- 1
for (merge.level in 1:nrow(tree$merge)) {
  child1 <- tree$merge[merge.level,1]
  child2 <- tree$merge[merge.level,2]

  ## Simplest case: merging between two leaves
  if ((child1 < 0) && (child2 < 0)) {

    ## Identify the two motifs
    n1 <- min(-child1,-child2) ## row number of the first motif in the description table
    n2 <- max(-child1,-child2) ## row number of the second motif in the description table

    id1 <- as.vector(description.table[n1,"id"]) ## Id of the first motif
    id2 <- as.vector(description.table[n2,"id"]) ## Id of the second motif


    compa.nb <- which(compare.matrices.table[,"id1"] == id1 & compare.matrices.table[,"id2"] == id2)
    
    ## Choose the relative orientation of the two motifs
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])
    consensus1 <- as.vector(description.table[n1,"consensus"]) ## Consensus of the first motif
    if (strand=="R") {
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
  }
}
  


###############################################################################
## Plot the tree horzontally
dev.new(width=10, height=7)
par(mar=c(3,2,1,20),family="mono")
plot(as.dendrogram(tree), horiz=TRUE)

#### HERE I AM ####

### Cluster ("kkcyTTTGTTATGCAAATGvarkc" "byATTGTcATGCAAATGcaaky")
clusters <- cut(as.dendrogram(tree), h=0.2)
plot(clusters$upper)
plot(clusters$lower[[9]])  
str(clusters$lower[[9]])
labels(clusters$lower[[9]])

### Cluster ("kkcyTTTGTTATGCAAATGvarkc" "byATTGTcATGCAAATGcaaky" "dbhYTbvTTATGCATAAbvARdvh" )
clusters <- cut(as.dendrogram(tree), h=0.4)
plot(clusters$upper)
plot(clusters$lower[[7]]) 
str(clusters$lower[[7]])
labels(clusters$lower[[7]])

### Cluster ("kkcyTTTGTTATGCAAATGvarkc" "byATTGTcATGCAAATGcaaky" "dbhYTbvTTATGCATAAbvARdvh" "TtTGCATgACAATrr" )
clusters <- cut(as.dendrogram(tree), h=0.63)
plot(clusters$upper)
plot(clusters$lower[[5]]) 
str(clusters$lower[[5]])
labels(clusters$lower[[5]])




### Cluster ("tsatATGCAAATgwry" "cyhcATTTGCATAACAAwrr")
clusters <- cut(as.dendrogram(tree), h=0.3)
plot(clusters$upper)
plot(clusters$lower[[5]]) 
str(clusters$lower[[5]])
labels(clusters$lower[[5]])

### Cluster ("tsatATGCAAATgwry" "cyhcATTTGCATAACAAwrr" "hrydcATTTGCATATGcAAATgwr")
clusters <- cut(as.dendrogram(tree), h=0.45)
plot(clusters$upper)
plot(clusters$lower[[4]]) 
str(clusters$lower[[4]])
labels(clusters$lower[[4]])




### Cluster ("hrydcATTTGCATATGcAAATgwr" "tsatATGCAAATgwry" "cyhcATTTGCATAACAAwrr" "TtTGCATgACAATrr" "dbhYTbvTTATGCATAAbvARdvh" "kkcyTTTGTTATGCAAATGvarkc" "byATTGTcATGCAAATGcaaky")
clusters <- cut(as.dendrogram(tree), h=0.75)
plot(clusters$upper)
plot(clusters$lower[[3]]) 
str(clusters$lower[[3]])
labels(clusters$lower[[3]])

### Cluster ("wtATGCTAATww" "ygsATATGCGCATATGCArATrwr" "hrydcATTTGCATATGcAAATgwr" "tsatATGCAAATgwry" "cyhcATTTGCATAACAAwrr" "TtTGCATgACAATrr" "dbhYTbvTTATGCATAAbvARdvh" "kkcyTTTGTTATGCAAATGvarkc" "byATTGTcATGCAAATGcaaky")
clusters <- cut(as.dendrogram(tree), h=0.90)
plot(clusters$upper)
plot(clusters$lower[[2]]) 
str(clusters$lower[[2]])
labels(clusters$lower[[2]])


##########################
### Find next neighbor ###
##########################



Output.path  <-  "Subtrees.PDF"
pdf(Output.path)


current.num.leaves  <-  NULL
prev.num.leaves  <-  NULL
current.branch.sizes  <-   NULL
prev.branch.sizes  <-  NULL
Clusters.found  <-  NULL


for (distance in seq(0,0.3, by= 0.02)){

  ### Obtain the clusters at determined height
  clusters <- cut(as.dendrogram(tree), h=distance)
  current.num.leaves  <-  length(clusters$lower)

  
  ### Initialize values of leaves and branch sizes
  if(distance == 0){

    ### Update the number of leaves
    current.num.leaves  <-  length(clusters$lower)
    prev.num.leaves  <-  current.num.leaves

    current.branch.sizes  <-  lapply(clusters$lower, function(x){length(labels(x))})
    names(current.branch.sizes)  <-  sapply(clusters$lower, function(X){ paste(sort(labels(X)), collapse = "_")})
    
    prev.branch.sizes  <-  current.branch.sizes
  }

  
  ### Check if finds a cluster
  if(current.num.leaves  != prev.num.leaves){
    
    ### Calculate the size of each branch
    current.branch.sizes  <-  lapply(clusters$lower, function(x){length(labels(x))})
    names(current.branch.sizes)  <-  sapply(clusters$lower, function(X){ paste(sort(labels(X)), collapse = "_")})
    
    print(names(current.branch.sizes))
    print(names(prev.branch.sizes))
    
    New.cluster.found  <-  Get.New.Cluster.Number(current.branch.sizes, prev.branch.sizes)
    
    prev.branch.sizes  <-  current.branch.sizes

    Clusters.found  <-  append(Clusters.found, labels(clusters$lower[[New.cluster.found]]))
    plot(clusters$lower[[New.cluster.found]])
  }

  prev.num.leaves  <-  current.num.leaves
  
}
dev.off()


###################################################
###################################################
Get.New.Cluster.Number  <-  function(curr, prev){
  
  cluster.index  <-  NULL
  cluster.index  <-  append(cluster.index, setdiff(names(current.branch.sizes), names(prev.branch.sizes)))
  return(cluster.index)
}


##########################################
### In progress 
Update.Vector  <-  function(curr, prev){

  Diff.pos  <-  Get.New.Cluster.Number(curr,prev)

  ### At the last position
  if(Diff.pos == length(prev) - 1){
    prev  <-  prev[c([1:Diff.pos])]
    return()
  }

  ### At the first position
  else if(Diff.pos == 1){
    prev  <-  prev[c([2:length(curr)])]
    return()
  }

  ### At middle position
  else{
    
  }
  
}


clusters <- cut(as.dendrogram(tree), h=0.20)
starting.nodes  <-  length(clusters$lower)
plot(clusters$upper)
plot(clusters$lower[[9]]) 
str(clusters$lower[[9]])
labels(clusters$lower[[2]])




########################################

