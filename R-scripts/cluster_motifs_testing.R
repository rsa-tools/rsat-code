################################################################
## Load paameters for the demo 1: motifs discovered with peak-motif demo.
demo.nb <- 1

## How to chance this path to make it general? 
rsat.dir <- Sys.getenv("RSAT")
source(file.path(rsat.dir, "R-scripts", "cluster_motifs_lib.R"))
dir.demo <- file.path(rsat.dir, "public_html", "demo_files")

demo.prefix <- "matrix-clustering_demo_peak-motifs"
file.prefix.peakmo <- file.path(dir.demo, demo.prefix)

if (demo.nb == 1) {
  file.prefix <- file.prefix.peakmo
}



## Define parameters
score <- "Ncor";
hclust.method <- "average"

## RDB matrices
#dir.results <- "/home/jaimecm/Documents/TAGC/Clustering_test/Test_different_hclust_methods/results/RDB_TF_Fam"
#file.prefix <- file.path(dir.results, "E_coli_LacI_LysR_")

infile <- paste(sep="", file.prefix, "_pairwise_compa.tab")
description.file <-  paste(sep="", file.prefix, "_matrix_descriptions.tab")

## Create a temporary result dir if required
dir.results <- file.path(Sys.getenv("HOME"), "rsat_demo", demo.prefix)
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)
setwd(dir.results)
out.prefix <- file.path(dir.results, "clustered_motifs")

stop("OK")

## Jaime: please check if all the subsequent code can be suppressed 
## (already incorporated in cluster_motifs.R)


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

## new.labels <- sapply(as.vector(description.table$label), function(X){

##   unlist(strsplit(X,"_"))[[length(unlist(strsplit(X,"_")))]]
## })
## names(new.labels) <- NULL

## tree$labels <- new.labels

#write(hc2Newick(tree, flat = TRUE), file="test_output_hclust_average.newick")




###########################
### Traversing the tree ###
###########################

#tree$labels <- as.vector(description.table$consensus)
tree$labels <- paste(as.vector(description.table$consensus), 1:length(description.table$consensus))

#############################################################
## Bottom-up traversal of the tree to orientate the logos
merge.level <- 1
motifs.info <- list()
#for (merge.level in 1:nrow(tree$merge)) {

merge.levels.leaves <- leaves.per.node(tree)
for (merge.level in 1:9) {

  child1 <- tree$merge[merge.level,1]
  child2 <- tree$merge[merge.level,2]

  ########################################
  ## Case 1: merging between two leaves ##
  ########################################
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
    motifs.info[[id1]][["strand"]] <- "D" 
    if (strand == "R") {
      consensus2 <- as.vector(description.table[n2,"rc_consensus"]) ## Consensus of the second motif
      motifs.info[[id2]][["strand"]] <- "R"
    } else {
      consensus2 <- as.vector(description.table[n2,"consensus"]) ## Consensus of the second motif
      motifs.info[[id2]][["strand"]] <- "D"
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
    
    ## Store the strand of the cluster, the consensuses and numbers
    motifs.info[[id1]][["consensus"]] <- consensus1
    motifs.info[[id1]][["number"]] <- n1
    
    motifs.info[[id2]][["consensus"]] <- consensus2
    motifs.info[[id2]][["number"]] <- n2

    motifs.info[[id1]][["spacer"]] <- length(unlist(strsplit(motifs.info[[id1]][["consensus"]], "-")))-1
    motifs.info[[id2]][["spacer"]] <- length(unlist(strsplit(motifs.info[[id2]][["consensus"]], "-")))-1
  }


  ############################################
  ## Case 2: merging a motif with a cluster ##
  ############################################
  if(((child1 < 0) && (child2 > 0)) || ((child1 > 0) && (child2 < 0))){
    
    n1 <- abs(min(child1, child2))
    n2 <- merge.levels.leaves[[merge.level]][which(merge.levels.leaves[[merge.level]] != n1)]
    N2 <- abs(max(child1, child2))
    
    ## Get ids
    motifs.info[[get.id(n1)]][["number"]] <- as.numeric(n1)
    motifs.info[[get.id(n1)]][["spacer"]] <- 0
    ids2 <- get.id(n2)
    
    ## Find the central motif of the cluster
    central.motifs <- central.motifs.ids(get.id(n1), ids2)
    id1 <- central.motifs[1]
    id2 <- central.motifs[2]
    
    ## Get the comparison number in the compare-matrices table
    compa.nb <- get.compa.nb(id1,id2)
    
    ## Get the strand
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])
    
    ## Identified the new and the aligned motif
    switch.ids <- 0
    if(length(motifs.info[[id1]]) > 2){
      aligned <- id1
      new <- id2
      temporal <- n1
      n1 <- n2
      n2 <- temporal
      switch.ids <- 1
    } else{
      aligned <- id2
      new <- id1
    }

    ## Get the offset
    offset <- as.vector(compare.matrices.table[compa.nb, "offset"])

    ## Assign values for the cases
    case <- 0
    if(switch.ids == 1){
      
      if(strand == "D"){
        
        if(motifs.info[[aligned]][["strand"]] == "D"){
          case <- 1
        } else{
          case <- 2
        }
      } else{
       if(motifs.info[[aligned]][["strand"]] == "D"){
          case <- 3
        } else{
          case <- 4
        } 
      }
    }else{
      if(strand == "D"){
        if(motifs.info[[aligned]][["strand"]] == "D"){
          case <- 5
        } else{
          case <- 6
        }
      } else{   
        if(motifs.info[[aligned]][["strand"]] == "D"){
          case <- 7
        } else{
          case <- 8
        } 
      } 
    }

    ## Get the spacers
    aligned.spacer <- as.numeric(motifs.info[[aligned]][["spacer"]])
    new.spacer <- as.numeric(motifs.info[[new]][["spacer"]])

    ## Choose the consensus for the new motif
    consensus.new <-as.vector(description.table[as.numeric(motifs.info[[new]][["number"]]),"consensus"])
    motifs.info[[new]][["strand"]] <- "D"
    motifs.info[[new]][["consensus"]] <- consensus.new
    
    ## Reset the offset
    if(case %in% c(1,3,5,8)){

      if(case == 3){
        ids <- get.id(n2)
        inverted.aligment.ids <- inverted.aligment(ids)
        for(id in names(inverted.aligment.ids)){
          motifs.info[[id]] <- inverted.aligment.ids[[id]]
        }
        new.spacer <- as.numeric(motifs.info[[new]][["spacer"]])
      }

      if(case %in% c(1,3)){
        spacer.diff <- (aligned.spacer - new.spacer)
      } else if (case %in% c(5,8)){
        spacer.diff <- (new.spacer - aligned.spacer)
      }
      
      offset <- offset + spacer.diff
      
    } else if(case %in% c(2,4,6,7)){

      ## Get the ids of the aligment that will be inverted
      if(case %in% c(2,4)){
        ids <- get.id(n2)
      } else if(case %in% c(6,7)){
        ids <- get.id(n1)
      }
      
      ## Invert the aligment and store the information in a list
      inverted.aligment.ids <- inverted.aligment(ids)
      
      ## Change the information in motifs.info list
      for(id in names(inverted.aligment.ids)){
        motifs.info[[id]] <- NULL
        motifs.info[[id]] <- inverted.aligment.ids[[id]]
      }
      
      aligned.spacer <- as.numeric(motifs.info[[aligned]][["spacer"]])
      new.spacer <- as.numeric(motifs.info[[new]][["spacer"]])

      if(case %in% c(6,7)){
        length.diff <- nchar(as.vector(description.table[as.numeric(motifs.info[[new]][["number"]]), "consensus"])) - nchar(as.vector(description.table[as.numeric(motifs.info[[aligned]][["number"]]), "consensus"]))
        spacer.diff <- (new.spacer - aligned.spacer)
      } else if(case %in% c(2,4)){
        length.diff <- nchar(as.vector(description.table[as.numeric(motifs.info[[aligned]][["number"]]), "consensus"])) - nchar(as.vector(description.table[as.numeric(motifs.info[[new]][["number"]]), "consensus"]))
        spacer.diff <- (aligned.spacer - new.spacer)
      }

      offset <- length.diff - offset + spacer.diff
      tree$labels[as.numeric(motifs.info[[new]][["number"]])] <- paste(motifs.info[[new]][["consensus"]], as.numeric(motifs.info[[new]][["number"]]))
    }

    ## Create the spacer
    spacer <- paste(collapse="",rep(x="-",times = abs(offset)))

    ##
    if(offset <= 0){
    
      for (id in get.id(n1)){
        motifs.info[[id]][["consensus"]] <- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
    else{
      
      for (id in get.id(n2)){
        motifs.info[[id]][["consensus"]] <- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
  }


  ##########################################
  ## Case 3: merging between two clusters ##
  ##########################################
  if ((child1 > 0) && (child2 > 0)) {

    N1 <- abs(min(child1, child2))
    N2 <- abs(max(child1, child2))
    
    n1 <- merge.levels.leaves[[N1]]
    n2 <- merge.levels.leaves[[N2]]

    ## Get ids
    ids1 <- get.id(n1)
    ids2 <- get.id(n2)
    
    ## Find the central motif of the cluster
    central.motifs <- central.motifs.ids(ids1, ids2)
    id1 <- central.motifs[1]
    id2 <- central.motifs[2]
    
    ## Get the comparison number in the compare-matrices table
    compa.nb <- get.compa.nb(id1,id2)
    
    ## Get the strand
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])

    ## Get the previous orientation of the aligned motifs
    prev.strand.1 <- motifs.info[[id1]][["strand"]]
    prev.strand.2 <- motifs.info[[id2]][["strand"]]
    change.offset <- 0

    ## Get the offset
    offset <- as.vector(compare.matrices.table[compa.nb, "offset"])

    ## Switch the ids
    if(id1 %in% ids2 == TRUE){
      temporal <- ids2
      ids1 <- ids2
      ids2 <- temporal
    }

    ## Get the current spacer of both motifs
    cluster.1.spacer <- as.numeric(motifs.info[[id1]][["spacer"]])
    cluster.2.spacer <- as.numeric(motifs.info[[id2]][["spacer"]])
    case <- 0

    ## Assign the value for the cases
    if(strand == "D"){
      if(prev.strand.1 == "R" && prev.strand.2 == "R"){
        case <- 1	
      } else if(prev.strand.1 == "R" && prev.strand.2 == "D"){
        case <- 2	
      } else if(prev.strand.1 == "D" && prev.strand.2 == "R"){
        case <- 3	
      } else if(prev.strand.1 == "D" && prev.strand.2 == "D"){
        case <- 4	
      }
    }else{
      if(prev.strand.1 == "R" && prev.strand.2 == "R"){
        case <- 5	
      } else if(prev.strand.1 == "R" && prev.strand.2 == "D"){
        case <- 6	
      } else if(prev.strand.1 == "D" && prev.strand.2 == "R"){
        case <- 7	
      } else if(prev.strand.1 == "D" && prev.strand.2 == "D"){
        case <- 8	
      }
    }
    
    ## Cases in which is required invert the aligment
    if(case %in% c(2,3,8)){
      
      ## Get the ids of the aligment that will be inverted
      ids <- get.id(n2)
      
      ## Invert the aligment and store the information in a list
      inverted.aligment.ids <- inverted.aligment(ids)
      
      ## Change the information in motifs.info list
      for(id in names(inverted.aligment.ids)){
        motifs.info[[id]] <- inverted.aligment.ids[[id]]
      }
      
      cluster.1.spacer <- as.numeric(motifs.info[[id1]][["spacer"]])
      cluster.2.spacer <- as.numeric(motifs.info[[id2]][["spacer"]])  
    }
    
    ## According to the cases, reset the offset
    if(case %in% c(1,2,5,6)){
      offset <- nchar(as.vector(description.table[as.numeric(motifs.info[[id1]][["number"]]), "consensus"])) - nchar(as.vector(description.table[as.numeric(motifs.info[[id2]][["number"]]), "consensus"]))  - offset + (cluster.1.spacer - cluster.2.spacer)
    } else if(case %in% c(3,4,7,8)){
      offset <- offset + (cluster.1.spacer - cluster.2.spacer)
    }
    
    ## Add the spacer to the motifs
    if(offset <= 0){
      spacer <- paste(collapse="",rep(x="-",times = abs(offset)))
      
      for (id in ids1){
        motifs.info[[id]][["consensus"]] <- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
    else{
      spacer <- paste(collapse="",rep(x="-",times = abs(offset)))
      
      for (id in ids2){
        motifs.info[[id]][["consensus"]] <- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
  }
}

#jpeg(filename = "/home/jaimecm/Documents/TAGC/Clustering_test/Test_different_hclust_methods/results/testing.jpeg")
dev.new(width=10, height=7)
#par(mar=c(3,2,1,30),family="mono")
par(mar=c(3,2,1,20),family="mono")
plot(as.dendrogram(tree), horiz=TRUE)
#dev.off()




## Produce the aligment table
alignment.table <- lapply(motifs.info, function(X){
  return(c(X[["number"]], X[["consensus"]], X[["strand"]], X[["spacer"]]))
})
alignment.table <- as.data.frame(t(data.frame(alignment.table)))
alignment.table$id <- as.vector(rownames(alignment.table))
alignment.table <- alignment.table[,c(5,1:4)]
colnames(alignment.table) <- c("#id", "number", "consensus", "strand", "spacer")
write.table(alignment.table, "/home/jaimecm/Documents/TAGC/Clustering_test/Test_different_hclust_methods/results/testing.txt", sep = "\t", quote = FALSE, row.names = FALSE)
