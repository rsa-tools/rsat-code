###################################################################
## Library with R functions required in R-scripts/cluster_motifs.R

################################################################
## Check parameter values
check.param <- function() {
  ## Check that input file has been specified
  if (!exists("infile")) {
    stop("Missing mandatory argument: infile=[matrix_comparison_table] ")
  }
  verbose(paste("Input file", infile), 1)
  
  ## Check that description file
  if (!exists("description.file")) {
    stop("Missing mandatory argument: description.file=[matrix_description_table] ")
  }
  verbose(paste("Description file", description.file), 1)
  
  ## Check that output file has been specified
  if (!exists("out.prefix")) {
    stop("Missing mandatory argument: out.prefix=[output_prefix] ")
  }
  verbose(paste("Output prefix", out.prefix), 1)
  
  ## Check that distance table file has been specified
  if (!exists("distance.table")) {
    distance.table <- paste(sep="", out.prefix, "_dist_table.tab")
  }
  verbose(paste("Distance table", distance.table), 1)
  
  ## Default score is the normalized correlation
  if (!exists("score")) {
    score <- "Ncor";
  }
  
  ## Default hclust method is the complete method
  if (!exists("hclust.method")) {
    hclust.method <- "average";
  }
}

################################################################
## Convert user-selected score to distance metrics
##
## Each score requires to be treated according to its nature
## (similarity or distance) plus some specificities (correlation goes
## from -1 to 1, ...).
convert.scores <- function (score.values,
                            score) {
  ## Similarity sores bounded to 1
  if ((score == "Ncor")
      || (score=="cor")
      || (score=="logocor")
      || (score=="Nlocogor")
      || (score=="Icor")
      || (score=="NIcor")
  ) {
    ## cor   		Pearson correlation (computed on residue occurrences in aligned columns)
    ## Ncor 			Relative width-normalized Pearson correlation
    ## logocor 			correlation computed on sequence logos
    ## Nlogocor 			Relative width-normalized logocor
    ## Icor 			Pearson correlation computed on Information content
    ## NIcor 			Relative width-normalized Icor
    score.dist <- 1 - score.values
    
  } else if ((score == "logoDP")
             || (score == "cov")) {
    ## logoDP 			dot product of sequence logos
    ## cov 			covariance between residues in aligned columns
    
    stop("logoDP and cov scores are not supported yet")
    
  } else if ((score == "dEucl")
             || (score == "NdEucl")
             || (score == "NdEucl")
             || (score == "NsEucl")
             || (score == "SSD")
             || (score == "SW")
             || (score == "NSW")
  ) {
    ## dEucl 			Euclidian distance between residue occurrences in aligned columns
    ## NdEucl 			Relative width-normalized dEucl
    ## NsEucl 			similarity derived from Relative width-normalized Euclidian distance
    ## SSD 			Sum of square deviations
    ## SW 			Sandelin-Wasserman
    ## NSW 			Relative width-normalized Sandelin-Wasserman
    
    score.dist <- score.values
    
  } else if (score == "match_rank") {
    ## match_rank rank of current match among all sorted matches
    stop("match_rank score is not supported yet")
    
  } else {
    stop(paste(score, "is not a valid score", sep="\t"))
  }
}


#######################################################
## Extract the tree from an hclust object 
createLeafNode <- function(hclust, i) {
  list(label = hclust$labels[[i]],
       order = hclust$order[[i]])
}

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
## Get the number of comparison between the input IDs
## in the compare-matrices results table
get.compa.nb <- function(ID1,ID2){
  
  return(which( (compare.matrices.table[,"id2"] == ID2 & compare.matrices.table[,"id1"] == ID1) | (compare.matrices.table[,"id1"] == ID2 & compare.matrices.table[,"id2"] == ID1) ))
}


########################################################
## Given the number of leaf, get the id of the motif
get.id <- function(num){
  
  return(as.vector(sapply(num, function(X){
    description.table[X,"id"]
  })))
}


########################################################
## Given a vector with IDs, return a list with the 
## information (consensus, number, id, strand, spacer)
## of the inverted alignment
inverted.aligment <- function(ids, motifs.info){

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

  #names(temporal.list) <- sapply(strsplit(names(temporal.list), "\\."), function(X){return(X[1])})
  names(temporal.list) <- ids
  return(temporal.list)
}


##################################################
## Align two leaves: creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.two.leaves <- function(child1, child2, motifs.info, tree){

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

    export.list <- list()
    export.list[["info"]] <- motifs.info
    export.list[["tree"]] <- tree
    return(export.list)
 }


##############################################################
## Align a leaf and one cluster: align the single motif relative
## to the already aligned cluster; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.leave.and.cluster <- function(child1, child2, motifs.info, tree){

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
      print(paste("##### ", ids, " #####"))
      inverted.aligment.ids <- inverted.aligment(ids, motifs.info)
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
    inverted.aligment.ids <- inverted.aligment(ids, motifs.info)
    
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

  export.list <- list()
  export.list[["info"]] <- motifs.info
  export.list[["tree"]] <- tree
  return(export.list)
}


##############################################################
## Align two clusters: align the two cluster, in some cases is
## necessary invert the alignment; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.clusters <- function(child1, child2, motifs.info, tree){

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

  export.list <- list()
  export.list[["info"]] <- motifs.info
  export.list[["tree"]] <- tree
  return(export.list)
}


###################################################
## Set the sizes of all consensuses to the largest
## one, adding "-" at the end of the consensuses
fill.downstream <- function(motifs.list){

  ## Get the highest length among the consensuses
  consensuses <- sapply(motifs.list, function(X){
    nchar(X[["consensus"]])
  })
  max.cons.length <- max(consensuses)

  for(id in 1:length(motifs.list)){

    spacer.length <-  max.cons.length - nchar(motifs.list[[id]][["consensus"]])
    spacer <- paste(rep("-", times = spacer.length), collapse = "")
    motifs.list[[id]][["consensus"]] <- paste( motifs.list[[id]][["consensus"]], spacer, sep = "")
    motifs.list[[id]][["offset_down"]] <- spacer.length
    
  }
  return(motifs.list)
}
