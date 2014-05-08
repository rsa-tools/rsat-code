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
    distance.table <<- paste(sep="", out.prefix, "_dist_table.tab")
  }
  verbose(paste("Distance table", distance.table), 1)
  
  ## Default score is the normalized correlation
  if (!exists("score")) {
    score <<- "Ncor";
  }
  
  ## Default hclust method is the complete method
  if (!exists("hclust.method")) {
    hclust.method <<- "average";
  }

  ## Define the kind of metric used: scores or distances
  supported.scores <- c("cor", "Ncor")
  supported.distances <- c(NULL)
  
  if(score %in% supported.scores){
    metric <<- "similarity"
    
    ## Default lower and upper thresholds equals to zero
    if (!exists("lth")) {
      lth <<- list()
      lth[["Ncor"]] <<- 0;
      lth[["cor"]] <<- 0;

      lth.values <<- unlist(lth)
      lth.scores <<- names(lth.values)
    }
    if (!exists("uth")) {
      uth <<- list()
      uth[["Ncor"]] <<- 1;
      uth[["cor"]] <<- 1;

      uth.values <<- unlist(uth)
      uth.scores <<- names(uth.values)
    }
    
  } else if(score %in% supported.distances){
    metric <<- "distances"

    ## ## Default lower and upper thresholds equals to zero
    ## if (!exists("lth")) {
    ##   lth <<- 1;
    ## }
    ## if (!exists("uth")) {
    ##   uth <<- 0;
    ## }
  }

  
  #####################
  ##
  if(exists("lthsp")){
    lthsp <- unlist(strsplit(lthsp, "-"))
    lth.scores <<- lthsp[seq(1,length(lthsp), by = 2)]
    lth.values <<- as.numeric(lthsp[seq(2,length(lthsp), by = 2)])
    
    supported <- c("Ncor", "cor")
    if(length(setdiff(supported, lth.scores)) > 0){
      for(add in setdiff(supported, lth.scores)){
        lth.scores <<- append(lth.scores, add)
        lth.values <<- append(lth.values, 0) 
      }
    }
  }

  ## Set the labels
  labels <<- unique(labels)
  print(paste("### ", labels, " ###"))
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
leaves.per.node <- function (tree) {
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
inverted.alignment <- function(ids){

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
    inverted.aligment.list[[X]][["consensus"]] <- new.consensus
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
align.two.leaves <- function(child1, child2){
  
  ## Identify the two motifs
  n1 <- min(-child1,-child2) ## row number of the first motif in the description table
  n2 <- max(-child1,-child2) ## row number of the second motif in the description table
  
  id1 <- get.id(n1) ## Id of the first motif
  id2 <- get.id(n2) ## Id of the second motif

  ## Nodes attributes
  attributes <- attributes.among.clusters(id1, id2)
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["cluster_1"]] <<- paste(n1, collapse = " ")
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["cluster_2"]] <<- paste(n2, collapse = " ")
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["min_score"]] <<- attributes[1]
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["max_score"]] <<- attributes[2]
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["median_score"]] <<- attributes[3]
  
  ## Comparison number in the compare-matrices table
  compa.nb <- get.comparison.number(id1,id2)
  
  ## Check the threshold values for the corresponding
  ## metric used
  aligned.motif.flag <- 0
  aligned.motif.flag <- alignment.status(id1, id2, hclust.method)

  ## In case the motifs should not be aligned
  ## fill the motifs.info list with the default parameters
  if(aligned.motif.flag == 0){
    
    internal.nodes.attributes[[paste("merge_level_", merge.level)]][["alignment_status"]] <<- "Non-aligned"
    
    for(n in c(n1, n2)){
      motifs.info[[get.id(n)]][["strand"]] <<- "D"
      motifs.info[[get.id(n)]][["number"]] <<- n
      motifs.info[[get.id(n)]][["spacer"]] <<- 0
      motifs.info[[get.id(n)]][["consensus"]] <<- as.vector(description.table[as.numeric(motifs.info[[get.id(n)]][["number"]]),"consensus"])
    }

    forest.nb <<- forest.nb + 1

  ## Conversely align the motifs
  }else{
    
    internal.nodes.attributes[[paste("merge_level_", merge.level)]][["alignment_status"]] <<- "Aligned"
    
    ## Choose the relative orientation of the two motifs
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])
    consensus1 <- as.vector(description.table[n1,"consensus"]) ## Consensus of the first motif
    motifs.info[[id1]][["strand"]] <<- "D" 
    if (strand == "R") {
      consensus2 <- as.vector(description.table[n2,"rc_consensus"]) ## Consensus of the second motif
      motifs.info[[id2]][["strand"]] <<- "R"
    } else {
      consensus2 <- as.vector(description.table[n2,"consensus"]) ## Consensus of the second motif
      motifs.info[[id2]][["strand"]] <<- "D"
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
    tree$labels[n1] <<- paste(consensus1, n1)      
    tree$labels[n2] <<- paste(consensus2, n2)
    
    ## Store the strand of the cluster, the consensuses and numbers
    motifs.info[[id1]][["consensus"]] <<- consensus1
    motifs.info[[id1]][["number"]] <<- n1
    
    motifs.info[[id2]][["consensus"]] <<- consensus2
    motifs.info[[id2]][["number"]] <<- n2
    
    motifs.info[[id1]][["spacer"]] <<- length(unlist(strsplit(motifs.info[[id1]][["consensus"]], "-")))-1
    motifs.info[[id2]][["spacer"]] <<- length(unlist(strsplit(motifs.info[[id2]][["consensus"]], "-")))-1
  }
}


##############################################################
## Align a leaf and one cluster: align the single motif relative
## to the already aligned cluster; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.leave.and.cluster <- function(child1, child2){
  
  n1 <- abs(min(child1, child2))
  n.aligned <- merge.levels.leaves[[merge.level]][which(merge.levels.leaves[[merge.level]] != n1)]
  n2 <- n.aligned

  ## Nodes attributes
  attributes <- attributes.among.clusters(get.id(n1), get.id(n2))
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["cluster_1"]] <<- paste(n1, collapse = " ")
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["cluster_2"]] <<- paste(n2, collapse = " ")
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["min_score"]] <<- attributes[1]
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["max_score"]] <<- attributes[2]
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["median_score"]] <<- attributes[3]
  
  ## Get ids
  motifs.info[[get.id(n1)]][["number"]] <<- as.numeric(n1)
  motifs.info[[get.id(n1)]][["spacer"]] <<- 0
  ids.aligned <- get.id(n.aligned)
  
  ## Find the central motif of the cluster
  central.motifs <- central.motifs.ids(get.id(n1), ids.aligned)
  id1 <- central.motifs[1]
  id2 <- central.motifs[2]
  
  ## Get the comparison number in the compare-matrices table
  compa.nb <- get.comparison.number(id1,id2)
  
  ## Check the threshold values for the corresponding
  ## metric used
  aligned.motif.flag <- 0  
  aligned.motif.flag <- alignment.status(id1, id2, hclust.method)

  ## In case the motifs should not be aligned
  ## fill the motifs.info list with the default parameters
  if(aligned.motif.flag == 0){
    
    motifs.info[[get.id(n1)]][["strand"]] <<- "D"
    motifs.info[[get.id(n1)]][["consensus"]] <<- as.vector(description.table[as.numeric(motifs.info[[get.id(n1)]][["number"]]),"consensus"])
    
    internal.nodes.attributes[[paste("merge_level_", merge.level)]][["alignment_status"]] <<- "Non-aligned"

    forest.nb <<- forest.nb + 1
    
  ## Conversely align the motifs
  } else{
    
    internal.nodes.attributes[[paste("merge_level_", merge.level)]][["alignment_status"]] <<- "Aligned"
    
    ## Get the strand
    strand <- as.vector(compare.matrices.table[compa.nb, "strand"])
    
    ## Identified the new and the aligned motif
    switch.ids <- 0
    new <- NULL
    aligned <- NULL 
    if(id1 != get.id(n1)){
      aligned <- id1
      new <- id2
      temporal <- NULL
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
    
    ## Invert the aligned motifs
    ## This is just required in some cases
    ## See table at the final comments
    if(case %in% c(2,4,6,7)){
      ids <- get.id(n.aligned)
      inverted.alignment.ids <- inverted.alignment(ids)
      for(id in names(inverted.alignment.ids)){
        motifs.info[[id]] <<- NULL
        motifs.info[[id]] <<- inverted.alignment.ids[[id]]
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <<- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
    
    ## Get the spacers
    aligned.spacer <- as.numeric(motifs.info[[aligned]][["spacer"]])
    new.spacer <- as.numeric(motifs.info[[new]][["spacer"]])
    
    ## Choose the consensus for the new motif
    if(case %in% c(1,2,5,6,7,8)){
      consensus.new <- as.vector(description.table[as.numeric(motifs.info[[new]][["number"]]),"consensus"])
      motifs.info[[new]][["strand"]] <<- "D"
      motifs.info[[new]][["consensus"]] <<- consensus.new
      
    } else if(case %in% c(3,4)){
      consensus.new <- as.vector(description.table[as.numeric(motifs.info[[new]][["number"]]),"rc_consensus"])
      motifs.info[[new]][["strand"]] <<- "R"
      motifs.info[[new]][["consensus"]] <<- consensus.new
    }
    tree$labels[as.numeric(motifs.info[[new]][["number"]])] <<- paste(motifs.info[[new]][["consensus"]], as.numeric(motifs.info[[new]][["number"]]))
    
    ## Reset the offset
    if(case %in% c(1:4)){
      spacer.diff <- (aligned.spacer - new.spacer)
    } else if(case %in% c(5:8)){
      spacer.diff <- (new.spacer - aligned.spacer)
    }
    offset <- offset + spacer.diff
    
    ## Create the spacer
    spacer <- paste(collapse="",rep(x="-",times = abs(offset)))
    
    ## Add the gaps
    if(offset < 0){
      for (id in get.id(n1)){
        motifs.info[[id]][["consensus"]] <<- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <<- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <<- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    } else{
      for (id in get.id(n2)){
        
        motifs.info[[id]][["consensus"]] <<- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <<- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <<- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
  }
}


##############################################################
## Align two clusters: align the two cluster, in some cases is
## necessary invert the alignment; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.clusters <- function(child1, child2){
  
  N1 <- abs(min(child1, child2))
  N2 <- abs(max(child1, child2))
    
  n1 <- merge.levels.leaves[[N1]]
  n2 <- merge.levels.leaves[[N2]]
  
  ## Get ids
  ids1 <- get.id(n1)
  ids2 <- get.id(n2)

  ## Nodes attributes
  attributes <- attributes.among.clusters(ids1, ids2)
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["cluster_1"]] <<- paste(n1, collapse = " ")
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["cluster_2"]] <<- paste(n2, collapse = " ")
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["min_score"]] <<- attributes[1]
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["max_score"]] <<- attributes[2]
  internal.nodes.attributes[[paste("merge_level_", merge.level)]][["median_score"]] <<- attributes[3]
  
  ## Find the central motif of the cluster
  central.motifs <- central.motifs.ids(ids1, ids2)
  id1 <- central.motifs[1]
  id2 <- central.motifs[2]
  
  ## Get the comparison number in the compare-matrices table
  compa.nb <- get.comparison.number(id1,id2)
  
  ## Get the strand
  strand <- as.vector(compare.matrices.table[compa.nb, "strand"])

  ## Check the threshold values for the corresponding
  ## metric used
  aligned.motif.flag <- 0
  score.val <- compare.matrices.table[compa.nb, score]
  aligned.motif.flag <- alignment.status(id1, id2, hclust.method)
  if(aligned.motif.flag == 1){
    
    internal.nodes.attributes[[paste("merge_level_", merge.level)]][["alignment_status"]] <<- "Aligned"
    
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
      inverted.alignment.ids <- inverted.alignment(ids)
      
      ## Change the information in motifs.info list
      for(id in names(inverted.alignment.ids)){
        motifs.info[[id]] <<- NULL 
        motifs.info[[id]] <<- inverted.alignment.ids[[id]]
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
        motifs.info[[id]][["consensus"]] <<- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <<- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
    else{
      spacer <- paste(collapse="",rep(x="-",times = abs(offset)))
      
      for (id in ids2){
        motifs.info[[id]][["consensus"]] <<- paste(spacer, motifs.info[[id]][["consensus"]], sep="")
        motifs.info[[id]][["spacer"]] <<- length(unlist(strsplit(motifs.info[[id]][["consensus"]], "-")))-1
        tree$labels[as.numeric(motifs.info[[id]][["number"]])] <<- paste(motifs.info[[id]][["consensus"]], as.numeric(motifs.info[[id]][["number"]]))
      }
    }
  }else {
    internal.nodes.attributes[[paste("merge_level_", merge.level)]][["alignment_status"]] <<- "Non-aligned"
    forest.nb <<- forest.nb + 1
  }
}
  

###################################################
## Set the sizes of all consensuses to the largest
## one, adding "-" at the end of the consensuses
fill.downstream <- function(motifs.list, forest.ids){

  ## Saves the max width of the motifs belongig to each forest
  max.width.forest <- sapply(forest.ids, function(X){

    ## Get the highest length among the consensuses
    consensuses <- sapply(X, function(Y){
      nchar(motifs.list[[Y]][["consensus"]])
    })
    max.cons.length <- max(consensuses)
  })

  ## Fill the spaces
  for(i in 1:length(max.width.forest)){
    for(id in forest.ids[[paste("forest_", i, sep = "")]]){
      spacer.length <-  max.width.forest[i] - nchar(motifs.list[[id]][["consensus"]])
      spacer <- paste(rep("-", times = spacer.length), collapse = "")
      motifs.list[[id]][["consensus"]] <- paste( motifs.list[[id]][["consensus"]], spacer, sep = "")
      motifs.list[[id]][["offset_down"]] <- spacer.length
    }
  }
  return(motifs.list)
}


################################################################
## Depending on the agglomeratiom method, single, complete or
## average, indicates if the motif or cluster will be aligned
alignment.status <- function(ids1, ids2, aglomeration.method){

  alignment.status <- 0
  if(aglomeration.method == "single"){
    alignment.status <- alignment.test.method.single(ids1, ids2)
  } else if(aglomeration.method == "complete"){
    alignment.status <- alignment.test.method.complete(ids1, ids2)
  } else if(aglomeration.method == "average"){
    alignment.status <- alignment.test.method.average(ids1, ids2)
  }
  return(alignment.status)
}


###################################################
## Get the comparison numbers for all the pairs
## between the leaves of the two clusters
get.comparison.numbers.between.ids <- function(ids1, ids2){

  compa.numbers <- sapply(ids2, function(X){
    sapply(ids1, function(Y){
      get.comparison.number(X, Y)
    })
  })
  compa.numbers <- as.vector(compa.numbers)
  return(compa.numbers)
}


#####################################################################
## Evaluate if at least on pair of leaves between the two clusters
## is upper/lower than the specified thresholds,
## if so, align the two clusters
alignment.test.method.single <- function(ids1, ids2){

  compa.numbers <- get.comparison.number (ids1, ids2)

  ## Get the scores of the comparisons
  scores <- compare.matrices.table[compa.numbers,score]

  if(metric == "similarity"){
    if(pairs.satisfy.lth() > 0){
      return(1)
    }else{
      return(0)
    }
  } else if (metric == "distance"){
    if(pairs.satisfy.uth() > 0){
      return(1)
    }else{
      return(0)
    }
  }
}


#####################################################################
## Evaluate if all pairs of leaves between the two clusters have an
## score upper/lower than the threshold,
##if so, align the two clusters
alignment.test.method.complete <- function(ids1, ids2){

  compa.numbers <- get.comparison.number (ids1, ids2)

  ## Get the scores of the comparisons
  scores <- compare.matrices.table[compa.numbers,score]
  #scores <- sapply(scores, function(X){ if(X == 0){X <- NA}else {X <- X}})  

  ## Count the number of pairs which satisfied the condition
  if(metric == "similarity"){
    sat.cond <- pairs.satisfy.lth()
  } else if(metric == "distances"){
    sat.cond <- pairs.satisfy.uth
  }

  ## If all the pairs satisfied the condition,
  ## then align the clusters
  if(sat.cond == length(scores)){
    return(1)
  } else{
    return(0)
  }
}


#####################################################################
## Evaluate if the mean of scores for all pairs of leaves between
## the two clusters is upper/lower than the threshold,
##if so, align the two clusters
alignment.test.method.average <- function(ids1, ids2){

  compa.numbers <- get.comparison.number(ids1, ids2)
  
  ## Get the scores of the comparisons
  if (metric == "similarity"){
    scores <- compare.matrices.table[compa.numbers, lth.scores]
  } else if(metric == "distance"){
    scores <- compare.matrices.table[compa.numbers, uth.scores]
  }
  #scores <- sapply(scores, function(X){ if(X == 0){X <- NA}else {X <- X}})  
  
  ## Calculate the median of the data
  median.scores<- apply(scores, 2, median)

  ## According to the kind of metric selected
  ## evaluates if the clusters will be aligned
  if(metric == "distance"){
    if(sum(median.scores <= uth.values) == length(uth.values)){
      return(1)
    }else{
      return(0)
    }
  } else if (metric == "similarity"){
    if(sum(median.scores >= lth.values) == length(lth.values)){
      return(1)
    }else{
      return(0)
    }
  }
}


########################################################
## Get the number of comparison between the input IDs
## in the compare-matrices results table
get.comparison.number <- function(ids1, ids2){

  compa.numbers <- sapply(ids2, function(X){
    sapply(ids1, function(Y){
      which( (compare.matrices.table[,"id2"] == X & compare.matrices.table[,"id1"] == Y) | (compare.matrices.table[,"id1"] == X & compare.matrices.table[,"id2"] == Y) )
    })
  })
  compa.numbers <- as.vector(compa.numbers)
  return(compa.numbers)
}


################################################################
## Identify the "central" leaf (motif) of a subtre (cluster),
## i.e. the motif with the smallest average distance to all
## other motifs of this cluster.
central.motifs.ids <- function(ids1, ids2){
  
  ## Get the comparison numbers
  compa.numbers <- sapply(ids2, function(X){
    sapply(ids1, function(Y){
      get.comparison.number(X, Y)
    })
  })
  compa.numbers <- as.vector(compa.numbers)

  ## Get the ids of the less distant nodes
  compa.info <- compare.matrices.table[compa.numbers,][which(compare.matrices.table[compa.numbers,score] == max(compare.matrices.table[compa.numbers,score])),c("id1","id2")][1,]
  return(as.vector(unlist(compa.info)))
}


##############################################
## Given two set of IDs belonging to two
## different clusters: return the min score
attributes.among.clusters <- function(ids1, ids2){

  compa.numbers <- get.comparison.number(ids1, ids2)
  scores <- compare.matrices.table[compa.numbers,score]
  min.score <- min(scores)
  max.score <- max(scores)
  median.score <- median(scores)
  
  return(c(min.score, max.score, median.score))
}


#################################
## Build a distance matrix from
## the distance score list
build.distance.matrix <- function(comparison.table){

  dist.table <- NULL
  distances.objects <- list()

  ## Extract score values
  score.values <- comparison.table[,score] 
  score.dist <- 1 - score.values
  
  ## Add a column with score column to the compare matrices table, will
  ## be used to generate a cross-table
  comparison.table$score.dist <- score.dist
  
  dist.table <- t(xtabs(score.dist ~ name1+name2, comparison.table) )
  ## Ensure that symmetrical distances are defined
  for (i in 1:nrow(dist.table)) {
    for (j in i:ncol(dist.table)) {
      if (i==j) {next}
      
      dist.max <- max(dist.table[i,j], dist.table[j,i])
      dist.table[i,j] <- dist.max
      dist.table[j,i] <- dist.max
    }
  }
  
  dist.table <- dist.table
  dist.matrix <- as.dist(dist.table)

  return(list(dist.table, dist.matrix))
}


######################################################
## Fill with gaps those alignments within  a forest
fill.downstream.forest <- function(motifs.list){

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


###############################################################
## Count the numbers of pairs which satisfied the thresholds
pairs.satisfy.lth <- function(){
  
  values.list <- list()
  for(i in 1:length(lth.scores)){
    values.list[[lth.scores[i]]] <- which(scores >= lth.values[i])
  }
  return(length(Reduce(intersect, values.list)))
}

pairs.satisfy.uth <- function(){
  
  values.list <- list()
  for(i in 1:length(lth.scores)){
    values.list[[lth.scores[i]]] <- which(scores <= lth.values[i])
  }
  return(length(Reduce(intersect, values.list)))
}


################
### Inefficient (but probably more robust) way to build a distance matrix, involving two loops
## matrix.names <- sort(unique(as.vector(as.matrix((compare.matrices.table[,c(1,2)])))))
## matrix.nb <- length(matrix.names)
## dist.matrix <- data.frame(matrix(NA, nrow=matrix.nb, ncol=matrix.nb))
## names(dist.matrix) <- matrix.names
## rownames(dist.matrix) <- matrix.names
## for (i in 1:matrix.nb) {
##   for (j in i:matrix.nb) {
##     ## Identify the relevant rows in matrix comparison table
##     matrix.id1 <- matrix.names[i]
##     matrix.id2 <- matrix.names[j]
##     matches <- ((compare.matrices.table[,1] == matrix.id1 & compare.matrices.table[,2] == matrix.id2) |
##                 (compare.matrices.table[,1] == matrix.id2 & compare.matrices.table[,2] == matrix.id1))
##     if (sum(matches) >= 1) {
##       current.score <- compare.matrices.table[which(matches),score]
##     } else {
##       current.score <- "NA"
##     }    
##     dist.matrix[i,j] <- current.score
##     dist.matrix[j,i] <- current.score
##   }
## }


