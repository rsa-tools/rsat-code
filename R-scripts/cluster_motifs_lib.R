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

    for(i in 1:length(lth.scores)){
      thresholds[[lth.scores[i]]] <<- lth.values[i]
    }
    
    supported <- c("cor", "Ncor")
    if(length(setdiff(supported, lth.scores)) > 0){
      for(add in setdiff(supported, lth.scores)){
        lth.scores <<- append(lth.scores, add)
        lth.values <<- append(lth.values, 0) 
      }
    }
  }

  ## Set the labels
  labels <<- unique(labels)
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

    ## Depending on whether the right branch points to a leave or an
    ## iternal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch2 < 0) {
      nodes2 <- -branch2 ## branch one only contains one leave
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

  id1.test <- get.id(n1)
  id2.test <- get.id(n2)

  ## Get the id of each node
  central.motifs <- central.motifs.ids(get.id(n1), get.id(n2))
  id1 <- central.motifs[1]
  id2 <- central.motifs[2]

  if(id1.test == id2){
    temporal <- id1
    id1 <- id2
    id2 <- temporal
    temporal <- n1
    n1 <- n2
    n2 <- temporal
  }
  
  
  ## Comparison number in the compare-matrices table
  compa.nb <- get.comparison.number(id1, id2)[1]
  
  ## Check the threshold values for the corresponding
  ## metric used
  aligned.motif.flag <- 0
  aligned.motif.flag <- alignment.status(id1.test, id2.test, hclust.method)
  print(paste("### FLAG: ", aligned.motif.flag, " ### METHOD: ", hclust.method," ###", sep = ""))
  
  ## In case the motifs should not be aligned
  ## fill the motifs.info list with the default parameters
  if(aligned.motif.flag == 0){
    
    for(n in c(n1, n2)){
      motifs.info[[get.id(n)]][["strand"]] <<- "D"
      motifs.info[[get.id(n)]][["number"]] <<- n
      motifs.info[[get.id(n)]][["spacer"]] <<- 0
      motifs.info[[get.id(n)]][["consensus"]] <<- as.vector(description.table[as.numeric(motifs.info[[get.id(n)]][["number"]]),"consensus"])
    }

  ## Conversely align the motifs
  }else{
    
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
  rm(n1, n2)
}


##############################################################
## Align a leaf and one cluster: align the single motif relative
## to the already aligned cluster; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.leave.and.cluster <- function(child1, child2){
  
  n1 <- abs(min(child1, child2))
  n.aligned <- merge.levels.leaves[[merge.level]][which(merge.levels.leaves[[merge.level]] != n1)]
  n2 <- n.aligned
  
  ## Get ids
  motifs.info[[get.id(n1)]][["number"]] <<- as.numeric(n1)
  motifs.info[[get.id(n1)]][["spacer"]] <<- 0
  ids.aligned <- get.id(n.aligned)

  
  print(paste("##################### IDS:", paste(get.id(n2), collapse = ", "), "###", sep = ""))
  print(paste("##################### IDS:", paste(ids.aligned, collapse = ", "), "###", sep = ""))
  
  ## Check the threshold values for the corresponding
  ## metric used
  aligned.motif.flag <- 0  
  aligned.motif.flag <- alignment.status(get.id(n1), get.id(n2), hclust.method)

  print(paste("### FLAG: ", aligned.motif.flag, " ### Aligned: ", ids.aligned," ###", sep = ""))
  
  ## In case the motifs should not be aligned
  ## fill the motifs.info list with the default parameters
  if(aligned.motif.flag == 0){
    
    motifs.info[[get.id(n1)]][["strand"]] <<- "D"
    motifs.info[[get.id(n1)]][["consensus"]] <<- as.vector(description.table[as.numeric(motifs.info[[get.id(n1)]][["number"]]),"consensus"])
    
  ## Conversely align the motifs
  } else{

    ## Find the central motif of the cluster
     central.motifs <- central.motifs.ids(get.id(n1), get.id(n2))
    ##central.motifs <- central.motifs.ids(get.id(n1), ids.aligned)
    id1 <- central.motifs[1]
    id2 <- central.motifs[2]
    
    ## Get the comparison number in the compare-matrices table
    compa.nb <- get.comparison.number(id1,id2)
    
    print(paste("##################### COMPA:", compa.nb, " ###", sep = ""))
    
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

    ## print(paste("##################### OFFSET:", offset, "###", sep = ""))
    
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
    rm(n.aligned, n1, n2)
}


##############################################################
## Align two clusters: align the two cluster, in some cases is
## necessary invert the alignment; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.clusters <- function(child1, child2){

  flag1 <- 0
  
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
  compa.nb <- get.comparison.number(id1,id2)
  
  ## Get the strand
  strand <- as.vector(compare.matrices.table[compa.nb, "strand"])

  ## Check the threshold values for the corresponding
  ## metric used
  aligned.motif.flag <- 0
  aligned.motif.flag <- alignment.status(ids1, ids2, hclust.method)

  print(paste("### FLAG: ", aligned.motif.flag, " ### METHOD: ", hclust.method," ###", sep = ""))
  
  if(aligned.motif.flag == 1){
    
    ## Get the offset
    offset <- as.vector(compare.matrices.table[compa.nb, "offset"])
    
    ## Switch the ids
    if(id1 %in% ids2){
      temporal <- NULL
      temporal <- ids1
      ids1 <- ids2
      ids2 <- temporal
    }

    ## Get the previous orientation of the aligned motifs
    prev.strand.1 <- motifs.info[[id1]][["strand"]]
    prev.strand.2 <- motifs.info[[id2]][["strand"]]
    change.offset <- 0

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
      ids <- ids2
      
      ## Invert the aligment and store the information in a list
      inverted.alignment.ids <- inverted.alignment(ids)
      
      ## Change the information in motifs.info list
      for(id in names(inverted.alignment.ids)){
        motifs.info[[id]] <<- NULL 
        motifs.info[[id]] <<- inverted.alignment.ids[[id]]
      }

      if(id1 %in% ids2){
        cluster.1.spacer <- as.numeric(motifs.info[[id1]][["spacer"]])
        cluster.2.spacer <- as.numeric(motifs.info[[id2]][["spacer"]])
      } else {
        cluster.1.spacer <- as.numeric(motifs.info[[id2]][["spacer"]])
        cluster.2.spacer <- as.numeric(motifs.info[[id1]][["spacer"]])
      }
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
  }
  rm(N1, N2, n1, n2)
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
    for(id in forest.ids[[paste("cluster_", i, sep = "")]]){
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

  closest.motifs <- central.motifs.ids(ids1, ids2)
  id1.close <- closest.motifs[1]
  id2.close <- closest.motifs[2]
  
  compa.numbers <- get.comparison.number(id1.close, id2.close)

  ## Get the scores of the comparisons
  closest.scores <- compare.matrices.table[compa.numbers, names(thresholds)]

  metrics.up.threshold <- NULL
  for(i in 1:length(thresholds)){ 
    metrics.up.threshold <- append(metrics.up.threshold,sum(closest.scores[,names(thresholds[i])] >= thresholds[[i]]))
  }

  if(sum(metrics.up.threshold) == length(thresholds)){
    return(1)
  }else{
   return(0)
  }

  ## TO DO: implement thresholds for distance metrics
  ## else if (metric == "distance"){
  ##   sat.cond <- pairs.satisfy.uth()
  ## }

}


#####################################################################
## Evaluate if all pairs of leaves between the two clusters have an
## score upper/lower than the threshold,
## if so, align the two clusters
alignment.test.method.complete <- function(ids1, ids2){

  farthest.motifs <- farthest.motifs.ids(ids1, ids2)
  id1.far <- farthest.motifs[1]
  id2.far <- farthest.motifs[2]
  
  compa.numbers <- get.comparison.number(id1.far, id2.far)

  ## Get the scores of the comparisons
  farthest.scores <- compare.matrices.table[compa.numbers, names(thresholds)]

  metrics.up.threshold <- NULL
  for(i in 1:length(thresholds)){ 
    metrics.up.threshold <- append(metrics.up.threshold,sum(farthest.scores[,names(thresholds[i])] >= thresholds[[i]]))
  }

  if(sum(metrics.up.threshold) == length(thresholds)){
    return(1)
  }else{
   return(0)
  }

  ## TO DO: implement thresholds for distance metrics
  ## else if (metric == "distance"){
  ##   sat.cond <- pairs.satisfy.uth()
  ## }

}


################################################################
## Identify the "farthest" leaf (motif) of a subtre (cluster),
## i.e. the motif with the largest average distance to all
## other motifs of this cluster.
farthest.motifs.ids <- function(ids1, ids2){
  
  ## Get the comparison numbers
  compa.numbers <- sapply(ids2, function(X){
    sapply(ids1, function(Y){
      get.comparison.number(X, Y)
    })
  })
  compa.numbers <- as.vector(compa.numbers)

  ## Get the ids of the less distant nodes
  compa.info <- compare.matrices.table[compa.numbers,][which(compare.matrices.table[compa.numbers,score] == min(compare.matrices.table[compa.numbers,score])),c("id1","id2")][1,]
  return(as.vector(unlist(compa.info)))
}

#####################################################################
## Evaluate if the mean of scores for all pairs of leaves between
## the two clusters is upper/lower than the threshold,
##if so, align the two clusters
alignment.test.method.average <- function(ids1, ids2){

  compa.numbers <- get.comparison.number(ids1, ids2)
  
  ## Get the scores of the comparisons
  if (metric == "similarity"){
    scores <- compare.matrices.table[compa.numbers, names(thresholds)]
  } else if(metric == "distance"){
    scores <- compare.matrices.table[compa.numbers, names(thresholds)]
  }
  #scores <- sapply(scores, function(X){ if(X == 0){X <- NA}else {X <- X}})  
  
  ## Calculate the median of the data
  median.scores <- apply(scores, 2, mean)

  ## According to the kind of metric selected
  ## evaluates if the clusters will be aligned
  if(metric == "distance"){
    if(sum(median.scores <= unlist(thresholds)) == length(unlist(thresholds))){
      return(1)
    }else{
      return(0)
    }
  } else if (metric == "similarity"){
    if(sum(median.scores >= unlist(thresholds)) == length(unlist(thresholds))){
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

  compa.numbers <- get.comparison.number(ids1, ids2)[1]
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


########################################
## Call the program convert-matrix to
## add empty columns to the matrices  
add.empty.columns <- function(id){

  strand <- motifs.info[[id]][["strand"]]

  if(strand == "D"){
    system(paste(dir.rsat, "/perl-scripts/convert-matrix -i ", single.mat.files[[id]], " -from tf -to tf -logo_format png -return counts,consensus,parameters -insert_col_left ", merge.consensus.info[[id]][["spacer"]], " -insert_col_right ", merge.consensus.info[[id]][["offset_down"]], " -o ", cluster.folder, "/merged_consensuses/merge_level_", merge.level, "/merged_consensus_", id, ".tf", sep = ""))
  } else{

    ## First convert the matrix to reverse complement
    system(paste(dir.rsat, "/perl-scripts/convert-matrix -i ", single.mat.files[[id]], " -from tf -to tf -return counts,consensus -rc -o ", cluster.folder, "/merged_consensuses/merge_level_", merge.level, "/merged_consensus_", id, "_temp.tf", sep = ""))

    temp.mat <- paste(cluster.folder, "/merged_consensuses/merge_level_", merge.level, "/merged_consensus_", id, "_temp.tf", sep = "")

    ## Then add the gaps
    system(paste(dir.rsat, "/perl-scripts/convert-matrix -i ", temp.mat, " -from tf -to tf -logo_format png -return counts,consensus,parameters -insert_col_left ", merge.consensus.info[[id]][["spacer"]], " -insert_col_right ", merge.consensus.info[[id]][["offset_down"]], " -o ", cluster.folder, "/merged_consensuses/merge_level_", merge.level, "/merged_consensus_", id, ".tf", sep = ""))

    system(paste("rm ", temp.mat, sep = ""))
    rm(temp.mat)    
  }
}

#####################################
## Generate the aligned matrices
## for each merge level of the tree
aligned.matrices.to.merge <- function(level){

  ## Create the folder with the merged consensuses
  system(paste("mkdir -p ", cluster.folder, "/merged_consensuses/merge_level_", level, sep = ""))
  flag <- system(paste("ls ", cluster.folder, "/merged_consensuses/merge_level_", level, "/ | wc -l", sep = ""), intern = TRUE)
  if(flag >= 1){
    system(paste("rm -r ", cluster.folder, "/merged_consensuses/merge_level_", level, "/*", sep = ""))
  }
  ids <- get.id(merge.levels.leaves[[level]])
  
  ## Add the spacer to the consensuses
  merge.consensus.info <<- consensus.internal.merge(motifs.info, ids)

  ## Get the single matrices file names
  single.mat.files <<- sapply(ids, function(X){
    system(paste("ls ", out.prefix, "* | grep ", X, "| grep -v 'merged' | grep '.tf'", sep = ""), intern = TRUE)
  })
  single.mat.files <<- as.list(single.mat.files)

  ## For each level, add the empty columns to the
  ## corresponding matrices
  sapply(ids, add.empty.columns)
}


###################################################
## Set the sizes of all consensuses to the largest
## one, adding "-" at the end of the consensuses
consensus.internal.merge <- function(motifs.list, ids){

  ## Temporal list with ids info
  temporal <- list()
  for(id in ids){
    temporal[[id]] <- motifs.list[[id]]
  }

  ## Get the highest length among the consensuses
  width <- sapply(temporal, function(X){
    consensuses <- sapply(X[["consensus"]], nchar)
  })
  max.width <- unique(max(width))

  ## Fill the spaces
  for(id in ids){
    spacer.length <-  max.width - nchar(temporal[[id]][["consensus"]])
    spacer <- paste(rep("-", times = spacer.length), collapse = "")
    temporal[[id]][["consensus"]] <- paste(temporal[[id]][["consensus"]], spacer, sep = "")
    temporal[[id]][["offset_down"]] <- spacer.length
  }
  return(temporal)
}


########################################################
## Creates a temporal file with the information of
## the clusters on each partition of the JSON tree
## this information is used to create the logos tree
## with JavaScript, at the end this file is deleted
JSON.clusters <- function(){
  tree2 <<- tree
  tree2$labels <- NULL
  tree2$labels <- as.vector(global.description.table$n)
  halfway.tree2 <- hclustToTree(tree2)
  jsonTree2 <- toJSON(halfway.tree2)
  jsonTree2 <- gsub("\\],", "\\]", jsonTree2, perl = TRUE)
  jsonTree2 <- paste("{\n\"name\": \"\",\n\"children\":", jsonTree2, "}", sep = "")
  jsonTree2 <- gsub("\n\"order\":\\s+\\d+", "", jsonTree2, perl = TRUE)
  copy.Json <- jsonTree2
  copy.Json <- gsub("[^ \\d+ \\[ \\]  \\.]", "", copy.Json, perl = TRUE)
  copy.Json <- gsub("(\\d+)", "|\\1", copy.Json, perl = TRUE)
  
  copy.Json.vector <- unlist(strsplit(copy.Json, " "))
  copy.Json.vector <- copy.Json.vector[2:length(copy.Json.vector)]
  copy.Json.vector <- copy.Json.vector[which(copy.Json.vector != "")]
  
  symbol.counter <- length(which(copy.Json.vector == "[")) 
  symbol <- which(copy.Json.vector == "[")[2:symbol.counter]
  
  up.pos <- 1
  col1 <- NULL
  col2 <- NULL
  for(j in 1:(symbol.counter-1)){
    
    col1 <- append(col1, j)
    
    up.count <- 0
    down.count <- 0
    for(i in symbol[j]:length(copy.Json.vector)){
      
      ## Counters
      if(copy.Json.vector[i] == "["){
        up.count <- up.count + 1
        
        ## Positions up
        if(up.count == 1){
          up.pos <- i
        }
      } else if(copy.Json.vector[i] == "]"){
        down.count <- down.count + 1
      }
      
      ## Condition
      if(up.count == down.count){
        
        ## Position down
        down.pos <- i
        
        ## Print the number of the clusters members
                                        #print(paste(" ### ", up.pos, "   ", down.pos, " ###"))
        temp <- NULL
        temp <- gsub("[^\\d+ \\|]", "" ,copy.Json.vector[up.pos:down.pos], perl = TRUE)
        numb <- NULL
        numb <- as.integer(unlist((strsplit(paste(temp, collapse = ""), "\\|"))))
        numb <- numb[which(numb != "NA")]
                                        #print(paste(" ### ", paste(numb, collapse = " "), " ###"))
        col2 <- append(col2, paste(numb, collapse = " "))
        break
      }
    }
  }
  JSON.clusters.table <- data.frame(col1)
  JSON.clusters.table$cluster <- col2
  colnames(JSON.clusters.table) <- c(";level", "cluster")
  
  cluster <- NULL
  for(x in 1:nrow(JSON.clusters.table)){
    leaves.JSON <- sort(as.integer(unlist(strsplit(JSON.clusters.table[x,2], " "))))
    for(y in 1:length(merge.levels.leaves)){
      leaves.merge <- sort(merge.levels.leaves[[y]])
      if(length(leaves.JSON) == length(leaves.merge)){
        if(sum(leaves.merge == leaves.JSON) == length(leaves.JSON)){
          cluster <- append(cluster, paste("merge_level_", y, sep = ""))
          break
        }
      }
    } 
  }
  JSON.clusters.table$mergelevel <- cluster
  JSON.clusters.table.file <- paste(sep="", cluster.folder, "/levels_JSON_cluster_", cluster.nb,"_table.tab")
  write.table(JSON.clusters.table, file = JSON.clusters.table.file, sep = "\t", quote = FALSE, row.names = FALSE)
}


##########################################################
## Return the name of the cluster, used to define
## the clusters with the bottom-up approach
return.cluster.name <- function(nodes, clusters.list){

  matches <- sapply(clusters.list, function(X){
    sum(nodes %in% X)
  })
  
  if(sum(matches[matches == length(nodes)]) == 0){
    return(0)
  } else{
    matches <- matches[matches == length(nodes)]
    return(names(matches))
  }
}


###########################################
## Get the nodes number of a given level
nodes.by.level <- function(level.nb){
  level.nodes.string <- paste(internal.nodes.attributes[[paste("merge_level_", level.nb, sep = "")]][["cluster_1"]], internal.nodes.attributes[[paste("merge_level_", level.nb, sep = "")]][["cluster_2"]], collapse = "")
  level.nodes.numbers <- as.numeric(unlist(strsplit(level.nodes.string, " ")))
  return(level.nodes.numbers)
}



###################################################
## Fill internal nodes attributes list. Traverse
## the tree and evaluates the nodes in the branches
## according to the clustering method. Flag = 1 if
## the branch are aligned, Flag = 0, if don't.
## NOTE: this is not the alignment step
## After this step, the clusters are defined
fill.internal.nodes.attributes <- function(){

  for (merge.level in 1:nrow(tree$merge)) {
    merge.level <<- merge.level 
  #for (merge.level in 1:2) {
    print(paste("#### merge level: ", merge.level, " ####"))
    
    child1 <- tree$merge[merge.level,1]
    child2 <- tree$merge[merge.level,2]
    
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merge_level"]] <<- merge.level
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["method"]] <<- hclust.method
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["node_1"]] <<- child1
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["node_2"]] <<- child2
    
    
    ########################################
    ## Case 1: merging between two leaves ##
    ########################################
    if ((child1 < 0) && (child2 < 0)) {

      align.two.leaves(child1, child2)

      ## Identify the nodes
      n1 <- min(-child1,-child2) 
      n2 <- max(-child1,-child2)
      
      internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merging_type"]] <<- 1
    }
    
    ############################################
    ## Case 2: merging a motif with a cluster ##
    ############################################  
    if(((child1 < 0) && (child2 > 0)) || ((child1 > 0) && (child2 < 0))){

      align.leave.and.cluster(child1, child2)

      ## Identified the nodes
      n1 <- abs(min(child1, child2))
      n.aligned <- merge.levels.leaves[[merge.level]][which(merge.levels.leaves[[merge.level]] != n1)]
      n2 <- n.aligned
      
      internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merging_type"]] <<- 2
    }
    
    
    ##########################################
    ## Case 3: merging between two clusters ##
    ########################################## 
    if ((child1 > 0) && (child2 > 0)) {
      align.clusters(child1, child2)
      
      ## Identified the nodes
      N1 <- abs(min(child1, child2))
      N2 <- abs(max(child1, child2))
      
      n1 <- merge.levels.leaves[[N1]]
      n2 <- merge.levels.leaves[[N2]]
      
      internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merging_type"]] <<- 3
    }
    
    
    if(hclust.method == "single"){
      ## Get the id of each node
      central.motifs <- central.motifs.ids(get.id(n1), get.id(n2))
      id1 <- central.motifs[1]
      id2 <- central.motifs[2]
    } else if(hclust.method == "complete"){
      farthest.motifs <- farthest.motifs.ids(get.id(n1), get.id(n2))
      id1 <- farthest.motifs[1]
      id2 <- farthest.motifs[2]
    }
    
    ## Nodes attributes
    attributes <- attributes.among.clusters(get.id(n1), get.id(n2))
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["cluster_1"]] <<- paste(n1, collapse = " ")
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["cluster_2"]] <<- paste(n2, collapse = " ")
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["min_score"]] <<- attributes[1]
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["max_score"]] <<- attributes[2]
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["median_score"]] <<- attributes[3]
    
    if(hclust.method == "average"){
      aligned.motif.flag <- alignment.status(get.id(n1), get.id(n2), hclust.method)
    }else{
      aligned.motif.flag <- alignment.status(id1, id2, hclust.method)
    }
    print(paste("### INTERNAL-> FLAG: ", aligned.motif.flag, " ###", paste = ""))
    internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["flag"]] <<- aligned.motif.flag
    if(aligned.motif.flag == 0){
      internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["alignment_status"]] <<- "Non-aligned"
    }else{
      internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["alignment_status"]] <<- "Aligned"
    }
  }
}


###########################################################
## The clusters are selected by traversing the tree with
## a bottom-up approach, rather than top-down as is used
## with the R-defined  function cutree.
define.clusters.bottom.up <- function(){
  cluster.number <- 0
  ## Define the clusters: bottom-up approach
  for(level in 1:length(internal.nodes.attributes)){
     #for(level in 1:2){
    print(paste("### ", level, " ###", paste = ""))
    
    #####################################################
    ## When the merge level is between to single nodes
    ## (merging_type = 1) 
    if(as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["merging_type"]]) == 1){
      cluster.number <- cluster.number + 1
      
      if(as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]]) == 1){
        
        ## Get the numbers of the nodes
        nodes.numbers <- nodes.by.level(level)
        
        ## Store the numbers in the list
        clusters[[paste("cluster_", cluster.number, sep = "")]] <<- nodes.numbers
        
        ## Conversely, two clusters are created
      }else{
        clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(internal.nodes.attributes[[level]][["cluster_1"]])
        cluster.number <- cluster.number + 1
        clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["cluster_2"]])
        internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]] <<- 0
      }
      next
    }
    
    ##############################################
    ## When the merge level is between a single
    ## nodes and a cluster(merging_type = 2)
    if(as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["merging_type"]]) == 2){
      
      if(as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]]) == 1){
        
        ## Obtain the status of the cluster to merge
        prev.cluster.nb <- as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["node_2"]])
        prev.status <- as.numeric(internal.nodes.attributes[[paste("merge_level_", prev.cluster.nb, sep = "")]][["flag"]])
        
        if(prev.status == 1){
          cluster.name <- return.cluster.name(nodes.by.level(prev.cluster.nb), clusters)
          clusters[[cluster.name]] <<- append(clusters[[cluster.name]], as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["cluster_1"]]))
        }else{
          cluster.number <- cluster.number + 1
          clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["cluster_1"]])
          internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]] <<- 0
        }
        
      }else{
        cluster.number <- cluster.number + 1
        clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["cluster_1"]])
        internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]] <<- 0
      }
      next
    }
    
    ## ##############################################
    ## ## When the merge level is between a two 
    ## ## clusters (merging_type = 3)
    if(as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["merging_type"]]) == 3){
      
      ## Obtain the status of the two clusters before merge them
      prev.cluster.1.nb <- as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["node_1"]])
      prev.status.1 <- as.numeric(internal.nodes.attributes[[paste("merge_level_", prev.cluster.1.nb, sep = "")]][["flag"]])
      
      prev.cluster.2.nb <- as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["node_2"]])
      prev.status.2 <- as.numeric(internal.nodes.attributes[[paste("merge_level_", prev.cluster.2.nb, sep = "")]][["flag"]])
      
      if(as.numeric(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]]) == 1){
        if(prev.status.1 == 0 | prev.status.2 == 0){
          internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]] <<- 0
        } else{
          cluster.name <- return.cluster.name(nodes.by.level(prev.cluster.1.nb), clusters)
          cluster.2.nodes <- as.numeric(unlist(strsplit(internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["cluster_2"]], " ")))
          cluster.name.erase <- return.cluster.name(cluster.2.nodes, clusters)
          clusters[[cluster.name]] <<- append(clusters[[cluster.name]], cluster.2.nodes)
          clusters[[cluster.name.erase]] <<- NULL
        }
      }else{
        internal.nodes.attributes[[paste("merge_level_", level, sep = "")]][["flag"]] <<- 0
      }
      next
    }
  }  
  names(clusters) <<- paste("cluster_", seq(1:length(clusters)), sep = "")
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
