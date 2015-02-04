##################################################
## Align two leaves: creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.two.leaves <- function(child1, child2, desc.table, compa.table, score = "Ncor", thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", hclust.tree, nodes.attributes = TRUE){

  ## Identify the node number
  n1 <- min(-child1,-child2) ## row number of the first motif in the description table
  n2 <- max(-child1,-child2) ## row number of the second motif in the description table
  merge.level <- which(hclust.tree$merge == child1)

  ## Saves the order of the IDs in the hclust object
  id1.hclust <- get.id(n1, desc.table)
  id2.hclust <- get.id(n2, desc.table)

  ## Get the id of each node on the description table
  n1.id <- get.id(n1, desc.table)
  n2.id <- get.id(n2, desc.table)

  ## Check the if the node shall be aligned
  aligned.motif.flag <- 0
  aligned.motif.flag <- alignment.test(id1.hclust, id2.hclust, compa.table, thresholds, hclust.method = method, hclust.metric = metric)

  ## Fill the attributes table
  if(nodes.attributes == TRUE){

    ## Save the merge clas: (1) two leaves, (2) one leaf and one cluster, (3) two clusters
    internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["merge_class"]] <<- 1

    ## Save the status of the alignment
    if(aligned.motif.flag == 0){
      internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["alignment_status"]] <<- "Non-aligned"
    } else{
      internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["alignment_status"]] <<- "Aligned"
    }
    internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["alignment_flag"]] <<- aligned.motif.flag

    ## Save the number of the motifs on each cluster, for each merge level
    internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["cluster_1"]] <<- paste(n1, collapse = " ")
    internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["cluster_2"]] <<- paste(n2, collapse = " ")
  }

  ## In case the motifs should not be aligned export the default parameters
  if(aligned.motif.flag == 0){

      n1.id <- get.id(n1, desc.table)
      motifs.info[[n1.id]][["name"]] <<- get.name(n1.id,desc.table)
      motifs.info[[n1.id]][["consensus_d"]] <<- get.consensus(n1.id, desc.table, RC = FALSE)
      motifs.info[[n1.id]][["consensus_rc"]] <<- get.consensus(n1.id, desc.table, RC = TRUE)
      motifs.info[[n1.id]][["strand"]] <<- "D"
      motifs.info[[n1.id]][["number"]] <<- n1
      motifs.info[[n1.id]][["spacer.up"]] <<- 0
      motifs.info[[n1.id]][["spacer.dw"]] <<- 0

      n2.id <- get.id(n2, desc.table)
      motifs.info[[n2.id]][["name"]] <<- get.name(n2.id,desc.table)
      motifs.info[[n2.id]][["consensus_d"]] <<- get.consensus(n2.id, desc.table, RC = FALSE)
      motifs.info[[n2.id]][["consensus_rc"]] <<- get.consensus(n2.id, desc.table, RC = TRUE)
      motifs.info[[n2.id]][["strand"]] <<- "D"
      motifs.info[[n2.id]][["number"]] <<- n2
      motifs.info[[n2.id]][["spacer.up"]] <<- 0
      motifs.info[[n2.id]][["spacer.dw"]] <<- 0

  ## Conversely align the motifs
  }else{

    ## Find the central motif of the cluster
    central.motifs <- closest.or.farthest.motifs.ids(n1.id, n2.id, compa.table, score = score, closest = TRUE)
    id1 <- central.motifs[1]
    id2 <- central.motifs[2]

    ## NOTE: the order of the ids (ID1 or ID2) should be the same as in the comparison table
    ## If the order in the hclust tree is the opposite, then the order of the numbers is inverted
    if(id1.hclust == id2){
      temporal <- n1
      n1 <- n2
      n2 <- temporal
    }
    rm(id1.hclust, id2.hclust, n1.id, n2.id)

    ## Comparison number in the compare-matrices table
    compa.nb <- get.comparison.number(id1, id2, compa.table)[1]

    ## Choose the relative orientation of the two motifs
    strand <- as.vector(compa.table[compa.nb, "strand"])

    consensus1a <- NULL
    consensus1b <- NULL
    consensus2a <- NULL
    consensus2b <- NULL

    ## Consensus of the first motif
    consensus1a <- get.consensus(id1, desc.table, RC = FALSE)
    consensus1b <- get.consensus(id1, desc.table, RC = TRUE)

    id1.strand <- "D"
    if (strand == "R") {
      ## Consensuses of the second motif
      consensus2a <- get.consensus(id2, desc.table, RC = TRUE)
      consensus2b <- get.consensus(id2, desc.table, RC = FALSE)
      id2.strand <- "R"
    } else {
      ## Consensuses of the second motif
      consensus2a <- get.consensus(id2, desc.table, RC = FALSE)
      consensus2b <- get.consensus(id2, desc.table, RC = TRUE)
      id2.strand <- "D"
    }

    ## Add the offset to the logos
    offset <- as.vector(compa.table[compa.nb, "offset"])
    spacer <- paste(collapse="",rep(x="-",times=abs(offset)))

    if (offset < 0) {
      consensus1a <- paste(spacer, consensus1a, sep = "")
      consensus1b <- paste(consensus1b, spacer, sep = "")
    } else {
      consensus2a <- paste(spacer, consensus2a, sep = "")
      consensus2b <- paste(consensus2b, spacer, sep = "")
    }


    ## Update the motifs information
    motifs.info[[id1]][["name"]] <<- get.name(id1, desc.table)
    motifs.info[[id1]][["consensus_d"]] <<- consensus1a
    motifs.info[[id1]][["consensus_rc"]] <<- consensus1b
    motifs.info[[id1]][["strand"]] <<- id1.strand
    motifs.info[[id1]][["number"]] <<- n1
    motifs.info[[id1]][["spacer.up"]] <<- get.spacer.nb(motifs.info[[id1]][["consensus_d"]])$up.spacer
    motifs.info[[id1]][["spacer.dw"]] <<- get.spacer.nb(motifs.info[[id1]][["consensus_d"]])$dw.spacer

    motifs.info[[id2]][["name"]] <<- get.name(id2, desc.table)
    motifs.info[[id2]][["consensus_d"]] <<- consensus2a
    motifs.info[[id2]][["consensus_rc"]] <<- consensus2b
    motifs.info[[id2]][["strand"]] <<- id2.strand
    motifs.info[[id2]][["number"]] <<- n2
    motifs.info[[id2]][["spacer.up"]] <<- get.spacer.nb(motifs.info[[id2]][["consensus_d"]])$up.spacer
    motifs.info[[id2]][["spacer.dw"]] <<- get.spacer.nb(motifs.info[[id2]][["consensus_d"]])$dw.spacer

    motifs.info.temp <- fill.downstream(get.id(leaves.per.node(hclust.tree)[[merge.level]], desc.table), motifs.info)
    motifs.info[names(motifs.info.temp)] <<- motifs.info.temp[names(motifs.info.temp)]
  }
}
