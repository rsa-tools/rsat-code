##################################################
## Align two leaves: creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.two.leaves <- function(child1,
                             child2,
                             thresholds = list(Ncor = 0.4, cor = 0.6, w = 5),
                             hclust.method = "average",
                             metric = "Ncor",
                             align = TRUE,
                             nodes.attributes = TRUE,
                             motif.at.tree.level = motif.at.tree.level){


  ## Identify the node numbers and merge level
  n1 <- min(-child1,-child2)
  n2 <- max(-child1,-child2)
  merge.level <- which(tree$merge == child1)

  ## Saves the order of the IDs in the hclust$merge object
  id1.hclust <- get.id(n1)
  id2.hclust <- get.id(n2)

  # message("; Aligning: ", id1.hclust, " - ", id2.hclust)
  if(align == TRUE){
    verbose(paste("Aligning motifs: ", id1.hclust, " - ", id2.hclust), 1)
  }


  ## According to the hierarchical clustering method selected,
  ## Check if the motifs corresponding to the current level shall be aligned
  aligned.motif.flag <- 0
  aligned.motif.flag <- check.alignment(id1.hclust,
                                        id2.hclust,
                                        thresholds,
                                        hclust.method = hclust.method,
                                        metric = metric)



  ## Fill the attributes table
  if(nodes.attributes == TRUE){

    ## Save the merge class: (class 1) two leaves, (class 2) one leaf and one cluster, (class 3) two clusters
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["merge_class"]] <<- 1

    ## Save the status of the alignment 0 = Non-Aligned , 1 = Aligned
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["alignment_flag"]] <<- aligned.motif.flag

    ## Save the number of the motifs on each cluster, for each merge level
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["cluster_1"]] <<- paste(n1, collapse = " ")
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["cluster_2"]] <<- paste(n2, collapse = " ")
  }

  if(align == TRUE){

    ## In case the motifs should not be aligned export the default parameters
    if(aligned.motif.flag == 0){

        for(n in c(n1,n2)){
          n.id <- get.id(n)
          motifs.info[[n.id]][["name"]] <<- get.name(n.id)
          motifs.info[[n.id]][["consensus_d"]] <<- get.consensus(n, RC = FALSE)
          motifs.info[[n.id]][["consensus_rc"]] <<- get.consensus(n, RC = TRUE)
          motifs.info[[n.id]][["strand"]] <<- "D"
          motifs.info[[n.id]][["number"]] <<- n
          motifs.info[[n.id]][["spacer.up"]] <<- 0
          motifs.info[[n.id]][["spacer.dw"]] <<- 0
        }

    ## Conversely align the motifs
    }else{

      ## Find the central pair of motifs of the cluster
      central.motifs <- closest.or.farthest.motifs.ids(id1.hclust,
                                                       id2.hclust,
                                                       metric = metric,
                                                       closest = TRUE)
      id1 <- central.motifs[1]
      id2 <- central.motifs[2]

      ## NOTE: the order of the ids (ID1 or ID2) should be the same as in the comparison table
      ## If the order in the hclust tree is the opposite, then the order of the numbers is inverted
      if(id1.hclust == id2){
        temporal <- n1
        n1 <- n2
        n2 <- temporal
      }
      rm(id1.hclust, id2.hclust)

      ## Comparison number in the compare-matrices table
      compa.nb <- get.comparison.number(id1, id2)[1]

      ## Choose the relative orientation of the two motifs
      strand <- as.vector(global.compare.matrices.table[compa.nb, "strand"])

      consensus1a <- NULL
      consensus1b <- NULL
      consensus2a <- NULL
      consensus2b <- NULL

      ## Consensus of the first motif
      consensus1a <- get.consensus(id1, RC = FALSE)
      consensus1b <- get.consensus(id1, RC = TRUE)

      id1.strand <- "D"
      if (strand == "R") {
        ## Consensuses of the second motif
        consensus2a <- get.consensus(id2, RC = TRUE)
        consensus2b <- get.consensus(id2, RC = FALSE)
        id2.strand <- "R"
      } else {
        ## Consensuses of the second motif
        consensus2a <- get.consensus(id2, RC = FALSE)
        consensus2b <- get.consensus(id2, RC = TRUE)
        id2.strand <- "D"
      }

      ## Add the offset to the consensuses
      offset <- as.vector(global.compare.matrices.table[compa.nb, "offset"])
      spacer <- paste(collapse="",rep(x="-",times=abs(offset)))

      if (offset < 0) {
        consensus1a <- paste(spacer, consensus1a, sep = "")
        consensus1b <- paste(consensus1b, spacer, sep = "")
      } else {
        consensus2a <- paste(spacer, consensus2a, sep = "")
        consensus2b <- paste(consensus2b, spacer, sep = "")
      }

      ## Update the motifs information
      motifs.info[[id1]][["name"]] <<- get.name(id1)
      motifs.info[[id1]][["consensus_d"]] <<- consensus1a
      motifs.info[[id1]][["consensus_rc"]] <<- consensus1b
      motifs.info[[id1]][["strand"]] <<- id1.strand
      motifs.info[[id1]][["number"]] <<- n1
      motifs.info[[id1]][["spacer.up"]] <<- get.spacer.nb(motifs.info[[id1]][["consensus_d"]])$up.spacer
      motifs.info[[id1]][["spacer.dw"]] <<- get.spacer.nb(motifs.info[[id1]][["consensus_d"]])$dw.spacer
      motifs.info[[id2]][["name"]] <<- get.name(id2)
      motifs.info[[id2]][["consensus_d"]] <<- consensus2a
      motifs.info[[id2]][["consensus_rc"]] <<- consensus2b
      motifs.info[[id2]][["strand"]] <<- id2.strand
      motifs.info[[id2]][["number"]] <<- n2



      motifs.info[[id2]][["spacer.up"]] <<- get.spacer.nb(motifs.info[[id2]][["consensus_d"]])$up.spacer
      motifs.info[[id2]][["spacer.dw"]] <<- get.spacer.nb(motifs.info[[id2]][["consensus_d"]])$dw.spacer

      motifs.info.temp <- fill.downstream(get.id(motif.at.tree.level[[merge.level]]), motifs.info)
      motifs.info[names(motifs.info.temp)] <<- motifs.info.temp[names(motifs.info.temp)]
    }
  }
}