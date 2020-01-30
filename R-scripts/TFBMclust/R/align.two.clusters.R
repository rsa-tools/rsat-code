##############################################################
## Align two clusters: align the two cluster, in some cases is
## necessary invert the alignment; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.two.clusters <- function(child1,
                               child2,
                               thresholds = list(Ncor = 0.4, cor = 0.6, w = 5),
                               method = "average",
                               metric = "Ncor",
                               align = TRUE,
                               nodes.attributes = TRUE,
                               motif.at.tree.level = motif.at.tree.level){

  ## Identify the merge level and the node numbers of each cluster
  merge.level <- which(tree$merge == child1)
  N1 <- abs(min(child1, child2))
  N2 <- abs(max(child1, child2))
  n1 <- motif.at.tree.level[[N1]]
  n2 <- motif.at.tree.level[[N2]]

  ## Get the id of each node on the description table
  ids1.hclust <- get.id(n1)
  ids2.hclust <- get.id(n2)

  ## Check the if the node shall be aligned
  aligned.motif.flag <- 0
  aligned.motif.flag <- check.alignment(ids1.hclust,
                                        ids2.hclust,
                                        thresholds,
                                        hclust.method = method,
                                        metric = metric)

  ## Fill the attributes table
  if(nodes.attributes == TRUE){

    ## Save the merge clas: (1) two leaves, (2) one leaf and one cluster, (3) two clusters
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["merge_class"]] <<- 3
    ## Save the status of the alignment 0 = Non-Aligned , 1 = Aligned
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["alignment_flag"]] <<- aligned.motif.flag

    ## Save the number of the motifs on each cluster, for each merge level
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["cluster_1"]] <<- paste(n1, collapse = " ")
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["cluster_2"]] <<- paste(n2, collapse = " ")
  }

  if(align == TRUE){

    ## Align the motifs
    if(aligned.motif.flag == 1){

      ## Find the central motif of the cluster
      ## Get the id of each node on the description table
      central.motifs <- closest.or.farthest.motifs.ids(ids1.hclust, ids2.hclust, metric = metric, closest = TRUE)
      id1 <- central.motifs[1]
      id2 <- central.motifs[2]

      ## NOTE: the order of the ids (ID1 or ID2) should be the same as in the comparison table
      ## If the order in the hclust tree is the opposite, then the order of the numbers is inverted
      switch.ids <- 0
      if(id1 %in% ids2.hclust){
        temporal <- NULL
        temporal <- ids1.hclust
        ids1.hclust <- ids2.hclust
        ids2.hclust <- temporal
        switch.ids <- 1
      }

      ## Get the comparison number in the compare-matrices table
      compa.nb <- get.comparison.number(id1, id2)[1]

      ## Get the strand
      strand <- as.vector(global.compare.matrices.table[compa.nb, "strand"])

      ## Get the offset
      offset <- as.vector(global.compare.matrices.table[compa.nb, "offset"])

      ## Get the previous orientation of the closest motifs which will be used as
      ## base to align the two cluster
      prev.strand.1 <- motifs.info[[id1]][["strand"]]
      prev.strand.2 <- motifs.info[[id2]][["strand"]]

      ## Assign values for the cases (1-8)
      ## Each case depends if the strand of the comparison between the closest motifs
      ## is R or D, and if the strand of the closest motifs is in D or D orientation.
      case <- 0
      if(strand == "D"){

        ## Case 1: current comparison strand = 'D', strand of motif 1 = 'R', strand of motif 2 = 'R'
        if(prev.strand.1 == "R" && prev.strand.2 == "R"){
          case <- 1

          ## Case 2: current comparison strand = 'D', strand of motif 1 = 'R', strand of motif 2 = 'D'
        } else if(prev.strand.1 == "R" && prev.strand.2 == "D"){
          case <- 2

          ## Case 3: current comparison strand = 'D', strand of motif 1 = 'D', strand of motif 2 = 'R'
        } else if(prev.strand.1 == "D" && prev.strand.2 == "R"){
          case <- 3

          ## Case 4: current comparison strand = 'D', strand of motif 1 = 'D', strand of motif 2 = 'D'
        } else if(prev.strand.1 == "D" && prev.strand.2 == "D"){
          case <- 4
        }
      } else{

        ## Case 5: current comparison strand = 'R', strand of motif 1 = 'R', strand of motif 2 = 'R'
        if(prev.strand.1 == "R" && prev.strand.2 == "R"){
          case <- 5

          ## Case 6: current comparison strand = 'R', strand of motif 1 = 'R', strand of motif 2 = 'D'
        } else if(prev.strand.1 == "R" && prev.strand.2 == "D"){
          case <- 6

          ## Case 7: current comparison strand = 'R', strand of motif 1 = 'D', strand of motif 2 = 'R'
        } else if(prev.strand.1 == "D" && prev.strand.2 == "R"){
          case <- 7

          ## Case 8: current comparison strand = 'R', strand of motif 1 = 'D', strand of motif 2 = 'D'
        } else if(prev.strand.1 == "D" && prev.strand.2 == "D"){
          case <- 8
        }
      }

      ## In some cases it is required to invert the alignment (change its orientation)
      ## in order to align properly the two clusters.
      ## Cases in which is required invert the aligment
      if(case %in% c(2,3,5,8)){

        if(case %in% c(3,8)){
          ids.inv <- ids2.hclust
        } else if (case %in% c(2,5)){
          ids.inv <- ids1.hclust
        }

        ## Invert the aligment and store the information in a list
        inverted.alignment.ids <- inverted.alignment(ids.inv, motifs.info)
        motifs.info[names(inverted.alignment.ids)] <<- inverted.alignment.ids[names(inverted.alignment.ids)]
      }

      ## Get the (upstream) spacer length of the closest motifs which will be used as
      ## base to align the two clusters
      cluster.1.spacer <- as.numeric(motifs.info[[id1]][["spacer.up"]])
      cluster.2.spacer <- as.numeric(motifs.info[[id2]][["spacer.up"]])

      ## According to the cases, reset the offset
      if(case %in% c(1,6)){
        offset <- nchar(get.consensus(id1, RC = FALSE)) - nchar(get.consensus(id2, RC = FALSE)) - offset + (cluster.1.spacer - cluster.2.spacer)
      } else if(case %in% c(2,3,4,5,7,8)){
        offset <- offset + (cluster.1.spacer - cluster.2.spacer)
      }

      ## Add the spacer to the motifs
      ids.mod <- NULL
      if(offset <= 0){
        ids.mod <- ids1.hclust
      } else{
        ids.mod <- ids2.hclust
      }
      spacer <- paste(collapse="",rep(x="-",times = abs(offset)))
      temp.motifs.info <- motifs.info[ids.mod]

      sapply(ids.mod, function(y){
        temp <- list()
        temp[[y]][["strand"]] <- temp.motifs.info[[y]][["strand"]]
        temp[[y]][["consensus_d"]] <- paste(spacer, temp.motifs.info[[y]][["consensus_d"]], sep = "")
        temp[[y]][["consensus_rc"]] <- paste(temp.motifs.info[[y]][["consensus_rc"]], spacer, sep = "")
        temp[[y]][["name"]] <- get.name(y)
        temp[[y]][["number"]] <- as.numeric(temp.motifs.info[[y]][["number"]])
        temp[[y]][["spacer.up"]] <- get.spacer.nb(temp.motifs.info[[y]][["consensus_d"]])$up.spacer
        temp[[y]][["spacer.dw"]] <- get.spacer.nb(temp.motifs.info[[y]][["consensus_d"]])$dw.spacer

        ## Update the motifs list
        temp.motifs.info[names(temp)] <<- temp[names(temp)]
      })

      ## Fill the downstream end
      motifs.info[names(temp.motifs.info)] <<- temp.motifs.info[names(temp.motifs.info)]
      temp.motifs.info <- fill.downstream(get.id(motif.at.tree.level[[merge.level]]), motifs.info)
      motifs.info[names(temp.motifs.info)] <<- temp.motifs.info[names(temp.motifs.info)]

    }
  }
}