##############################################################
## Align a leaf and one cluster: align the single motif relative
## to the already aligned cluster; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.leaf.and.cluster <- function(child1,
                                   child2,
                                   thresholds = list(Ncor = 0.4, cor = 0.6, w = 5),
                                   method = "average",
                                   metric = "Ncor",
                                   align = TRUE,
                                   nodes.attributes = TRUE,
                                   motif.at.tree.level = motif.at.tree.level){

  ## Identify the node numbers: the new node leaf and the aligned leaves
  merge.level <- which(tree$merge == child1)
  n1 <- abs(min(child1, child2))
  n.aligned <- leaves.per.node(tree)[[merge.level]][which(leaves.per.node(tree)[[merge.level]] != n1)]
  n2 <- n.aligned
  ids.aligned <- get.id(n.aligned)

  ## Get the id of each node
  id1.hclust <- get.id(n1)
  id2.hclust <- get.id(n2)

  ## Initialize the information for the n1 (new) motif
  motifs.info[[id1.hclust]][["name"]] <<- n1
  motifs.info[[id1.hclust]][["spacer.up"]] <<- 0
  motifs.info[[id1.hclust]][["spacer.dw"]] <<- 0

  ## According to the hierarchical clustering method selected,
  ## Check if the motifs corresponding to the current level shall be aligned
  aligned.motif.flag <- 0
  aligned.motif.flag <- check.alignment(id1.hclust,
                                        id2.hclust,
                                        thresholds,
                                        hclust.method = method,
                                        metric = metric)

  ## Fill the attributes table
  if(nodes.attributes == TRUE){

    ## Save the merge class: (class 1) two leaves, (class 2) one leaf and one cluster, (class 3) two clusters
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["merge_class"]] <<- 2

    ## Save the status of the alignment 0 = Non-Aligned , 1 = Aligned
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["alignment_flag"]] <<- aligned.motif.flag

    ## Save the number of the motifs on each cluster, for each merge level
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["cluster_1"]] <<- paste(n1, collapse = " ")
    internal.nodes.attributes[[paste("node_", merge.level, sep = "")]][["cluster_2"]] <<- paste(n2, collapse = " ")
  }

  if(align == TRUE){

    ## In case the motifs should not be aligned
    ## fill the motifs.info list with the default parameters
    if(aligned.motif.flag == 0){

      n1.id <- get.id(n1)

      motifs.info[[n1.id]][["name"]] <<- get.name(n1.id)
      motifs.info[[n1.id]][["consensus_d"]] <<- get.consensus(n1.id, RC = FALSE)
      motifs.info[[n1.id]][["consensus_rc"]] <<- get.consensus(n1.id, RC = TRUE)
      motifs.info[[n1.id]][["strand"]] <<- "D"
      motifs.info[[n1.id]][["number"]] <<- n1
      motifs.info[[n1.id]][["spacer.up"]] <<- 0
      motifs.info[[n1.id]][["spacer.dw"]] <<- 0

      ## Conversely align the motifs
    } else{

      ## Set the number of the motifs
      ## This is done here because the variable storing this number can be changed
      motifs.info[[get.id(n1)]][["number"]] <<- as.numeric(n1)

      ## Find the central motif of the cluster
      ## Get the id of each node on the description table
      central.motifs <- closest.or.farthest.motifs.ids(id1.hclust, id2.hclust, metric = metric, closest = TRUE)
      id1 <- central.motifs[1]
      id2 <- central.motifs[2]

      ## NOTE: the order of the ids (ID1 or ID2) should be the same as in the comparison table
      ## If the order in the hclust tree is the opposite, then the order of the numbers is inverted
      switch.ids <- 0
      new.node <- NULL
      aligned <- NULL
      if(id1 != id1.hclust){
        aligned <- id1
        new.node <- id2
        temporal <- n1
        n1 <- n2
        n2 <- temporal
        switch.ids <- 1
      } else{
        aligned <- id2
        new.node <- id1
      }
      rm(id1.hclust, id2.hclust)

      ## Get the comparison number in the compare-matrices table
      compa.nb <- get.comparison.number(id1, id2)[1]

      ## Get the strand
      strand <- as.vector(global.compare.matrices.table[compa.nb, "strand"])

      ## Get the offset
      offset <- as.vector(global.compare.matrices.table[compa.nb, "offset"])

      ## Assign values for the cases (1-8)
      ## Each case depends if the Ids were inverted, if the strand of the comparison between the closest motifs
      ## is R or D, and in the orientation of the previously aligned motif.
      case <- 0
      if(switch.ids == 1){
        if(strand == "D"){
          if(motifs.info[[aligned]][["strand"]] == "D"){

            ## Case 1: inverted IDs, current comparison strand = 'D' and strand of previously aligned motif = 'D'
            case <- 1
          } else{
            ## Case 2: inverted IDs, current comparison strand = 'D' and strand of previously aligned motif = 'R'
            case <- 2
          }
        } else{
          if(motifs.info[[aligned]][["strand"]] == "D"){

            ## Case 3: inverted IDs, current comparison strand = 'R' and strand of previously aligned motif = 'D'
            case <- 3
          } else{
            ## Case 4: inverted IDs, current comparison strand = 'R' and strand of previously aligned motif = 'R'
            case <- 4
          }
        }
      }else{
        if(strand == "D"){
          if(motifs.info[[aligned]][["strand"]] == "D"){

            ## Case 5: non-inverted IDs, current comparison strand = 'D' and strand of previously aligned motif = 'D'
            case <- 5
          } else{
            ## Case 6: non-inverted IDs, current comparison strand = 'D' and strand of previously aligned motif = 'R'
            case <- 6
          }
        } else{
          if(motifs.info[[aligned]][["strand"]] == "D"){

            ## Case 7: non-inverted IDs, current comparison strand = 'R' and strand of previously aligned motif = 'D'
            case <- 7
          } else{
            ## Case 8: non-inverted IDs, current comparison strand = 'R' and strand of previously aligned motif = 'R'
            case <- 8
          }
        }
      }

      ## In some cases it is required to invert the alignment (change its orientation)
      ## in order to align the new.node leaf and store in the motifs.info variable
      if(case %in% c(2,4,6,7)){
        inverted.alignment.aligned.ids <- inverted.alignment(ids.aligned, motifs.info)
        motifs.info[names(inverted.alignment.aligned.ids)] <<- inverted.alignment.aligned.ids[names(inverted.alignment.aligned.ids)]
      }

      ## Choose the consensus for the new.node motif
      if(case %in% c(1,2,5,6,7,8)){
        consensus.1.new.node <- get.consensus(new.node, RC = TRUE)
        consensus.2.new.node <- get.consensus(new.node, RC = FALSE)
        motifs.info[[new.node]][["strand"]] <<- "D"
        motifs.info[[new.node]][["consensus_d"]] <<- consensus.1.new.node
        motifs.info[[new.node]][["consensus_rc"]] <<- consensus.2.new.node

      } else if(case %in% c(3,4)){
        consensus.1.new.node <- get.consensus(new.node, RC = TRUE)
        consensus.2.new.node <- get.consensus(new.node, RC = FALSE)
        motifs.info[[new.node]][["strand"]] <<- "R"
        motifs.info[[new.node]][["consensus_d"]] <<- consensus.1.new.node
        motifs.info[[new.node]][["consensus_rc"]] <<- consensus.2.new.node
      }
      motifs.info[[new.node]][["name"]] <<- get.name(new.node)

      ## Reset the offset
      aligned.spacer.up <- as.numeric(motifs.info[[aligned]][["spacer.up"]])
      if(case %in% c(1:4)){
        spacer.diff <- aligned.spacer.up
      } else if(case %in% c(5:8)){
        spacer.diff <- -aligned.spacer.up
      }
      offset <- offset + spacer.diff

      ## Create the spacer
      spacer <- paste(collapse="",rep(x="-",times = abs(offset)))

      ## Add the gaps
      ids.mod <- NULL
      if(offset <= 0){
        ids.mod <- get.id(n1)
      }else{
        ids.mod <- get.id(n2)
      }

      temp.motifs.info <- motifs.info[ids.mod]
      sapply(ids.mod, function(x){
        temp <- list()
        temp[[x]][["strand"]] <- temp.motifs.info[[x]][["strand"]]
        temp[[x]][["consensus_d"]] <- paste(spacer, temp.motifs.info[[x]][["consensus_d"]], sep="")
        temp[[x]][["consensus_rc"]] <- paste(temp.motifs.info[[x]][["consensus_rc"]], spacer, sep="")
        temp[[x]][["name"]] <- temp.motifs.info[[x]][["name"]]
        temp[[x]][["number"]] <- as.numeric(temp.motifs.info[[x]][["number"]])
        temp[[x]][["spacer.up"]] <- get.spacer.nb(temp.motifs.info[[x]][["consensus_d"]])$up.spacer
        temp[[x]][["spacer.dw"]] <- get.spacer.nb(temp.motifs.info[[x]][["consensus_d"]])$dw.spacer
        temp.motifs.info[names(temp)] <<- temp[names(temp)]
      })

      ## Fill the downstream end
      motifs.info[names(temp.motifs.info)] <<- temp.motifs.info[names(temp.motifs.info)]
      temp.motifs.info <- fill.downstream(get.id(motif.at.tree.level[[merge.level]]), motifs.info)
      motifs.info[names(temp.motifs.info)] <<- temp.motifs.info[names(temp.motifs.info)]
    }
  }
}