##############################################################
## Align a leaf and one cluster: align the single motif relative
## to the already aligned cluster; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.leave.and.cluster <- function(child1, child2, desc.table, compa.table, score = "Ncore", thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", hclust.tree){

  ## Identify the node numbers, the new leaf and the aligned leaves
  merge.level <- which(hclust.tree$merge == child1)
  n1 <- abs(min(child1, child2))
  n.aligned <- leaves.per.node(hclust.tree)[[merge.level]][which(leaves.per.node(hclust.tree)[[merge.level]] != n1)]
  n2 <- n.aligned
  ids.aligned <- get.id(n.aligned, desc.table)

  ## Get the id of each node on the description table
  id1.hclust <- get.id(n1, desc.table)
  id2.hclust <- get.id(n2, desc.table)

  ## Check the if the node shall be aligned
  aligned.motif.flag <- 0
  aligned.motif.flag <- alignment.test(id1.hclust, id2.hclust, compa.table, thresholds, hclust.method = method, hclust.metric = metric)

  ## In case the motifs should not be aligned
  ## fill the motifs.info list with the default parameters
  if(aligned.motif.flag == 0){

    motifs.info.temp <-  sapply(c(n1, n2), function(X){
      list(
        name = get.name(get.id(X, desc.table),desc.table),
        consensus = get.consensus(get.id(X, desc.table), desc.table, RC = FALSE),
        strand = "D",
        number = X,
        spacer.up = 0,
        spacer.dw = 0
      )
    })
    motifs.info.temp <-  data.frame(t(motifs.info.temp))
    rownames(motifs.info.temp) <- motifs.info.temp$name
    motifs.info.temp <- apply(motifs.info.temp, 1, as.list)
    motifs.info[names(motifs.info.temp)] <- motifs.info.temp[names(motifs.info.temp)]

  ## Conversely align the motifs
  } else{

#     ## Get ids
#     motifs.info[[get.id(n1)]][["number"]] <<- as.numeric(n1)
#     motifs.info[[get.id(n1)]][["spacer"]] <<- 0
#     ids.aligned <- get.id(n.aligned)

    ## Find the central motif of the cluster
    ## Get the id of each node on the description table
    central.motifs <- closest.or.farthest.motifs.ids(n1.id, n2.id, compa.table, score = "Ncor", closest = TRUE)
    id1 <- central.motifs[1]
    id2 <- central.motifs[2]

    ## NOTE: the order of the ids (ID1 or ID2) should be the same as in the comparison table
    ## If the order in the hclust tree is the opposite, then the order of the numbers is inverted
    switch.ids <- 0
    new <- NULL
    aligned <- NULL
    if(id1 != id1.hclust){
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
    rm(id1.hclust, id2.hclust)

    ## Get the comparison number in the compare-matrices table
    compa.nb <- get.comparison.number(id1, id2, compa.table)[1]

    ## Get the strand
    strand <- as.vector(compa.table[compa.nb, "strand"])

    ## Get the offset
    offset <- as.vector(compa.table[compa.nb, "offset"])

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
    ## in order to align the new leaf
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
    motifs.info[[new]][["name"]] <<- get.name(new)

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
  rm(n.aligned, n1, n2)
}
