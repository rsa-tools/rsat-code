##############################################################
## Align a leaf and one cluster: align the single motif relative
## to the already aligned cluster; creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.leaf.and.cluster <- function(child1, child2, desc.table, compa.table, score = "Ncore", thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", hclust.tree, nodes.attributes = TRUE){

  print("Case 2")

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

  ## Fill the attributes table
  if(nodes.attributes == TRUE){

    ## Save the merge clas: (1) two leaves, (2) one leaf and one cluster, (3) two clusters
    internal.nodes.attributes[[paste("level_", merge.level, sep = "")]][["merge_class"]] <<- 2

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
    motifs.info[names(motifs.info.temp)] <<- motifs.info.temp[names(motifs.info.temp)]

  ## Conversely align the motifs
  } else{

    ## Set the number of the motifs
    ## This is done here because the variable storing this number can be changed
    motifs.info[[get.id(n1, desc.table)]][["number"]] <<- as.numeric(n1)

    ## Find the central motif of the cluster
    ## Get the id of each node on the description table
    central.motifs <- closest.or.farthest.motifs.ids(id1.hclust, id2.hclust, compa.table, score = "Ncor", closest = TRUE)
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
    ## in order to align the new leaf and store in the motifs.info variable
    if(case %in% c(2,4,6,7)){
      inverted.alignment.aligned.ids <- inverted.alignment(ids.aligned, motifs.info, desc.table)
      motifs.info[names(inverted.alignment.aligned.ids)] <<- inverted.alignment.aligned.ids[names(inverted.alignment.aligned.ids)]
    }

    ## Choose the consensus for the new motif
    if(case %in% c(1,2,5,6,7,8)){
      consensus.new <- get.consensus(new, desc.table, RC = FALSE)
      motifs.info[[new]][["strand"]] <<- "D"
      motifs.info[[new]][["consensus"]] <<- consensus.new

    } else if(case %in% c(3,4)){
      consensus.new <- get.consensus(new, desc.table, RC = TRUE)
      motifs.info[[new]][["strand"]] <<- "R"
      motifs.info[[new]][["consensus"]] <<- consensus.new
    }
    motifs.info[[new]][["name"]] <<- get.name(new, desc.table)

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
    if(offset <= 0){
        ids.mod <- get.id(n1, desc.table)
    }else{
      ids.mod <- get.id(n2, desc.table)
    }
    temp.motifs.info <- motifs.info
    sapply(ids.mod, function(id){
      temp <- list()
      temp[[id]][["strand"]] <- temp.motifs.info[[id]][["strand"]]
      temp[[id]][["consensus"]] <- paste(spacer, temp.motifs.info[[id]][["consensus"]], sep="")
      temp[[id]][["name"]] <- temp.motifs.info[[id]][["name"]]
      temp[[id]][["number"]] <- as.numeric(temp.motifs.info[[id]][["number"]])
      temp[[id]][["spacer.up"]] <- get.spacer.nb(motifs.info[[id]][["consensus"]])$up.spacer
      temp[[id]][["spacer.dw"]] <- get.spacer.nb(motifs.info[[id]][["consensus"]])$dw.spacer
      temp.motifs.info[names(temp)] <<- temp[names(temp)]
    })

    ## Fill the downstream end
    motifs.info[names(temp.motifs.info)] <<- temp.motifs.info[names(temp.motifs.info)]
    temp.motifs.info <- fill.downstream(get.id(leaves.per.node(tree)[[merge.level]], desc.table), motifs.info)
    motifs.info[names(temp.motifs.info)] <<- temp.motifs.info[names(temp.motifs.info)]
  }
}