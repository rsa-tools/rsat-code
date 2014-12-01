##################################################
## Align two leaves: creates a list with the info
## (strand, consensus, offset) of the aligned motifs
align.two.leaves <- function(child1, child2, desc.table, compa.table, score = "Ncore", thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", hclust.tree){

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
  aligned.motif.flag <- alignment.test(id1, id2, compa.table, thresholds, hclust.method = method, hclust.metric = metric)

  ## In case the motifs should not be aligned export the default parameters
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
  }else{

    ## Find the central motif of the cluster
    central.motifs <- closest.or.farthest.motifs.ids(n1.id, n2.id, compa.table, score = "Ncor", closest = TRUE)
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
    consensus1 <- as.vector(desc.table[n1,"consensus"]) ## Consensus of the first motif
    id1.strand <- "D"
    #motifs.info[[id1]][["strand"]] <<- "D"
    if (strand == "R") {
      consensus2 <- as.vector(desc.table[n2,"rc.consensus"]) ## Consensus of the second motif
      id2.strand <- "R"
      #motifs.info[[id2]][["strand"]] <<- "R"
    } else {
      consensus2 <- as.vector(desc.table[n2,"consensus"]) ## Consensus of the second motif
      id2.strand <- "R"
      #motifs.info[[id2]][["strand"]] <<- "D"
    }

    ## Add the offset to the logos
    offset <- as.vector(compa.table[compa.nb, "offset"])
    spacer <- paste(collapse="",rep(x="-",times=abs(offset)))

    if (offset < 0) {
      consensus1 <- paste(sep="", spacer, consensus1)
    } else {
      consensus2 <- paste(sep="", spacer, consensus2)
    }

    ## Reset the consensus with the aligned and re-oriented consensus
#     tree$labels[n1] <<- paste(consensus1, n1)
#     tree$labels[n2] <<- paste(consensus2, n2)

    ## Update the motifs information
    motifs.info[[id1]][["name"]] <<- get.name(id1, desc.table)
    motifs.info[[id1]][["consensus"]] <<- consensus1
    motifs.info[[id1]][["strand"]] <<- id1.strand
    motifs.info[[id1]][["number"]] <<- n1
    motifs.info[[id1]][["spacer.up"]] <<- length(unlist(strsplit(motifs.info[[id1]][["consensus"]], "-")))-1
    motifs.info[[id1]][["spacer.dw"]] <<- 0

    motifs.info[[id2]][["name"]] <<- get.name(id2, desc.table)
    motifs.info[[id2]][["consensus"]] <<- consensus2
    motifs.info[[id2]][["strand"]] <<- id2.strand
    motifs.info[[id2]][["number"]] <<- n2
    motifs.info[[id2]][["spacer.up"]] <<- length(unlist(strsplit(motifs.info[[id2]][["consensus"]], "-")))-1
    motifs.info[[id2]][["spacer.dw"]] <<- 0

    motifs.info.temp <- fill.downstream(get.id(leaves.per.level[[merge.level]], desc.table), motifs.info)
    motifs.info[names(motifs.info.temp)] <<- motifs.info.temp[names(motifs.info.temp)]
  }
}
