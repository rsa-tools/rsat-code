#################################################################
## Define the clusters: Using a bottom-up approach, and using
## the order of assembly of the hierarchcial tree and the
## attributes information with the alignment status at each
## node of the tree.
find.clusters <- function(attributes.list, tree){

  #######################################################################
  ## Given a level of one hierarchical tree and the alignment status
  ## at each level of the tree, search the levels thar are chained.
  ## This means, search the chain of levels that were grouped together
  find.chained.levels <- function(x){

    chained.levels <- NULL
    current.level <- NULL
    alignment.flag <- 1
    chained.levels <- append(chained.levels, x)

    while(alignment.flag == 1){

      ## Search the next chained level
      current.level <- find.next.levels.in.tree(x)

      ## The case when the root tree is analyzed
      if(x == current.level){
        return(x)
      }

      ## Get the alignment status at the level
      alignment.flag <- as.numeric(attributes.list[[paste("node_", current.level, sep = "")]][["alignment_flag"]])

      ## If the status at the current level is 0 then return the chained levels
      if(alignment.flag == 0){
        return(chained.levels)

        ## Conversely, add the new level to the chain
      } else {
        chained.levels <- append(chained.levels, current.level)
        x <- current.level
      }
    }
  }


  ################################################
  ## Given a level of a hierarchical tree, find
  ## the next level pointing the current level
  ## NOTE: this function is only called within the function find.clusters
  find.next.levels.in.tree <- function(x){

    if (x == length(tree$merge)/2){
      return(x)

    } else {

      ## Get the level
      level <- which(tree$merge == x)

      if(level > length(tree$merge)/2 & level <= length(tree$merge)){
        return(level - length(tree$merge)/2)
      } else if(level <= length(tree$merge)/2){
        return(level)
      }
    }
  }


  ## Initialize the variables
  clusters <<- list()
  cluster.counter <<- 0
  motifs.at.tree.level <<- leaves.per.node(tree)
  checked.levels <<- rep(0, times = length(attributes.list))

  ## Attribute treatment. If a particular node with an alignment flag == 0
  ## is pointing to a node that potentially could be aligned (flag == 1)
  ## in a further step the last is set to 0, in order to avoid wrong alignments.
  if(length(attributes.list) > 2){
    sapply(1:(length(attributes.list)-1), function(x){

      ## Get the flag of the current level
      flag <- as.numeric(attributes.list[[x]][["alignment_flag"]])

      if(flag == 0){
        next.level <- find.next.levels.in.tree(x)
        next.level.flag <- as.numeric(attributes.list[[next.level]][["alignment_flag"]])
        if(next.level.flag == 1){
          attributes.list[[next.level]][["alignment_flag"]] <<- 0
        }
      }
    })
  }


  #################################################
  ## The three is traversed in a bottom-up way
  sapply(1:length(attributes.list), function(lvl){

    ## Only those unchecked levels will be analyzed
    if(checked.levels[lvl] == 0){

      ## Get the alignment status of the current level
      level.align.flag <- attributes.list[[paste("node_", lvl, sep = "")]][["alignment_flag"]]

      ## Skip to next level if the alignment flag is 0
      if(level.align.flag == 0){

        ## If the alignment flag is 0 but there is at least one single node
        ## Each single node will correspond to a new one cluster
        if( sum(tree$merge[lvl,] < 1) > 0 ){

          for (m in tree$merge[lvl,]){
            if(m < 1){
              cluster.counter <<- cluster.counter + 1
              clusters[[paste("cluster_", cluster.counter, sep = "")]] <<- abs(m)
            }
          }
        }

        ## Check the current level
        checked.levels[lvl] <<- 1

        ## Find the chained levels corresponding to the current level
      } else {

        chained.lvl <- find.chained.levels(lvl)

        ## Those levels corresponding to the chained levels are checked
        for (l in chained.lvl){
          checked.levels[l] <<- 1
        }

        ## Get the highest levels where the motifs were grouped
        max.chained.lvl <- max(chained.lvl)

        ## Get the motifs IDs corresponding to this level
        chained.nb <- motifs.at.tree.level[[max.chained.lvl]]

        ## Increase the number of clusters and store the new one
        cluster.counter <<- cluster.counter + 1
        clusters[[paste("cluster_", cluster.counter, sep = "")]] <<- chained.nb
      }
    }
  })

  ## Obtain the unique elements from the clusters gathered
  ## and rename them
  clusters <- unique(clusters)
  names(clusters) <- paste("cluster", 1:length(clusters), sep = "_")
  return(clusters)
}