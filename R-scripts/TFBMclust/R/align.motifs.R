align.motifs <- function(thresholds = list(Ncor = 0.4, cor = 0.6, w = 5),
                         method = "average",
                         metric = "Ncor",
                         align = FALSE,
                         nodes.attributes = TRUE,
                         intermediate.alignments = FALSE){


  ################################################
  ## Given a level of a hierarchical tree, find
  ## the next level pointing the current level
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


  motifs.info <<- list()
  motifs.info.levels <- list()
  internal.nodes.attributes <<- list()
  motif.at.tree.level <- leaves.per.node(tree)

  checked.levels <<- rep(0, times = length(tree$merge)/2)

  ## Iterate through the merge of the hierarchical tree
  apply(tree$merge, 1, function(x){

    ## Store the levels
    child1 <- x[1]
    child2 <- x[2]
    ## Store the level of the tree
    level <- which(tree$merge == child1)

#     print(paste("Aligning Level: ", level, " - Checked: ", checked.levels[level]))

    ## If it is indicated, save the levels attributes.
    if(nodes.attributes == TRUE){
      internal.nodes.attributes[[paste("node_", level, sep = "")]][["alignment_flag"]] <<- 0
      internal.nodes.attributes[[paste("node_", level, sep = "")]][["node"]] <<- level
      internal.nodes.attributes[[paste("node_", level, sep = "")]][["method"]] <<- method
      internal.nodes.attributes[[paste("node_", level, sep = "")]][["node_1"]] <<- child1
      internal.nodes.attributes[[paste("node_", level, sep = "")]][["node_2"]] <<- child2
    }


    if(checked.levels[level] == 0){

      ##############################
      ## Case 1: align two leaves
      if ((child1 < 0) && (child2 < 0)) {
        internal.nodes.attributes[[paste("node_", level, sep = "")]][["merge_class"]] <<- 1
        align.two.leaves(child1,
                         child2,
                         thresholds,
                         method,
                         metric = metric,
                         align = align,
                         nodes.attributes = nodes.attributes,
                         motif.at.tree.level = motif.at.tree.level)

      ###########################################################
      ## Case 2: align a leaf with a cluster (already aligned)
      } else if ((child1 < 0) && (child2 > 0)) {
        internal.nodes.attributes[[paste("node_", level, sep = "")]][["merge_class"]] <<- 2
        align.leaf.and.cluster(child1,
                               child2,
                               thresholds,
                               method,
                               metric = metric,
                               align = align,
                               nodes.attributes = nodes.attributes,
                               motif.at.tree.level = motif.at.tree.level)

      ###########################################################
      ## Case 3: align two (already aligned) clusters
      } else if ((child1 > 0) && (child2 > 0)) {
        internal.nodes.attributes[[paste("node_", level, sep = "")]][["merge_class"]] <<- 3
        align.two.clusters(child1,
                           child2,
                           thresholds,
                           method,
                           metric = metric,
                           align = align,
                           nodes.attributes = nodes.attributes,
                           motif.at.tree.level = motif.at.tree.level)
      }

     if(internal.nodes.attributes[[paste("node_", level, sep = "")]][["alignment_flag"]] == 0){
       next.level <- find.next.levels.in.tree(level)
       checked.levels[next.level] <<- 1
      }


      ############################################
      if(intermediate.alignments == TRUE){

        ## Export the intermediate-alignment information in order to create the branch-motifs
        motifs.info.levels[[paste("node_", level, sep = "")]] <<- motifs.info[get.id(motif.at.tree.level[[level]])]
        motifs.info.levels[[paste("node_", level, sep = "")]] <<- sapply(motifs.info.levels[[paste("node_", level, sep = "")]], function(x){
          return(x[c("name", "number", "strand", "consensus_d", "consensus_rc", "spacer.up", "spacer.dw")])
        })
      }

    } else {
      next.level <- find.next.levels.in.tree(level)
      checked.levels[next.level] <<- 1

      if(align == TRUE){
        if((child1 < 0) && (child2 > 0)){
          n1.id <- get.id(-child1)
          motifs.info[[n1.id]][["name"]] <<- get.name(n1.id)
          motifs.info[[n1.id]][["consensus_d"]] <<- get.consensus(n1.id, RC = FALSE)
          motifs.info[[n1.id]][["consensus_rc"]] <<- get.consensus(n1.id, RC = TRUE)
          motifs.info[[n1.id]][["strand"]] <<- "D"
          motifs.info[[n1.id]][["number"]] <<- -child1
          motifs.info[[n1.id]][["spacer.up"]] <<- 0
          motifs.info[[n1.id]][["spacer.dw"]] <<- 0
        }
      }
    }
  })
  return(list(motifs.alignment = motifs.info, node.attributes = internal.nodes.attributes, intermediate.alignments = motifs.info.levels))
}