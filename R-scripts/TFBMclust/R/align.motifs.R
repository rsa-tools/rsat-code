align.motifs <- function(hclust.tree,
                         desc.table,
                         compa.table,
                         thresholds = list(Ncor = 0.4, cor = 0.6, w = 5),
                         method = "average",
                         metric = "Ncor",
                         nodes.attributes = TRUE,
                         intermediate.alignments = FALSE){

  motifs.info <<- list()
  motifs.info.levels <- list()
  internal.nodes.attributes <<- list()
  motif.at.tree.level <- leaves.per.node(tree)

  ## Iterate through the merge of the hierarchical tree
  apply(hclust.tree$merge, 1, function(x){

    ## Store the nodes
    child1 <- x[1]
    child2 <- x[2]
    ## Store the level of the tree
    level <- which(tree$merge == child1)

    ## If it is indicated, save the nodes attributes.
    if(nodes.attributes == TRUE){
      internal.nodes.attributes[[paste("level_", level, sep = "")]][["level"]] <<- level
      internal.nodes.attributes[[paste("level_", level, sep = "")]][["method"]] <<- method
      internal.nodes.attributes[[paste("level_", level, sep = "")]][["node_1"]] <<- child1
      internal.nodes.attributes[[paste("level_", level, sep = "")]][["node_2"]] <<- child2
    }

    ##############################
    ## Case 1: align two leaves
    if ((child1 < 0) && (child2 < 0)) {
      align.two.leaves(child1,
                       child2,
                       desc.table,
                       compa.table,
                       thresholds,
                       method,
                       metric,
                       hclust.tree,
                       nodes.attributes = nodes.attributes,
                       motif.at.tree.level = motif.at.tree.level)

    ###########################################################
    ## Case 2: align a leaf with a cluster (already aligned)
    } else if ((child1 < 0) && (child2 > 0)) {
      align.leaf.and.cluster(child1,
                             child2,
                             desc.table,
                             compa.table,
                             thresholds,
                             method,
                             metric,
                             hclust.tree,
                             nodes.attributes = nodes.attributes,
                             motif.at.tree.level = motif.at.tree.level)

    ###########################################################
    ## Case 3: align two (already aligned) clusters
    } else if ((child1 > 0) && (child2 > 0)) {
      align.two.clusters(child1,
                         child2,
                         desc.table,
                         compa.table,
                         thresholds,
                         method,
                         metric,
                         hclust.tree,
                         nodes.attributes = nodes.attributes,
                         motif.at.tree.level = motif.at.tree.level)
    }


    ############################################
    if(intermediate.alignments == TRUE){

      ## Export the intermediate-alignment information in order to create the branch-motifs
      motifs.info.levels[[paste("merge_level_", level, sep = "")]] <<- motifs.info[get.id(motif.at.tree.level[[level]], desc.table)]
      motifs.info.levels[[paste("merge_level_", level, sep = "")]] <<- sapply(motifs.info.levels[[paste("merge_level_", level, sep = "")]], function(x){
        return(x[c("name", "number", "strand", "consensus_d", "consensus_rc", "spacer.up", "spacer.dw")])
      })
    }
  })

  return(list(motifs.alignment = motifs.info, node.attributes = internal.nodes.attributes, intermediate.alignments = motifs.info.levels))
}