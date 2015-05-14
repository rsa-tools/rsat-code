align.motifs <- function(hclust.tree, desc.table, compa.table, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", nodes.attributes = TRUE, intermediate.alignments = FALSE){

  motifs.info <<- list()
  motifs.info.levels <- list()
  internal.nodes.attributes <<- list()

  apply(hclust.tree$merge, 1, function(x){
    child1 <- x[1]
    child2 <- x[2]
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
      align.two.leaves(child1, child2, desc.table, compa.table, thresholds, method, metric, hclust.tree, nodes.attributes = nodes.attributes)
    }

    ###########################################################
    ## Case 2: align a leaf with a cluster (already aligned)
    if ((child1 < 0) && (child2 > 0)) {
      align.leaf.and.cluster(child1, child2, desc.table, compa.table, thresholds, method, metric, hclust.tree, nodes.attributes = nodes.attributes)
    }

    ###########################################################
    ## Case 3: align two (already aligned) clusters
    if ((child1 > 0) && (child2 > 0)) {
      align.two.clusters(child1, child2, desc.table, compa.table, thresholds, method, metric, hclust.tree, nodes.attributes = nodes.attributes)
    }


    ############################################
    if(intermediate.alignments == TRUE){

      ## Export the intermediate-alignment information in order to create the branch-motifs
      motifs.info.levels[[paste("merge_level_", level, sep = "")]] <<- motifs.info[get.id(leaves.per.node(tree)[[level]], desc.table)]
      motifs.info.levels[[paste("merge_level_", level, sep = "")]] <<- sapply(motifs.info.levels[[paste("merge_level_", level, sep = "")]], function(x){
        return(x[c("name", "number", "strand", "consensus_d", "consensus_rc", "spacer.up", "spacer.dw")])
      })
    }
  })

  return(list(motifs.alignment = motifs.info, node.attributes = internal.nodes.attributes, intermediate.alignments = motifs.info.levels))
}