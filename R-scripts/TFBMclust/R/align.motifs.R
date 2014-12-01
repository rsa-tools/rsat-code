align.motifs <- function(tree, desc.table, compa.table, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), score = "Ncor", method = "average", metric = "Ncor"){

  motifs.info <<- list()

  apply(tree$merge,1,function(x){
    child1 <- x[1]
    child2 <- x[2]
    level <- which(tree$merge == child1)

    ##############################
    ## Case 1: align two leaves
    if ((child1 < 0) && (child2 < 0)) {
      align.two.leaves(child1, child2, desc.table, compa.table, score, thresholds, method, metric, hclust.tree)
    }
  })
  return(motifs.info)
}


  for (merge.level in 1:nrow(tree$merge)) {
    ##for (merge.level in 1:5) {
    merge.level <<- merge.level

    child1 <- tree$merge[merge.level,1]
    child2 <- tree$merge[merge.level,2]

#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merge_level"]] <<- merge.level
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["method"]] <<- hclust.method
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["node_1"]] <<- child1
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["node_2"]] <<- child2


    ########################################
    ## Case 1: merging between two leaves ##
    ########################################
    if ((child1 < 0) && (child2 < 0)) {

      align.two.leaves(child1, child2)

      ## Identify the nodes
      n1 <- min(-child1,-child2)
      n2 <- max(-child1,-child2)

      internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merging_type"]] <<- 1
    }

#     ############################################
#     ## Case 2: merging a motif with a cluster ##
#     ############################################
#     if(((child1 < 0) && (child2 > 0)) || ((child1 > 0) && (child2 < 0))){
#
#       align.leave.and.cluster(child1, child2)
#
#       ## Identified the nodes
#       n1 <- abs(min(child1, child2))
#       n.aligned <- merge.levels.leaves[[merge.level]][which(merge.levels.leaves[[merge.level]] != n1)]
#       n2 <- n.aligned
#
#       internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merging_type"]] <<- 2
#     }
#
#
#     ##########################################
#     ## Case 3: merging between two clusters ##
#     ##########################################
#     if ((child1 > 0) && (child2 > 0)) {
#       align.clusters(child1, child2)
#
#       ## Identified the nodes
#       N1 <- abs(min(child1, child2))
#       N2 <- abs(max(child1, child2))
#
#       n1 <- merge.levels.leaves[[N1]]
#       n2 <- merge.levels.leaves[[N2]]
#
#       internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["merging_type"]] <<- 3
#     }
#
#
#     if(hclust.method == "single"){
#       ## Get the id of each node
#       central.motifs <- central.motifs.ids(get.id(n1), get.id(n2))
#       id1 <- central.motifs[1]
#       id2 <- central.motifs[2]
#     } else if(hclust.method == "complete"){
#       farthest.motifs <- farthest.motifs.ids(get.id(n1), get.id(n2))
#       id1 <- farthest.motifs[1]
#       id2 <- farthest.motifs[2]
#     }
#
#     ## Nodes attributes
#     attributes <- attributes.among.clusters(get.id(n1), get.id(n2))
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["cluster_1"]] <<- paste(n1, collapse = " ")
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["cluster_2"]] <<- paste(n2, collapse = " ")
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["min_score"]] <<- attributes[1]
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["max_score"]] <<- attributes[2]
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["median_score"]] <<- attributes[3]
#
#     if(hclust.method == "average"){
#       aligned.motif.flag <- alignment.status(get.id(n1), get.id(n2), hclust.method)
#     }else{
#       aligned.motif.flag <- alignment.status(id1, id2, hclust.method)
#     }
#
#     internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["flag"]] <<- aligned.motif.flag
#     if(aligned.motif.flag == 0){
#       internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["alignment_status"]] <<- "Non-aligned"
#     }else{
#       internal.nodes.attributes[[paste("merge_level_", merge.level, sep = "")]][["alignment_status"]] <<- "Aligned"
#     }
#   }
}