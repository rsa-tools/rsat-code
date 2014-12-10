###########################################################
## The clusters are selected by traversing the tree with
## a bottom-up approach, rather than top-down as is used
## with the R-defined function cutree.
define.clusters.bottom.up <- function(attributes.list, tree, desc.table){

  ########################################################################
  ## Search the name of the cluster which correspond to the given nodes
  ## NOTE: This function is called only inside the define.clusters.bottom.up() function
  return.cluster.name <- function(nodes, clusters.list){

    matches <- sapply(clusters.list, function(X){
      sum(nodes %in% X)
    })

    if(sum(matches[matches == length(nodes)]) == 0){
      return(0)
    } else{
      matches <- matches[matches == length(nodes)]
      return(names(matches))
    }
  }

  cluster.number <- 0
  level <<- 0
  clusters <<- list()
  attributes.cp <- attributes.list

  sapply(names(attributes.cp), function(X){

    level <<- level + 1

    class.level <- paste("class_", as.numeric(attributes.cp[[X]][["merge_class"]]), sep = "")
    switch(class.level,

      ## merge_class = 1
      ## Alignment of two leaves
      class_1 = {

        cluster.number <<- cluster.number + 1

        ## Creates one cluster with the motifs if the flag = 1
        if(as.numeric(attributes.cp[[X]][["alignment_flag"]]) == 1){

          ## Get the numbers of the nodes
          nodes.numbers <- leaves.per.node(tree)[[level]]

          ## Store the numbers in the list
          clusters[[paste("cluster_", cluster.number, sep = "")]] <<- nodes.numbers

        ## Conversely, two clusters are created
        } else{
          clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(attributes.cp[[X]][["cluster_1"]])
          cluster.number <<- cluster.number + 1


          clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(attributes.cp[[X]][["cluster_2"]])
          attributes.cp[[X]][["alignment_flag"]] <<- 0
        }
      },

      ## merge_class = 2
      ## Alignment of one leaf and one cluster
      class_2 = {

        if(as.numeric(attributes.cp[[X]][["alignment_flag"]]) == 1){

          ## Obtain the status of the cluster to merge
          prev.cluster.nb <- as.numeric(attributes.cp[[X]][["node_2"]])
          prev.status <- as.numeric(attributes.cp[[paste("level_", prev.cluster.nb, sep = "")]][["alignment_flag"]])

          ## Following the same order than the tree$merge
          ## If the cluster was succesfully aligned:
          ## add the number of the leaf to the number of the aligned cluster
          if(prev.status == 1){
            cluster.name <- return.cluster.name(leaves.per.node(tree)[[prev.cluster.nb]], clusters)
            clusters[[cluster.name]] <<- append(clusters[[cluster.name]], as.numeric(attributes.cp[[X]][["cluster_1"]]))

          ## Conversely, creates a new cluster with the single leaf
          }else{
            cluster.number <<- cluster.number + 1
            clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(attributes.cp[[X]][["cluster_1"]])
            attributes.cp[[X]][["flag"]] <<- 0
          }

        ## Conversely, creates a new cluster with the single leaf
        } else{

          cluster.number <<- cluster.number + 1
          clusters[[paste("cluster_", cluster.number, sep = "")]] <<- as.numeric(attributes.cp[[X]][["cluster_1"]])
          attributes.cp[[X]][["flag"]] <<- 0
        }
      },

      ## merge_class = 3
      ## Alignment of two clusters
      class_3 = {

        ## Obtain the status of the two clusters before merge them
        prev.cluster.1.nb <- as.numeric(attributes.cp[[X]][["node_1"]])
        prev.status.1 <- as.numeric(attributes.cp[[prev.cluster.1.nb ]][["alignment_flag"]])

        prev.cluster.2.nb <- as.numeric(attributes.cp[[X]][["node_2"]])
        prev.status.2 <- as.numeric(attributes.cp[[prev.cluster.2.nb ]][["alignment_flag"]])

        if(as.numeric(attributes.cp[[X]][["alignment_flag"]]) == 1){
          if(prev.status.1 == 0 | prev.status.2 == 0){
            attributes.cp[[X]][["flag"]] <<- 0
          } else{
            cluster.name <- return.cluster.name(leaves.per.node(tree)[[prev.cluster.1.nb]], clusters)
            cluster.2.nodes <- as.numeric(unlist(strsplit(attributes.cp[[X]][["cluster_2"]], " ")))
            cluster.name.erase <- return.cluster.name(cluster.2.nodes, clusters)
            clusters[[cluster.name]] <<- append(clusters[[cluster.name]], cluster.2.nodes)
            clusters[[cluster.name.erase]] <<- NULL
            #rm(clusters[[cluster.name.erase]])
          }
        }else{
          attributes.cp[[X]][["alignment_flag"]] <<- 0
        }
      }
    )
  })


  ## Re-name the clusters
  names(clusters) <<- paste("cluster_", seq(1:length(clusters)), sep = "")
  clusters <- lapply(clusters,function(x){get.id(x, desc.table)})
  return(clusters)
}