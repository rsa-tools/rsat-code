####################################################
## Draw a heatmap usig a particular allgomeration
## method.
draw.heatmap.motifs <- function(dist.table, method = "average", clusters.list, alignment.list){

  dist.table <- as.matrix(dist.table)

  ## Set the suported colors
  color <- rainbow(length(clusters.list))
  cluster.to.color <- list(length(clusters.list))
  color.counter <- 0


  ## Fill a list where each element correspond to a cluster name
  ## and the value its corresponding color
  sapply(names(clusters.list), function(x){
    color.counter <<- color.counter + 1
    cluster.to.color[[x]] <<- color[color.counter]
  })

  ## Create a vector with the the corresponding color of each motif
  color.order <- vector()
  sapply(colnames(dist.table), function(y){

    current.cluster <- names(unlist(sapply(clusters.list, function(x){ grep(y,x, value = TRUE)})))
    color.order <<- append(color.order, cluster.to.color[[current.cluster]])
  })

  ## Color gradient for the heatmap
  grad <- colorRampPalette(c("green", "black", "red"))(n = 299)


  ## Calculate the bottom border
  bottom <- 0.277 * max(nchar(colnames(dist.table)))
  par(oma=c(bottom,0,0,2), family="mono")

  # Get the aligned consensuses, which will be used as the Row names
  consensus <-sapply(colnames(dist.table), function(x){
    as.vector(alignment.list[[x]][["consensus"]])
  })
  consensus <- as.vector(consensus)

  # Get the aligned consensuses, which will be used as the Row names
  columns.heatmap <- sapply(colnames(dist.table), function(x){
    if(x == alignment.list[[x]][["name"]]){
      as.vector(x)
    } else{
      paste(x, alignment.list[[x]][["name"]], sep = "  ")
    }
  })
  columns.heatmap <- as.vector(columns.heatmap)

  ## Draw the heatmap
  heatmap(dist.table,
          hclustfun = function(d){hclust(d, method = method)},
          symm=TRUE,
          ColSideColors = color.order,
          RowSideColors = color.order,
          Rowv = as.dendrogram(tree),
          Colv =as.dendrogram(tree),
          col = grad,
          labRow = consensus,
          labCol = columns.heatmap,
          cexRow = 30/length(consensus) + 0.05,
          cexCol = 30/length(consensus) + 0.05
#           cexRow = 0.002132353*length(consensus),
#           cexCol = 0.002132353*length(consensus)
  )
}



