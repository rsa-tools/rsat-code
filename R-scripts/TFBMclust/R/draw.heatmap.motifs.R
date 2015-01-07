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
  #grad <- colorRampPalette(c("green", "black", "red"))(n = 299)
  grad <- greenred(200)

  ## Calculate the bottom border
  rigth <- round(170/length(alignment.list) + (length(alignment.list)/2 * 0.001), digits = 2)
  bottom <- round(190/length(alignment.list) + (length(alignment.list)/2 * 0.001), digits = 2)

  par(oma=c(bottom,0.5,0.5,rigth), family="mono")

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
  heatmap.2(dist.table,

          ## Display the tree symetrically
          trace = "none",
          symm=TRUE,

          ## Set the colors of columns, rows and cells
          ColSideColors = color.order,
          RowSideColors = color.order,
          col = grad,

          ## The order of the values is set according these dendrograms
          Rowv = as.dendrogram(tree),
          Colv =as.dendrogram(tree),

          ## To show only the dendrogram in row
          dendrogram = "column",

          ## Set the col and row labels
          labRow = consensus,
          labCol = columns.heatmap,

          ## Set the margins
          cexRow = 20/length(consensus) + 0.1,
          cexCol = 20/length(consensus) + 0.02,

          ## Set the key with the values
          key = TRUE,
          keysize = 1.5,
          key.xlab = "Distance",
          key.ylab = "",
          key.title = "Color key",
          density.info = "none"
  )
}



