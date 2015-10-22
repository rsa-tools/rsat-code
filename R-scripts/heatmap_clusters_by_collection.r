# library("gplots")
## Color palette
# grad <- colorRampPalette(c("white", "green", "red", "blue", "black"),space="rgb")
# grad <- colorRampPalette(rev(brewer.pal(11, "Spectral" ) ), space="Lab")
# grad <- colorRampPalette(c("white", "green", "blue"),space="rgb")
# ## Plot the heatmap
# heatmap.2(clusters.matrix,
#
#           ## Select clustering method
#           hclustfun=function(c){hclust(c, method='single')},
#
#           ## Hide/show the trace
#           trace = "none",
#
# #           Colv = FALSE,
#          Rowv = TRUE,
#
#           ## Separate the columns
#           colsep = 1:(length(collection.names)-1),
#           sepcolor = 'white',
#           sepwidth = 0.05,
#
#           ## Set the margins
# #           margins = c(5, 23),
#
#           ## Set the colors of columns, rows and cells
#           #           ColSideColors = color.order,
#           #           RowSideColors = color.order,
#           col = grad,
#           breaks = 1000,
#
#           ## Set the tree positions
#            dendrogram = "row",
#
#           ## Set the col and row labels
#           labCol = collection.names,
#           labRow = clusters.names,
#
#           ## Set the font size
#           cexCol = 0.4,
#           cexRow = 0.07,
#
#           ## Set the key with the values
#           key = FALSE,
# #           keysize = 1.5,
# #           key.xlab = "Nb of motifs",
# #           key.ylab = "",
# #           key.title = "",
# #           density.info = "none",
# )
#
#
# # ## Create matrix for D3 heatmap
# # x <- NULL
# #
# #   for(j in 1:dim(clusters)[2]){
# #     for(i in 1:dim(clusters)[1]){
# #      x <<- rbind(x, matrix(c(j,i, as.numeric(clusters[j,i])), nrow = 1))
# #   }
# # }
# # colnames(x) <- c("Row", "Col", "Value")
# # write.table(x, file = "/home/jcastro/rsat/R-scripts/matrix_heatmap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# #
# # cluster_id <- paste( paste("'", "cluster_", 1:dim(clusters)[1], "' ", sep = ""), collapse = ",")
# # col.row <- 1:dim(clusters)[1]


## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path
source(file.path(dir.rsat, 'R-scripts/config.R'))
library("RColorBrewer")

###########################################
## Read arguments from the command line.
##
## Arguments passed on the command line
## will over-write the default arguments
## specified above.
args <- commandArgs(trailingOnly=TRUE);
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  verbose(args, 3)
}


## Read cluster count table
clusters <- read.table(file = cluster.counts.file, sep = "\t", header = TRUE)
names(clusters) <- gsub("X.Cluster_ID", "Cluster_ID", names(clusters))

#################################
## Create the percentage table

## Count the number of motif per collection
nb.db <- dim(clusters)[2] - 2
motif.DB.counts <- apply(clusters[,3:(nb.db+2)], 2, sum)

percent.table <- NULL
coverage.contingency.table <- NULL
x <- sapply(names(motif.DB.counts), function(DB){

  ## Select those cluster with at least one motif corresponding
  ## to the current motifDB
  DB.motifs <- clusters[clusters[,DB] > 0,]

  #################################################################################
  ## Calculate the overlap between the databases

  ## Select those cluster with at least one motif corresponding
  ## to the current motifDB
  coverage <- apply(DB.motifs[3:dim(DB.motifs)[2]], 2, sum) / motif.DB.counts

  coverage.contingency.table <<- cbind(coverage.contingency.table, matrix(coverage, ncol = 1))

  #################################################################################
  ## Count the number of exclusive motifs of each database

  ## Count the number of motifs that correspond exclusively to a collection of motifs
  DB.motifs.exclusive <- apply(DB.motifs[,3:(nb.db+2)],1, sum)
  DB.motifs.exclusive <- length(DB.motifs.exclusive[DB.motifs.exclusive == 1])

  ## Calculate the percentage of the collection which is unique
  DB.percent <- round(DB.motifs.exclusive / motif.DB.counts[DB], digits = 4)

  ## Calculate the percentage of the total collection corresponding to the
  ## unique motifs of the analyzed motifDB
  Total.percent <- round(DB.motifs.exclusive / sum(motif.DB.counts), digits = 4)

  #   print(paste("Nb Unique motifs: ", DB.motifs.exclusive, " -  %(internal) :", DB.percent, " -  %(total): ", Total.percent))
  percent.table <<- cbind(percent.table, matrix(c(motif.DB.counts[DB], DB.motifs.exclusive, DB.percent, Total.percent), ncol = 1))
})

#########################################################
## Add a new column and re-order the percentage matrix
percent.table <- cbind(percent.table, c("DB_nb_motifs", "Nb_exclusive_motifs", "DB_percent", "Total_percent"))
percent.table <- percent.table[,c(dim(percent.table)[2],1:(dim(percent.table)[2]-1))]
colnames(percent.table) <- c("#Collection", names(clusters[,3:(nb.db+2)]))
percent.table <- t(percent.table)
write.table(percent.table, file = percent.table.file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

#########################################################
## Round and export the coverage contingency table
coverage.contingency.table <- round(coverage.contingency.table, digits = 3)
colnames(coverage.contingency.table) <- names(motif.DB.counts)
rownames(coverage.contingency.table) <- names(motif.DB.counts)
write.table(coverage.contingency.table, file = coverage.table.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Convert the coverage table to the format required in D3 heatmap
y <- NULL
for(j in 1:dim(coverage.contingency.table)[1]){
  for(i in 1:dim(coverage.contingency.table)[2]){
    y <<- rbind(y, matrix(c(j,i, as.numeric(coverage.contingency.table[j,i])), nrow = 1))
  }
}
colnames(y) <- c("Row", "Col", "Value")
verbose(paste("Exporting data with collection coverage for D3", coverage.table.d3), 2)
write.table(y, file = coverage.table.d3, sep = "\t", quote = FALSE, row.names = FALSE)

###########################################################
## Create attributes table to fill the D3 coverage fields
row.nb <- dim(coverage.contingency.table)[1]
col.nb <- dim(coverage.contingency.table)[2]
default.labels <- paste(paste("'", names(motif.DB.counts), "'", sep = ""), collapse = ",")
default.number <- paste(1:length(motif.DB.counts), collapse = ",")
left <- (max(as.vector(sapply(names(motif.DB.counts), nchar))) + 2) * 10
cell.size <- 20
bottom <- 120
legend.header <- bottom - 35
if(row.nb < 5){
  bottom <- 120
  legend.header <- bottom - 35
} else if(row.nb < 8){
  bottom <- 170
  legend.header <- bottom - 35
} else if(row.nb < 13){
  bottom <- 220
  cell.size <- 15
  legend.header <- bottom - 27
} else if(row.nb < 18){
  bottom <- 270
  cell.size <- 15
  legend.header <- bottom - 27
}
coverage.info <- matrix(c("Collection_labels", default.labels,
                       "Collection_number", default.number,
                       "Left_space", left,
                       "Bottom_space", bottom,
                       "Col_number", col.nb,
                       "Row_number", row.nb,
                       "Legend_Head", legend.header
), nrow = 2)
coverage.info.df <- t(data.frame(coverage.info))
write.table(coverage.info.df, file = coverage.heatmap.attributes.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



##################################################
## Create the collection's contribution heatmap

## Get the values + names
clusters.names <- as.vector(clusters[,1])
collection.names <- names(clusters[3:(dim(clusters)[2])])
clusters <- clusters[3:(dim(clusters)[2])]
clusters.matrix <- as.matrix(clusters)

step <- 5
if(max(clusters) < 10){
  step <- 1
} else if(max(clusters) < 21){
  step <- 2
} else if(max(clusters) < 31){
  step <- 3
} else if(max(clusters) < 41){
  step <- 4
} else if(max(clusters) > 51){
  step <- 5
}

###########################################
## Create Gradient Hexadecimal:
## Given X hexa colors creates a color
## palette.
## This is exported and will be read later
## in the D3 heatmap code.
# rgb.palette <- colorRampPalette(c("#75FD3A", "#1F6800"), space = "rgb")
# rgb.palette <- colorRampPalette(c("#FFE991", "#FEAD23", "#930047"), space = "rgb")
# rgb.palette <- colorRampPalette(c("#FFE991", "#930047"), space = "rgb")
rgb.palette <- colorRampPalette(brewer.pal(9, "YlOrRd"), space="Lab")
white <- "#FFFFFF"
white <- append(white,rgb.palette(ceiling((max(clusters)/step+1))))

###############################################################
## Run the hierarchical clustering with the three methods
## (average + complete + single). Save the order of the nodes
order.list <<- list()
order.list.names <<- list()
for(m in c("average", "complete", "single")){
  tree <- hclust(dist(clusters.matrix), method = m)
  order.list[[m]] <- paste(tree$order, collapse = ",")
  order.list.names[[m]] <- paste(paste("'cluster_", tree$order, "'", sep = ""), collapse = ",")
}

###############################################
## Parse the Heatmap table format used in D3
## This table is printed in a new file
x <- data.frame(t(clusters))
names(x) <- paste("cluster_", 1:dim(clusters)[1], sep = "")
y <- NULL
for(j in 1:dim(x)[1]){
  for(i in 1:dim(x)[2]){
    y <<- rbind(y, matrix(c(j,i, as.numeric(x[j,i])), nrow = 1))
  }
}
colnames(y) <- c("Row", "Col", "Value")
verbose(paste("Exporting heatmap with cluster by collection table for D3", heatmap.table.d3), 2)
write.table(y, file = heatmap.table.d3, sep = "\t", quote = FALSE, row.names = FALSE)


############################
## Output data (to print)

## Color palette in Hexa code
gradient <- paste("[", paste(paste("'", white, "'", sep=""), collapse=","), "];", sep = "")

## Get the clusters names orderer according the linkage method
cluster.names <- paste(paste("'cluster_", 1:dim(clusters)[1], "'", sep =""), collapse=",")
# average.names <- order.list.names[["average"]]
# complete.names <- order.list.names[["complete"]]
# single.names <- order.list.names[["single"]]

## Get the clusters number orderer according the linkage method
cluster.number <- paste(1:dim(clusters)[1], collapse=",")
average.number <- order.list[["average"]]
complete.number <- order.list[["complete"]]
single.number <- order.list[["single"]]

## Default names
default.names <- paste(paste("'cluster_", 1:dim(clusters)[1], "'", sep = ""), collapse = ",")
default.number <- paste(1:dim(clusters)[1], collapse = ",")

## Heatmap variables
col.nb <- dim(clusters)[1]
row.nb <- dim(clusters)[2]

## Row
heatmap.rows.nb <- paste(1:row.nb, collapse=",")
heatmap.rows.name <- paste(paste("'", collection.names, "'", sep = ""), collapse = ",")

## Collections
collections <- paste(paste("'", collection.names, "'", sep = ""), collapse = ",")

## Range to color the values
domain.nb <- seq(from = 0, to = max(clusters), by = step)
domain.nb[1] <- 1
domain <- paste(domain.nb, collapse=",")

## Legend
legend <- c(0,domain.nb)
legend <- paste(legend, collapse=",")

## Right space
left <- (max(as.vector(sapply(collection.names, nchar))) + 2) * 10

## Div bottom + Cell size
cell.size <- 20
bottom <- 120
legend.header <- bottom - 35
if(row.nb < 5){
  bottom <- 120
  legend.header <- bottom - 35
} else if(row.nb < 8){
  bottom <- 170
  legend.header <- bottom - 35
} else if(row.nb < 13){
  bottom <- 220
  cell.size <- 15
  legend.header <- bottom - 27
} else if(row.nb < 18){
  bottom <- 270
  cell.size <- 15
  legend.header <- bottom - 27
}

order.info <- matrix(c("Gradient", gradient,
                       "Cluster_names", cluster.names,
                       "Cluster_number", cluster.number,
                       "Average_number", average.number,
                       "Complete_number", complete.number,
                       "Single_number", single.number,
                       "Cell_size", cell.size,
                       "Col_number", col.nb,
                       "Row_number", row.nb,
                       "Row_order_default", heatmap.rows.nb,
                       "Domain", domain,
                       "Legend", legend,
                       "Legend_Head", legend.header,
                       "Left_space", left,
                       "Bottom_space", bottom,
                       "Collections", collections
                       ), nrow = 2)
order.info.df <- t(data.frame(order.info))
verbose(paste("Exporting table with the order of the clusters (Required in D3 Heatmap)", order.list.file), 2)
write.table(order.info.df, file = attributes.list.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
