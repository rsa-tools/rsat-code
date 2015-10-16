## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path
source(file.path(dir.rsat, 'R-scripts/config.R'))
library("RColorBrewer")
library("gplots")

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


# matrix.scan.results <- "/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/heatmap_matrix_matches/test/matrix_scan_output_test.tab"
# matrix.scan.results <- "/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/heatmap_matrix_matches/test/CapStarrseq_Active_Prom_K562_merge_IP_matrix_scan_pval_1e-4_HOCOMOCO_bg_mkv_2.tab"

##############################
## Read matrix-scan table
# matrix.scan.results <- "" ## Read from command-line arguments
scan.results <- read.csv(file = matrix.scan.results, sep = "\t", header = TRUE, comment.char = ";")
names(scan.results) <- gsub("X.seq_id", "seq_id", names(scan.results))

#################
## Set p-value
p.val <- as.numeric(p.val)

############################
## Get the sequences ID's
seq.id <- unique(as.vector(scan.results[scan.results$ft_type == "limit",1]))

#########################
## Remove the 'limits'
scan.results <- scan.results[scan.results$ft_type != "limit",]

#############################
## Get tha matrices name's
matrix.names <- unique(as.vector(scan.results$ft_name))

##############################
## Create the matches table
count.matches.tab <- NULL
count.matches.tab <- sapply(seq.id, function(seq){
     table(scan.results[scan.results$seq_id == seq & scan.results$Pval <= p.val,]$ft_name)
})
count.matches.tab <- t(count.matches.tab)
count.matches.tab <- count.matches.tab[,1:(dim(count.matches.tab)[2] - 1)]


#################

# ## Set the suported colors
# nb.clusters <- 25
# clusters.names <- paste("cluster_", 1:nb.clusters, sep = "")
# color <- rainbow(nb.clusters)
# cluster.to.color <- list()
# color.counter <- 0
# clusters <- cutree(tree, k = nb.clusters)
# 
# 
# ## Fill a list where each element correspond to a cluster name
# ## and the value its corresponding color
# sapply(1:nb.clusters, function(x){
#   color.counter <<- color.counter + 1
#   cluster.to.color[[x]] <<- color[color.counter]
# })
# 
# ## Create a vector with the the corresponding color of each motif
# current.cluster <- sapply(1:length(names(clusters)), function(y){
#    as.vector(clusters[names(clusters) == as.character(y)])
# })
# 
# color.order <- sapply(current.cluster, function(c){
#   cluster.to.color[[c]]
# })



########################

pdf(heatmap.pdf.file)
palette <- colorRampPalette(c("#FFE991", "#FF8000", "#930047"), space = "rgb")
white <- "#FFFFFF"
palette <- append(white, palette(300))
# palette <- colorRampPalette(brewer.pal(11, "Spectral"), space="Lab")

a <- (count.matches.tab)
# rownames(a) <- 1:length(rownames(a))
# colnames(a) <- 1:length(colnames(a))
plot(hclust(dist(a), method = "complete"))
tree <- hclust(dist(a), method = "complete")
#cutree(tree, k = 7)
heatmap.2(a,
          
          ## Display the tree symetrically
          trace = "none",
          
          ## Set the colors of columns, rows and cells
#           ColSideColors = color.order,
#           RowSideColors = color.order,
          col = palette,
          
          ## Set the font size
          cexRow = 0.065,
          cexCol = 0.065,
          
          ## Set the key with the values
          key = TRUE,
          keysize = 1.5,
          key.xlab = "Number of matches",
          key.ylab = "",
          density.info = "none"
          )
dev.off()

###################################################

###############################################################
## Run the hierarchical clustering with the three methods
## (average + complete + single). Save the order of the nodes
order.list <<- list()
order.list.names <<- list()

for(ft in c("col", "row")){
  
  if(ft == "col"){
    clusters.matches <- dist(as.matrix(t(count.matches.tab)))
  } else if(ft == "row"){
    clusters.matches <- dist(as.matrix(count.matches.tab))
  }
  
  for(m in c("average", "complete", "single")){
    tree <- hclust(dist(clusters.matches), method = m)
    order.list[[m]][[ft]] <- paste(tree$order, collapse = ",")
    order.list.names[[m]][[ft]] <- paste(paste("'cluster_", tree$order, "'", sep = ""), collapse = ",")
  }
}
average.number.col <- order.list[["average"]][["col"]]
complete.number.col <- order.list[["complete"]][["col"]]
single.number.col <- order.list[["single"]][["col"]]
average.number.row <- order.list[["average"]][["row"]]
complete.number.row <- order.list[["complete"]][["row"]]
single.number.row <- order.list[["single"]][["row"]]

###############################################
## Convert the table to the format required
## for the D3 heatmap
matches.tsv <<- NULL

x <- sapply(1:dim(count.matches.tab)[1], function(j){
  sapply(1:dim(count.matches.tab)[2], function(i){
    matches.tsv <<- rbind(matches.tsv, matrix(c(j,i, as.numeric(count.matches.tab[j,i])), nrow = 1))
  })
})
colnames(matches.tsv) <- c("Row", "Col", "Value")


# for(j in 1:dim(x)[1]){
#   for(i in 1:dim(x)[2]){
#     matches.tsv <<- rbind(matches.tsv, matrix(c(j,i, as.numeric(x[j,i])), nrow = 1))
#   }
# }

######################################################
## Export the table that will be read by D3 heatmap
# heatmap.table.d3 <- "" ## Read from command-line arguments
write.table(matches.tsv, file = matches.tsv.file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(count.matches.tab, file = matches.tab.file, sep = "\t", quote = FALSE)


###########################################
## Create Gradient Hexadecimal:
## Given X hexa colors creates a color palette.

step <- max(count.matches.tab)-1
if(max(count.matches.tab) < 10){
  step <- 1
} else if(max(count.matches.tab) < 21){
  step <- 2
} else if(max(count.matches.tab) < 31){
  step <- 3
} else if(max(count.matches.tab) < 41){
  step <- 4
} else if(max(count.matches.tab) > 51){
  step <- 5
}

rgb.palette <- colorRampPalette(c("#FFE991", "#930047"), space = "rgb")
white <- "#FFFFFF"
white <- append(white,rgb.palette(ceiling((max(count.matches.tab)/step+1))))

################################
## Color palette in Hexa code
gradient <- paste("[", paste(paste("'", white, "'", sep=""), collapse=","), "];", sep = "")

###############################
## Column names => Matrix ID
## Row names => Seq ID
col.nb <- dim(count.matches.tab)[2]
row.nb <- dim(count.matches.tab)[1]
column.heatmap <- paste(paste("'", matrix.names,"'", sep =""), collapse=",")
row.heatmap <- paste(paste("'", seq.id,"'", sep =""), collapse=",")

###############################
## Column Default Order
column.default.order <- paste(1:length(matrix.names), collapse=",")
row.default.order <- paste(1:length(seq.id), collapse=",")

############
## Legend
domain.nb <- seq(from = 1, to = max(count.matches.tab), by = step)
domain <- paste(domain.nb, collapse=",")
legend <- c(0,domain.nb)
legend <- paste(legend, collapse=",")

############
## Borders + Cell size
left <- (max(as.vector(sapply(seq.id, nchar))) + 2) * 10
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
                       "Matrix_name", column.heatmap,
                       "Matrix_number", column.default.order,
                       "Row_order_default", row.default.order,
                       "Average_col", average.number.col,
                       "Complete_col", complete.number.col,
                       "Single_col", single.number.col,
                       "Average_row", average.number.row,
                       "Complete_row", complete.number.row,
                       "Single_row", single.number.row,
                       "Cell_size", cell.size,
                       "Col_number", col.nb,
                       "Row_number", row.nb,
                       "Domain", domain,
                       "Legend", legend,
                       "Legend_Head", legend.header,
                       "Left_space", left,
                       "Bottom_space", bottom,
                       "Seq_IDs", row.heatmap
), nrow = 2)
order.info.df <- t(data.frame(order.info))

######################
## Export the table
# attributes.list.file <- "/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/heatmap_matrix_matches/matrix_scan_heatmap_attributes.tsv"
write.table(order.info.df, file = attributes.list.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
