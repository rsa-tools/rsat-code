## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path
source(file.path(dir.rsat, 'R-scripts/config.R'))
suppressPackageStartupMessages(library("RColorBrewer", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("gplots", warn.conflicts=FALSE, character.only = TRUE))

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


#  matrix.scan.results <- "/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/heatmap_matrix_matches/test/matrix_scan_output_test.tab"
# matrix.scan.results <- "/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/heatmap_matrix_matches/test/CapStarrseq_Active_Prom_K562_merge_IP_matrix_scan_pval_1e-4_HOCOMOCO_bg_mkv_2.tab"
# matrix.scan.results <- "/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/motifs/scan/Yeast_matches_mkv_1.tab"


##############################
## Read matrix-scan table
# matrix.scan.results <- "" ## Read from command-line arguments
scan.results <- read.csv(file = matrix.scan.results, sep = "\t", header = TRUE, comment.char = ";")
names(scan.results) <- gsub("X.seq_id", "seq_id", names(scan.results))

#################
## Set p-value
p.val <- as.numeric(p.val)

#################################
## Get the sequences + motif IDs 
#seq.id <- unique(as.vector(scan.results$seq_id))
seq.id <- unique(as.vector(scan.results[scan.results$ft_type == "limit",1]))

#########################
## Remove the 'limits'
scan.results <- scan.results[scan.results$ft_type != "limit",]
matrix.names.no.filter <- unique(as.vector(scan.results$ft_name))

#######################
## Filter by p-value 
scan.results <- scan.results[scan.results$Pval <= p.val,]
filtered.seq.id <- unique(as.vector(scan.results$seq_id))

#############################
## Get the matrices name's
matrix.names <- unique(as.vector(scan.results$ft_name))

# ## Get the name of the sequences that have no match
# ## after the p-value filtering 
# lost.motifs <- setdiff(matrix.names.no.filter, matrix.names)
# lost.motifs.df <- NULL
# if(length(lost.motifs) > 0){
#   
#   ## Create a DataFrame with 0's. 
#   ## Each row corresponds to the filtered sequences
#   lost.motifs.df <- data.frame(matrix(rep(rep(0, length(lost.motifs)), length(seq.id)), ncol = length(lost.motifs)))
#   colnames(lost.motifs.df) <- lost.motifs
# }

##############################
## Create the matches table
count.matches.tab <- NULL
count.matches.tab <- sapply(seq.id, function(seq){
     table(scan.results[scan.results$seq_id == seq,]$ft_name)
})
count.matches.tab <- t(count.matches.tab)
count.matches.tab <- count.matches.tab[,1:(dim(count.matches.tab)[2] - 1)]


###########################
## Set heatmap key color

# palette <- colorRampPalette(c("#FFE991", "#F7D358", "#930047"), space = "rgb")
# white <- "#FFFFFF"
# palette <- append(white, palette(300))

rgb.palette <- colorRampPalette(brewer.pal(9, "YlOrRd"), space="Lab")
palette <- "#FFFFFF"
palette <- append(palette, rgb.palette(300))


########################################################
## If it is required change the counts of the matches
## in a bollean expression.
## 0 == no match
## 1 == at least one match
if(count.mode == "presence"){
  
  c.names <- colnames(count.matches.tab)
  r.names <- rownames(count.matches.tab)
  col.nb <- dim(count.matches.tab)[2]

  count.matches.tab <- count.matches.tab > 0
  count.matches.tab <- matrix(as.numeric(count.matches.tab),ncol = col.nb)
  rownames(count.matches.tab) <- r.names
  colnames(count.matches.tab) <- c.names
  
  colors <- colorRampPalette(c("#FFFFFF", "#800026"))
  palette <- colors(2)
}

write.table(count.matches.tab, file = matches.tab.file, sep = "\t", quote = FALSE)

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

if(draw.heatmap == 1){

  heatmap.pdf.file <- paste(prefix, "_R_heatmap_matches.pdf", sep = "")
  pdf(heatmap.pdf.file)
  ########################
#   for (m in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){
#     for(d in c("euclidean", "maximum", "manhattan", "binary")){

  for (m in c("ward.D", "single", "complete", "average", "centroid")){
    for(d in c("euclidean", "manhattan")){

  # for (m in c("ward")){
  #   for(d in c("manhattan", "canberra")){
      
      # a <- (count.matches.tab)
      # rownames(a) <- 1:length(rownames(a))
      # colnames(a) <- 1:length(colnames(a))
      # plot(hclust(dist(a), method = "complete"))
      # tree <- hclust(dist(a), method = "complete")
      #cutree(tree, k = 7)
      
      heatmap.2(count.matches.tab,
                
                # plot labels
                main = paste("Link: ", m , " - Dist: ", d, sep = ""),
                ylab = "Sequence ID",
                xlab = "Motif ID",
                
                ## Set distance calculation method
                distfun = function(x){ dist(x, method = d) },
                
                ## Set hclust method
                hclustfun = function(x){ hclust(x, method = m) },
                
                ## Remove the trace
                trace = "none",
                
                ## Set the colors of columns, rows and cells
                #           ColSideColors = color.order,
                #           RowSideColors = color.order,
                col = palette,
                
                ## Set the font size
                cexRow = 0.06,
                cexCol = 0.06,
                
                ## Set the key with the values
                key = TRUE,
                keysize = 1,
                key.xlab = "Ocurrences",
                key.ylab = "",
                density.info = "none"
      ) 
    }
  }
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
    
    
    for(m in c("average", "centroid", "complete", "median", "single", "ward.D", "ward.D2", "mcquitty")){
      tree <- hclust(dist(clusters.matches), method = m)
      order.list[[m]][[ft]] <- paste(tree$order, collapse = ",")
      order.list.names[[m]][[ft]] <- paste(paste("'cluster_", tree$order, "'", sep = ""), collapse = ",")
    }
  }
  average.number.col <- order.list[["average"]][["col"]]
  centroid.number.col <- order.list[["centroid"]][["col"]]
  complete.number.col <- order.list[["complete"]][["col"]]
  median.number.col <- order.list[["median"]][["col"]]
  mcquitty.number.col <- order.list[["mcquitty"]][["col"]]
  single.number.col <- order.list[["single"]][["col"]]
  wardd.number.col <- order.list[["ward.D"]][["col"]]
  wardd2.number.col <- order.list[["ward.D2"]][["col"]]
  
  average.number.row <- order.list[["average"]][["row"]]
  centroid.number.row <- order.list[["centroid"]][["row"]]
  complete.number.row <- order.list[["complete"]][["row"]]
  median.number.row <- order.list[["median"]][["row"]]
  mcquitty.number.row <- order.list[["mcquitty"]][["row"]]
  single.number.row <- order.list[["single"]][["row"]]
  wardd.number.row <- order.list[["ward.D"]][["row"]]
  wardd2.number.row <- order.list[["ward.D2"]][["row"]]
  
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
  
  
  ######################################################
  ## Export the table that will be read by D3 heatmap
  # heatmap.table.d3 <- "" ## Read from command-line arguments
  write.table(matches.tsv, file = matches.tsv.file, sep = "\t", quote = FALSE, row.names = FALSE)
  
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
  
  rgb.palette <- colorRampPalette(brewer.pal(9, "YlOrRd"), space="Lab")
  white <- "#FFFFFF"
  white <- append(white,rgb.palette(ceiling((max(count.matches.tab)/step+1))))
  
  if(count.mode == "presence"){
    colors <- colorRampPalette(c("#FFFFCC", "#800026"))
    white <- colors(2)
  }
  
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
  left <- (max(as.vector(sapply(c(seq.id, matrix.names), nchar))) + 2) * 10
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
                         "Centroid_col", centroid.number.col,
                         "Median_col", median.number.col,
                         "Mcquitty_col", mcquitty.number.col,
                         "Single_col", single.number.col,
                         "Ward_d_col", wardd.number.col,
                         "Ward_d2_col", wardd2.number.col,
                         
                         "Average_row", average.number.row,
                         "Complete_row", complete.number.row,
                         "Single_row", single.number.row,
                         "Centroid_row", centroid.number.row,
                         "Median_row", median.number.row,
                         "Mcquitty_row", mcquitty.number.row,
                         "Single_row", single.number.row,
                         "Ward_d_row", wardd.number.row,
                         "Ward_d2_row", wardd2.number.row,
                              
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
}