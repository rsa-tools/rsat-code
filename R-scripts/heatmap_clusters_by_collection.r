## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path
source(file.path(dir.rsat, 'R-scripts/config.R'))


## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")

## Load required libraries
## List of packages to install
required.packages = c("RColorBrewer",
                      "gplots",
                      "amap")

# List of RSAT-specific packages to be compiled on the server
for (pkg in c(required.packages)) {
  suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE, lib.loc=c(dir.rsat.rlib, .libPaths())))
}

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
heatmap.color.classes <- as.numeric(heatmap.color.classes)


################################################
## JSON + JavaScript code to add to HTML file ##
################################################

#venn.sortAreas(div, d);--return--
JSON.intersection = '
var sets_--nb-- = [
  {"sets": [0], "label": "--collection_A--", "size": --size_collection_A--},
  {"sets": [1], "label": "--collection_B--", "size": --size_collection_B--},
  {"sets": [0, 1], "size": --intersection_A_B--}]; '

Venn.diagram.set = '
<script>--return--
var chart = venn.VennDiagram()--return--
.width(275)--return--
.height(275);--return--
--return--
var div = d3.select("#--venn--")--return--
div.datum(--set--).call(chart);--return--
--return--
var tooltip = d3.select("body").append("div")--return--
.attr("class", "venntooltip");--return--
--return--
div.selectAll("path")--return--
.style("stroke-opacity", 0)--return--
.style("stroke", "#fff")--return--
.style("stroke-width", 0);--return--
--return--
div.selectAll("g")--return--
.on("mouseover", function(d, i) {--return--
  --return--
  tooltip.transition().duration(400).style("opacity", .9);--return--
  tooltip.text(d.size + " %");--return--
  --return--
  var selection = d3.select(this).transition("tooltip").duration(400);--return--
  selection.select("path")--return--
  .style("stroke-width", 3)--return--
  .style("fill-opacity", d.sets.length == 1 ? .4 : .1)--return--
  .style("stroke-opacity", 1);--return--
})--return--
--return--
.on("mousemove", function() {--return--
  tooltip.style("left", (d3.event.pageX) + "px")--return--
  .style("top", (d3.event.pageY - 28) + "px");--return--
})--return--
--return--
.on("mouseout", function(d, i) {--return--
  tooltip.transition().duration(400).style("opacity", 0);--return--
  var selection = d3.select(this).transition("tooltip").duration(400);--return--
  selection.select("path")--return--
  .style("stroke-width", 0)--return--
  .style("fill-opacity", d.sets.length == 1 ? .25 : .0)--return--
  .style("stroke-opacity", 0);--return--
});--return--
</script>--return--'

# cluster.counts.file <- "/home/jaimicore/Documents/PhD/clusters_summary_table.tab"

## Read cluster count table
clusters <- read.table(file = cluster.counts.file, sep = "\t", header = TRUE)
names(clusters) <- gsub("X.Cluster_ID", "Cluster_ID", names(clusters))
cluster.names.original <- as.vector(clusters$Cluster_ID)

# Read the root motif table and save the path to the logos
# root.motifs.table <- "/home/jaimicore/Documents/PhD/Multi_algorithms_analysis_hclust-average_Ncor0.4_cor0.6_root_motifs_table.tab"
root.motifs.files <- read.table(file = root.motifs.table, sep = "\t", header = TRUE)
names(root.motifs.files) <- gsub("X.Cluster_ID", "Cluster_ID", names(root.motifs.files))

## Create the arrays with the paths to the logos
## This images will be displayed in the heatmap with the number of clusters
cl.names <- as.vector(root.motifs.files$Cluster_ID)
logos <- as.vector(root.motifs.files$Logo)
logos.rc <- as.vector(root.motifs.files$Logo_RC)
pic.logos <- paste("pics['", cl.names, "'] = '", logos, "';", sep = "")
pic.logos.rc <- paste("pics_rc['", cl.names, "'] = '", logos.rc, "';", sep = "")

pic.logos <- paste(pic.logos, collapse = " ")
pic.logos.rc <- paste(pic.logos.rc, collapse = " ")

#################################
## Create the percentage table

## Count the number of motif per collection
nb.db <- dim(clusters)[2] - 2
motif.DB.counts <- apply(clusters[,3:(nb.db+2)], 2, sum)

percent.table <- NULL
coverage.contingency.table <- NULL
intersect.counter <- 0
Diagrams <- NULL
JSON.intersect.file <- paste(coverage.json.folder, "/intersection_data.json", sep = "")
JSON.intersect.file.rel <- paste(basename(coverage.json.folder), "/intersection_data.json", sep = "")
intersection.clusters <- list()
x <- sapply(names(motif.DB.counts), function(DB){
  
  ## Select those cluster with at least one motif corresponding
  ## to the current motifDB
  DB.motifs <- clusters[clusters[,DB] > 0,]
  
  ## Get collection A information
  collection.A.name <- DB
  # collection.A.size <- as.numeric(motif.DB.counts[DB])
  collection.A.nb.clusters <- length(DB.motifs[,DB])
  collection.A.size <- 100
  
  #################################################################################
  ## Calculate the overlap between the databases
  
  ## Select those cluster with at least one motif corresponding
  ## to the current motifDB
  coverage <- apply(DB.motifs[3:dim(DB.motifs)[2]], 2, sum) / motif.DB.counts
  
  sapply(names(coverage), function(n){
    
    intersect.counter <<- intersect.counter + 1
    intersection.clusters[[intersect.counter]] <<- list()
    
    ## Count the number of motifs in Collections A and B
    collection.A.nb.motifs <- sum(DB.motifs[,DB])
    collection.B.intersection <- sum(DB.motifs[,n])
    collection.B.nb.motifs <- motif.DB.counts[n]
    
    ## Get collection B information
    collection.B.name <- n
    collection.B.size <- as.numeric(motif.DB.counts[n])
    collection.B.nb.clusters <- length(which(DB.motifs[,n] > 0))
    collection.B.size <- 100
    
    ## Get the name of the clusters in Collection A and B
    collection.A.clusters.names <- as.vector(DB.motifs[,c("Cluster_ID",collection.A.name)]$Cluster_ID)
    
    t <- DB.motifs[,c("Cluster_ID",collection.B.name)]
    positive.index <- which(t[,2] > 0)
    collection.B.clusters.names <- as.vector(DB.motifs[positive.index,c("Cluster_ID",collection.B.name)]$Cluster_ID)
    
    intersect.size <- round(coverage[n] * 100)
    
    Venn.diagram.set.cp <- Venn.diagram.set
    set.nb <- paste("sets_", intersect.counter, sep = "")
    venn.nb <- paste("venn", intersect.counter, sep = "")
    Venn.diagram.set.cp <- gsub("--set--", set.nb, Venn.diagram.set.cp)
    Venn.diagram.set.cp <- gsub("--venn--", venn.nb, Venn.diagram.set.cp)
    Venn.diagram.set.cp <- gsub("\n", "", Venn.diagram.set.cp)
    Diagrams <<- append(Diagrams, Venn.diagram.set.cp)
    
    ## Fill the data for the Venn diagrams and export each JSON 
    ## in a different file
    JSON.intersection.cp <- JSON.intersection
    JSON.intersection.cp <- gsub("--collection_A--", collection.A.name, JSON.intersection.cp)
    JSON.intersection.cp <- gsub("--collection_B--", collection.B.name, JSON.intersection.cp)
    JSON.intersection.cp <- gsub("--size_collection_A--", collection.A.size, JSON.intersection.cp)
    JSON.intersection.cp <- gsub("--size_collection_B--", collection.B.size, JSON.intersection.cp)
    JSON.intersection.cp <- gsub("--intersection_A_B--", intersect.size, JSON.intersection.cp)
    JSON.intersection.cp <- gsub("--nb--", intersect.counter, JSON.intersection.cp)
    
    intersection.clusters[[intersect.counter]][["Collections"]] <<- paste(unique(c(collection.A.name, collection.B.name)), collapse = "<br>")
    int.cl <- intersect(collection.A.clusters.names, collection.B.clusters.names)
    int.cl <- paste(int.cl, collapse = "<br>")
    intersection.clusters[[intersect.counter]][["Clusters"]] <<- int.cl
    intersection.clusters[[intersect.counter]][["Sentence"]] <<- paste(collection.A.name, " covers ", round(coverage[n] * 100, digits = 2), "% of ", collection.B.name, sep = "")
    intersection.clusters[[intersect.counter]][["CollectionA_size"]] <<- collection.A.nb.motifs
    intersection.clusters[[intersect.counter]][["CollectionB_size"]] <<- collection.B.nb.motifs
    intersection.clusters[[intersect.counter]][["CollectionB_intersection"]] <<- collection.B.intersection
    intersection.clusters[[intersect.counter]][["CollectionA_name"]] <<- DB
    intersection.clusters[[intersect.counter]][["CollectionB_name"]] <<- n
    
    if(intersect.counter == 1){
      file.remove(JSON.intersect.file, showWarnings = FALSE)
    }
    write(JSON.intersection.cp, file = JSON.intersect.file, append = TRUE)
        
  })
  
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

Diagrams <- paste(Diagrams, collapse = "")

coverage.counter <- 0
coverage.pics.buttons <- NULL
thrash <- sapply(names(motif.DB.counts), function(x){
  sapply(names(motif.DB.counts), function(y){

    coverage.counter <<- coverage.counter + 1
    
    coverage.pics.buttons <<- append(coverage.pics.buttons, paste("<div class='coverage_button button_click' id='d_", coverage.counter,"_link'> <strong>", x, " vs ", y,"</strong></div>", sep = ""))
    
    if(coverage.counter %% length(names(motif.DB.counts)) == 0){
      coverage.pics.buttons <<- append(coverage.pics.buttons, paste("--return--", sep = ""))
    }
  })
})
coverage.pics.buttons <- as.vector(coverage.pics.buttons)
coverage.pics.buttons <- paste(coverage.pics.buttons, collapse = "")


coverage.pics <- sapply(1:(length(names(motif.DB.counts)) ^2), function(x){
    coverage.counter <<- coverage.counter + 1
    
    paste("<div id='d_", x, "' class='coverage_pic' style='display:none;position:relative;float:left;font-size:8px;'><p style='text-align:center;padding-top:5px;' class='mono'><strong>Venn diagram</strong></p><div style='padding-left:175px' id='venn", x, "'></div><p style='text-align:center;padding-top:5px;' class='mono'>", intersection.clusters[[x]][["Sentence"]], "</p><table style='width:650px;text-align:left;' class='mono'><thead><tr><th>Collections</th><th>", intersection.clusters[[x]][["CollectionA_name"]]," size</th><th>", intersection.clusters[[x]][["CollectionB_name"]], "<br>(intersection)</th><th>", intersection.clusters[[x]][["CollectionB_name"]]," size</th><th>Intersection</th></tr></thead><tbody><td>", intersection.clusters[[x]][["Collections"]], "</td><td>", intersection.clusters[[x]][["CollectionA_size"]], "</td><td>", intersection.clusters[[x]][["CollectionB_intersection"]], "</td><td>", intersection.clusters[[x]][["CollectionB_size"]], "</td><td>", intersection.clusters[[x]][["Clusters"]],"</td></tbody></table></div> --return--", sep = "")

})
coverage.pics <- as.vector(coverage.pics)
coverage.pics <- paste(coverage.pics, collapse = "")
  
hide.show.coverage.pics <- sapply(1:(length(names(motif.DB.counts)) ^2), function(x){
  
  paste("$(document).ready(function() { $('#d_", x,"_link').click(function() { $('.coverage_button').removeClass('selected_coverage_button'); $('.coverage_pic').hide(); $('#d_", x, "').show(); $(this).toggleClass('selected_coverage_button') }); }); --return--", sep = "")
})
hide.show.coverage.pics <- as.vector(hide.show.coverage.pics)
hide.show.coverage.pics <- paste(hide.show.coverage.pics, collapse = "")

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
coverage.contingency.table <- coverage.contingency.table
colnames(coverage.contingency.table) <- names(motif.DB.counts)
rownames(coverage.contingency.table) <- names(motif.DB.counts)
write.table(coverage.contingency.table, file = coverage.table.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#######################################################
## Run the hierarchical clustering with four methods
## (average + complete + single + ward). 
## Save the order of the nodes
comp.order.list.columns <- list()
comp.order.list.rows <- list()

sapply(c("average", "complete", "single", "ward"), function(m){
  
  if(m == "ward"){
    temp <- m
    m <- "ward.D"
  }
  
  
  pfile <- paste(coverage.json.folder, "/coverage_clustering_", m,".json", sep = "")
  pdf(file = pfile)
  hm.collections <- heatmap.2(coverage.contingency.table,
                           hclustfun = function(x) hclust(x,method = m),
                           distfun = function(x) Dist(x,method = 'pearson')
  )
  dev.off()
  
  if(m == "ward.D"){
    m <- "ward"
  }
  
  comp.order.list.rows[[m]] <<- paste(rev(hm.collections[[1]]), collapse = ",")
  comp.order.list.columns[[m]] <<- paste(rev(hm.collections[[2]]), collapse = ",")
})

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
col.nb <- dim(coverage.contingency.table)[1]
row.nb <- dim(coverage.contingency.table)[2]
default.labels <- paste(paste("'", names(motif.DB.counts), "'", sep = ""), collapse = ",")
default.number <- paste(1:length(motif.DB.counts), collapse = ",")
left <- (max(as.vector(sapply(names(motif.DB.counts), nchar))) + 2) * 10
cell.size <- 20
bottom <- 80
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

## Save the order of the columns of the complementarity heatmap
comp.average.c.number <- comp.order.list.columns[["average"]]
comp.complete.c.number <- comp.order.list.columns[["complete"]]
comp.single.c.number <- comp.order.list.columns[["single"]]
comp.ward.c.number <- comp.order.list.columns[["ward"]]

## Save the order of the rows of the complementarity heatmap
comp.average.r.number <- comp.order.list.rows[["average"]]
comp.complete.r.number <- comp.order.list.rows[["complete"]]
comp.single.r.number <- comp.order.list.rows[["single"]]
comp.ward.r.number <- comp.order.list.rows[["ward"]]

heatmap.width <- 300 + (cell.size * col.nb)

coverage.info <- matrix(c("Collection_labels", default.labels,
                          "Collection_number", default.number,
                          "Left_space", left,
                          "Bottom_space", bottom,
                          "Col_number", col.nb,
                          "Row_number", row.nb,
                          "Legend_Head", legend.header,
                          "Average_r_number_comp", comp.average.r.number,
                          "Complete_r_number_comp", comp.complete.r.number,
                          "Single_r_number_comp", comp.single.r.number,
                          "Ward_r_number_comp", comp.ward.r.number,
                          "Average_c_number_comp", comp.average.c.number,
                          "Complete_c_number_comp", comp.complete.c.number,
                          "Single_c_number_comp", comp.single.c.number,
                          "Ward_c_number_comp", comp.ward.c.number,
                          "Coverage_pics", coverage.pics,
                          "Coverage_pics_buttons", coverage.pics.buttons,
                          "Hide_show_coverage_pics", hide.show.coverage.pics,
                          "Diagrams", Diagrams,
                          "Coverage_JSON", JSON.intersect.file.rel,
                          "Heatmap_width", heatmap.width
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
rgb.palette <- colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")
white <- "#FFFFFF"
white <- append(white,rgb.palette(ceiling((max(clusters)/step))))

###############################################################
## Run the hierarchical clustering with the three methods
## (average + complete + single). Save the order of the nodes
order.list.rows <- list()
order.list.columns <- list()
order.list.names <- list()

sapply(c("average", "complete", "single", "ward"), function(m){
  
  if(m == "ward"){
    temp <- m
    m <- "ward.D"
  }
  
  pfile <- paste(coverage.json.folder, "/collection_clustering_", m,".json", sep = "")
  pdf(file = pfile)


  hm.clusters <- heatmap.2(clusters.matrix,
                           hclustfun = function(x) hclust(x,method = m),
                           distfun = function(x) Dist(x,method = 'pearson')
  )
  dev.off()
  
  if(m == "ward.D"){
    m <- "ward"
  }
  
  order.list.rows[[m]] <<- paste(rev(hm.clusters[[1]]), collapse = ",")
  cluster.order <- as.numeric(unlist(strsplit(order.list.rows[[m]], ",")))
  order.list.columns[[m]] <<- paste(rev(hm.clusters[[2]]), collapse = ",")
  order.list.names[[m]] <<- paste(paste("'", cluster.names.original[cluster.order], "'", sep = ""), collapse = ",")
})

###############################################
## Parse the Heatmap table format used in D3
## This table is printed in a new file
x <- data.frame(t(clusters))
names(x) <- cluster.names.original
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

cluster.names <- paste(paste("'", cluster.names.original, "'", sep =""), collapse=",")

## Get the clusters number orderer according the linkage method
cluster.number <- paste(1:dim(clusters)[1], collapse=",")
average.r.number <- order.list.rows[["average"]]
complete.r.number <- order.list.rows[["complete"]]
single.r.number <- order.list.rows[["single"]]
ward.r.number <- order.list.rows[["ward"]]

average.c.number <- order.list.columns[["average"]]
complete.c.number <- order.list.columns[["complete"]]
single.c.number <- order.list.columns[["single"]]
ward.c.number <- order.list.columns[["ward"]]

## Default names
default.names <- paste(paste("'", cluster.names.original, "'", sep = ""), collapse = ",")
default.number <- paste(1:dim(clusters)[1], collapse = ",")

## Heatmap variables
col.nb <- dim(clusters)[1]
row.nb <- dim(clusters)[2]

## Row
heatmap.rows.nb <- paste(1:row.nb, collapse=",")
heatmap.rows.name <- paste(paste("'", collection.names, "'", sep = ""), collapse = ",")

## Color scale
color.scale <- append("#FFFFFF",rgb.palette(21))
color.scale <- paste("'", color.scale, "'",collapse=",")

## Collections
collections <- paste(paste("'", collection.names, "'", sep = ""), collapse = ",")

## Range to color the values
domain.nb <- seq(from = 1, to = max(clusters), by = step)

domain <- paste(domain.nb, collapse=",")

## Legend
# legend <- c(0,domain.nb)
legend <- 0
legend <- append(legend, seq(from = 1, to = max(clusters), by = step))
legend <- paste(legend, collapse=",")

## Right space
left <- (max(as.vector(sapply(c(collection.names, cluster.names.original), nchar))) + 2.5) * 10

## Right space
top <- (max(as.vector(sapply(collection.names, nchar))) + 2.5) * 10

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

html.body.size <- 200 + left + (col.nb*cell.size) + 30

order.info <- matrix(c("Gradient", gradient,
                       "Cluster_names", cluster.names,
                       "Cluster_number", cluster.number,
                       "Average_c_number", average.c.number,
                       "Complete_c_number", complete.c.number,
                       "Single_c_number", single.c.number,
                       "Ward_c_number", ward.c.number,
                       "Average_r_number", average.r.number,
                       "Complete_r_number", complete.r.number,
                       "Single_r_number", single.r.number,
                       "Ward_r_number", ward.r.number,
                       "Cell_size", cell.size,
                       "Col_number", col.nb,
                       "Row_number", row.nb,
                       "Row_order_default", heatmap.rows.nb,
                       "Domain", domain,
                       "Legend", legend,
                       "Legend_Head", legend.header,
                       "Left_space", left,
                       "Top_space", top,
                       "Bottom_space", left,
                       "Body", html.body.size,
                       "Collections", collections,
                       "Logos", pic.logos,
                       "Logos_RC", pic.logos.rc,
                       "Color_scale", color.scale
), nrow = 2)
order.info.df <- t(data.frame(order.info))
verbose(paste("Exporting table with the order of the clusters (Required in D3 Heatmap)", order.list.file), 2)
write.table(order.info.df, file = attributes.list.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# tfs <- c("DSP1", "Ez", "GAF", "PHOL", "PHO", "PH", "PSQ", "SPPS")
# paste("/home/mentrevan/polycomb_clusters/results/TF-and-PolycombSU_summits/", tfs,"_replicate2_summit_sorted_by_scores_pm300/", tfs, "_replicate2_summit_sorted_by_scores_pm300_peak-motifs_top1000peaks/results/discovered_motifs/", tfs, "_replicate2_summit_sorted_by_scores_pm300_top1000peaks_motifs_discovered.tf", sep = "")