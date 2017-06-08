## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

## Load some libraries
source(file.path(dir.rsat, 'R-scripts/config.R'))

## Load required libraries
required.packages = c("RColorBrewer","gplots")

## List of RSAT-specific packages to be compiled on the server
for (pkg in c(required.packages)) { #required.packages.bioconductor
  suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE, lib.loc=c(dir.rsat.rlib, .libPaths())))
}

###########################################
## Read arguments from the command line.
##
## Arguments passed on the command line
## will over-write the default arguments
## specified above.
# message("Reading arguments from command-line")
args <- commandArgs(trailingOnly=TRUE)
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
} else if (!exists("maxNWD.table.file")) {
  stop("Missing mandatory argument ( NWD table (matrix-quality output) ): maxNWD.table.file ")
} else if (!exists("html.template.file")) {
  stop("Missing mandatory argument (The path to the RSAT template for the heatmap): html.template.file ")
} else if (!exists("d3.base")) {
  stop("Missing mandatory argument (The path to the D3 source code):d3.base  ")
} else if (!exists("d3.array.base")) {
  stop("Missing mandatory argument (The path to the D3 Array library source code):d3.array.base  ")
} else if (!exists("maxNWD.tsv")) {
  stop("Missing mandatory argument (The path to the tsv file used as input for the dynamic heatmap): maxNWD_tsv ")
} else if (!exists("maxNWD.heatmap.html")) {
  stop("Missing mandatory argument (The path to the D3 Dynamic heatmap): maxNWD.heatmap.html  ")
} else if (!exists("Diff.maxNWD.tsv")) {
  stop("Missing mandatory argument (The path to the tsv file used as input for the dynamic heatmap): Diff.maxNWD_tsv ")
} else if (!exists("Diff.maxNWD.heatmap.html")) {
  stop("Missing mandatory argument (The path to the D3 Dynamic heatmap): Diff.maxNWD.heatmap.html  ")
} 
if (!exists("heatmap.color.palette")) {
  heatmap.color.palette <- "RdBu";
}
if (!exists("heatmap.color.classes")) {
  heatmap.color.classes <- as.numeric(9);
}
heatmap.color.classes <- as.numeric(heatmap.color.classes)


# cat /home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap/matrix-enrichment_heatmap.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " prefix = '/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap/second_test_auto'; maxNWD.table.file = '/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Motif_Enrichment_all_nwd_plot/maxNWDsignificantScore_heatmap_compare.txt'; html.template.file = 'motif_enrichment_dynamic_heatmap_d3.html'; maxNWD.tsv = '/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap/second_test_auto_matrix_heatmap.tsv'; maxNWD.heatmap.html = '/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap/second_test_auto_motif_enrichment_maxNWD_heatmap.html'; d3.base = '/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap/d3.v3.min.js'; d3.array.base = '/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap/d3-array.v0.6.min.js'" 

# prefix <- "test_motif_enrichment"
# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/t/temp/Human_motifs_Epromoters_vs_Inactive_Promoters_2/Dynamic_Heatmap")
# maxNWD.table.file <- "/home/jaimicore/Dropbox/matrix-clustering_article/NFKB_ChIP-seq/results/matrix-enrichment/NFKB_vs_Random_fragments/NFKB_vs_Random_fragments_all_nwd_plot/maxNWDsignificantScore_heatmap_compare.txt"
# html.template.file <- "motif_enrichment_dynamic_heatmap_d3.html"

base.name <- basename(prefix)

#######################################
## Read the input file: maxNWD table
max.NWD.table <- read.table(maxNWD.table.file, sep = "\t", header = TRUE)
max.NWD.table <- round(max.NWD.table, digits = 3)
max.NWD.table[is.na(max.NWD.table)] <- 0

## Calculate Differential
sets <- colnames(max.NWD.table)
for(n in 1:(length(sets)-1)){
  
  for(m in 2:(length(sets))){
    ## Get the set names
    set1 <- sets[n]
    set2 <- sets[m]
    
    ## Skip when the both sets are the same
    if(n < m){
      ## Assign the names to the differential
      diff.name.1 <- paste("Diff", set1, set2, sep = "__")
      diff.name.2 <- paste("Diff", set2, set1, sep = "__")
      
      print(diff.name.1)
      print(diff.name.2)
      
      max.NWD.table[,diff.name.1] <- max.NWD.table[,set1] - max.NWD.table[,set2]
      max.NWD.table[,diff.name.2] <- max.NWD.table[,diff.name.1] * -1
    }
  }
}

nb.sets <- length(sets)
nb.diff.columns <- dim(max.NWD.table)[2]

NWD.sub <- max.NWD.table[,(nb.sets+1):nb.diff.columns]
NWD.sub <- NWD.sub[order(NWD.sub[,1], decreasing = TRUE),]
NWD.sub <- as.matrix(round(NWD.sub, digits = 2))

#######################################################################
## Convert the DataFrame in a 'tsv' object which is the input format
## for D3 heatmap.
for(set in c("Normal", "Diff")){
  
  tsv.tab <- NULL
  
  if(set == "Normal"){
    tab <- max.NWD.table
    file.export <- maxNWD.tsv
    heatmap.html <- maxNWD.heatmap.html
  } else if (set == "Diff"){
    tab <- NWD.sub
    file.export <- Diff.maxNWD.tsv
    heatmap.html <- Diff.maxNWD.heatmap.html
  }
  
  for(j in 1:dim(tab)[1]){
    for(i in 1:dim(tab)[2]){
      tsv.tab <<- rbind(tsv.tab, matrix(c(j,i, as.numeric(tab[j,i])), nrow = 1))
    }
  }
  
  ## Export the TSV table
  colnames(tsv.tab) <- c("Row", "Col", "Value")
  write.table(tsv.tab, file = file.export, sep = "\t", quote = FALSE, row.names = FALSE)
  
  ###############################################################
  ## Open and modify the Heatmap D3 template with the new data
  html.report <- readLines(html.template.file)
  
  ## Get the Sequences (column) and Motifs (row) names
  sequences.names <- colnames(tab)
  motifs.names <- rownames(tab)
  
  ## Add the logo path to the template
  logo.path <- sapply(motifs.names, function(m){
    paste("'", m, "/", base.name, "_", m, "_logo_m1.png'", sep = "")
  })
  logo.path <- as.vector(logo.path)
  logo.path <- paste(logo.path, collapse = ",")
  html.report <- gsub("--logo_path--", logo.path, html.report)
  
  ## Add the link to the distrib comparison curves pictures
  curves.path <- sapply(motifs.names, function(m){
    paste("'", m, "/", base.name, "_", m, "_score_distrib_compa_logy.png'", sep = "")
  })
  curves.path <- as.vector(curves.path)
  curves.path <- paste(curves.path, collapse = ",")
  html.report <- gsub("--curves_path--", curves.path, html.report)
  
  
  ## Insert the number of rows and columns
  col.nb <- length(sequences.names)
  row.nb <- length(motifs.names)
  html.report <- gsub("--r_numb--", row.nb, html.report)
  html.report <- gsub("--c_numb--", col.nb, html.report)
  
  ## Insert the D3 paths
  html.report <- gsub("--d3_base--", d3.base, html.report)
  html.report <- gsub("--d3_array_base--", d3.array.base, html.report)
  
  ## Insert the Default order of the rows and columns
  seq.number <- paste(1:col.nb, collapse = ",")
  motif.number <- paste(1:row.nb, collapse = ",")
  html.report <- gsub("--col_order_default--", seq.number, html.report)
  html.report <- gsub("--row_order_default--", motif.number, html.report)
  
  ## Insert the Column and Row names
  sequences.names.cat <- paste("'" , sequences.names, "'", sep = "")
  sequences.names.cat <- paste(sequences.names.cat, collapse = ", ")
  motifs.names.cat <- paste("'" , motifs.names, "'", sep = "")
  motifs.names.cat <- paste(motifs.names.cat, collapse = ", ")
  html.report <- gsub("--col_names--", sequences.names.cat, html.report)
  html.report <- gsub("--row_names--", motifs.names.cat, html.report)
  
  ## Insert the TSV file
  html.report <- gsub("--file--", shortpath.maxNWD.tsv, html.report)
  
  ## Div bottom + Cell size
  cell.size <- 30
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
    cell.size <- 20
    legend.header <- bottom - 27
  } else if(row.nb < 18){
    bottom <- 270
    cell.size <- 15
    legend.header <- bottom - 27
  }
  html.report <- gsub("--cell_size--", cell.size, html.report)
  
  ## Left space
  left <- (max(as.vector(sapply(c(sequences.names, motifs.names), nchar))) + 2.5) * 10
  html.report <- gsub("--left--", left, html.report)
  
  ## D3 path
  D3 <- "d3js.org/d3.v3.min.js"
  html.report <- gsub("--d3--", D3, html.report)
  
  ## Body size
  html.body.size <- 250 + left + (col.nb*cell.size)
  html.report <- gsub("--body--", html.body.size, html.report)
  
  ## Calculate the legend names for the color scale
  # stp <- (max(max.NWD.table) - min(max.NWD.table))/9
  # legend.domain.values <- seq(from = min(max.NWD.table), to = max(max.NWD.table), by = 0.05)
  # limit <- max(abs(c(min(max.NWD.table), max(max.NWD.table))))
  # limit <- round(limit, digits = 1)
  # legend.domain.values <- seq( (min(max.NWD.table) - 0.2), limit, by = 0.05)
  
  ## Remove noisy extremities
  tab.range <- range(tab)
  tab.range.abs <- abs(tab.range)
  
  if (tab.range.abs[1] > tab.range.abs[2]){
    tab.range[2] <- -tab.range[1]
  } else {
    tab.range[1] <- -tab.range[2]
  }
  
  legend.domain.values <- quantile(tab.range, probs = seq(0, 1, 0.1))
  legend.length <- length(legend.domain.values)
  legend <- legend.domain.values
  legend <- round(legend, digits = 3)
  legend <- paste(rev(legend), collapse = ",")
  html.report <- gsub("--data_legend--", legend, html.report)
  
  ## Create Gradient Hexadecimal:
  ## Given X hexa colors creates a color
  palette.hexa <- colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")(11)
  palette.hexa <- rev(palette.hexa)
  
  palette <- paste("'" , palette.hexa, "'", sep = "")
  palette <- paste(palette, collapse = ", ")
  
  palette.rev <- paste("'" , rev(palette.hexa), "'", sep = "")
  palette.rev <- paste(palette.rev, collapse = ", ")
  
  html.report <- gsub("--gradient--", palette, html.report)
  html.report <- gsub("--gradient_rev--", palette.rev, html.report)
  
  ## Color domain
  domain <- legend.domain.values[1:(legend.length-1)]
  domain <- paste(domain, collapse = ",")
  html.report <- gsub("--domain--", domain, html.report)
  
  ## Export the report
  # maxNWD.heatmap.html <- paste(prefix, "_motif_enrichment_maxNWD_heatmap.html", sep = "")
  write(html.report, file = heatmap.html)

}
