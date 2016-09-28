## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}
dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")

## Load some libraries
source(file.path(dir.rsat, 'R-scripts/config.R'))

## Load required libraries
required.packages = c("IRanges",
                      "RColorBrewer",
                      "gplots",
                      "jpeg",
                      "amap",
                      "dynamicTreeCut",
                      "qvalue")

## List of RSAT-specific packages to be compiled on the server
for (pkg in c(required.packages)) { #required.packages.bioconductor
  suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE, lib.loc=c(dir.rsat.rlib, .libPaths())))
}


#################################################################################################
## Functions
create.html.tab <- function(tab, img = 0, plot = 0, link.text.covered = 0, link.text.not.covered = 0){
  
  full.tab <- NULL
  head.tab <- "<div id='individual_motif_tab' style='width:1500px;display:none' class='tab div_chart_sp'><p style='font-size:12px;padding:0px;border:0px'><b>Individual Motif View</b></p><table id='Motif_tab' class='hover compact stripe' cellspacing='0' width='1190px' style='padding:15px;align:center;'><thead><tr><th class=\"tab_col\"> Motif_ID </th><th class=\"tab_col\"> Motif_name </th> <th class=\"tab_col\"> P-value </th> <th class=\"tab_col\"> E-value </th> <th class=\"tab_col\"> Significance </th> <th class=\"tab_col\"> FDR </th> <th class=\"tab_col\"> Nb of hits </th><th class=\"tab_col\"> Nb of sequences </th><th class=\"tab_col\">Fraction of sequences</th><th class=\"tab_col\"> Chi-squared</th><th class=\"tab_col\">Profile cluster</th> <th class=\"tab_col\"> Profile </th> <th class=\"tab_col\"> TFBSs </th><th class=\"tab_col\"> TFBSs per seq </th> <th class=\"tab_col\"> Logo </th> <th class=\"tab_col\"> Logo (RC) </th> <th class=\"tab_col\"> Covered sequences </th> <th class=\"tab_col\"> Not Covered sequences </th> </tr></thead><tbody>"
  content.tab <- apply(tab, 1, function(row){
    
    row.length <- length(row)
    rows.nb <- 1:row.length
    
    ## Get the number of the columns with/without picture or plot
    ## This is done because the tab require different arguments
    rows.simple <- rows.nb[!(rows.nb %in% img)]
    rows.simple <- rows.simple[!(rows.simple %in% plot)]
    rows.simple <- rows.simple[!(rows.simple %in% link.text.covered)]
    rows.simple <- rows.simple[!(rows.simple %in% link.text.not.covered)]
    
    rows.pic <- rows.nb[rows.nb %in% img]
    rows.plot <- rows.nb[rows.nb %in% plot]
    rows.text.link.covered <- rows.nb[rows.nb %in% link.text.covered]
    rows.text.link.not.covered <- rows.nb[rows.nb %in% link.text.not.covered]
    
    ## Columns with simple text
    rows.text <- paste("<td>", row[rows.simple], "</td>", collapse = "")
    
    ## Columns with images
    rows.pic.text <- paste("<td><img class='logo_tab' src ='", as.character(row[rows.pic]), "'/></td>", collapse = "")
    
    rows.text.link.covered <- paste("<td><a href='", as.character(row[rows.text.link.covered]), "' target='_blank'>Covered</a></td>", collapse = "")
    rows.text.link.not.covered <- paste("<td><a href='", as.character(row[rows.text.link.not.covered]), "' target='_blank'>Not Covered</a></td>", collapse = "")
    
    ## Columns with plots and links
    rows.plot.pdf <- sapply(row[rows.plot], function(x){
      gsub("jpeg","pdf", x)
    })
    rows.plot.text <- paste("<td><a href='", rows.plot.pdf, "' target='_blank'><img class='plot_tab' src ='", row[rows.plot], "'/></a></td>", collapse = "")

    ## Head and tail tags
    row.head <- "<tr>"
    row.tail <- "</tr>"
    paste(row.head, rows.text, rows.plot.text, rows.pic.text, rows.text.link.covered, rows.text.link.not.covered, row.tail, sep = "")    

  })
  
  tail.tab <- "</tbody></table></div>"
  
  full.tab <- c(head.tab, content.tab , tail.tab)
  return(full.tab) 
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

# message("Checking mandatory arguments")
if (!exists("matrix.scan.file")) {
  stop("Missing mandatory argument (matrix-scan results table): matrix.scan.file ")
} else if (!exists("sequence.names.file")) {
  stop("Missing mandatory argument (sequence names table): sequence.names.file ")
} else if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
} else if (!exists("ID.to.names.correspondence.tab")) {
  stop("Missing mandatory argument (Correspondence of the motif IDs and names): ID.to.names.correspondence.tab ")
} else if (!exists("seq.length")) {
  
  ## For the moment we assume all the input sequences have the same size
  stop("Missing mandatory argument (Sequence length): seq.length ")
} else if (!exists("html.template.file")){
  stop("Missing mandatory argument (HTML template to draw the profiles): html.template.file ")
} else if (!exists("results.folder")){
  stop("Missing mandatory argument (Working directory): results.folder ")
} else if (!exists("basename")){
  stop("Missing mandatory argument (Basename): basename ")
}

if (!exists("p.val")) {
  p.val <- 1e-3
} 
if (!exists("bin")) {
  bin <- 25
}
if (!exists("draw.area")) {
  draw.area <- 0
}
if (!exists("off.set")) {
  off.set <- as.numeric(0)
}
if (!exists("off.set.type")) {
  off.set.type <- as.character("")
}
if (!exists("individual.plots")) {
  individual.plots <- 0
}
if (!exists("logo.folder")) {
  logo.folder <- paste(basename, "_logos", sep = "")
}
if (!exists("individual.plots")) {
  individual.plots <- 0
}
if (!exists("heatmap.dendo")) {
  heatmap.dendo <- "show"
}
if (!exists("heatmap.color.palette")) {
  heatmap.color.palette <- "RdBu";
}
if (!exists("heatmap.color.classes")) {
  heatmap.color.classes <- as.numeric(11);
}
heatmap.color.classes <- as.numeric(heatmap.color.classes)

## Heatmap dendogram position
if (heatmap.dendo == "show"){
  heatmap.dendo <- "row"
} else if(heatmap.dendo == "hide"){
  heatmap.dendo <- "none"
}
print(heatmap.dendo)

## Create a file to store the resulting tables
covered.tables.dir <- paste(prefix, "_covered_sequences_info", sep = "")
dir.create(covered.tables.dir, showWarnings = FALSE)

# matrix.scan.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/mkv_1/Jun_Chip_seq_bin_size_25_pval1e-3_mkv_1_matrix_scan_results_PARSED.tab"
# prefix <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/mkv_1/Jun_Chip_seq_bin_size_25_pval1e-3_mkv_1"
# ID.to.names.correspondence.tab <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/mkv_1/Jun_Chip_seq_bin_size_25_pval1e-3_mkv_1_TF_ID_name_correspondence.tab"
# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/mkv_1/")
# sequence.names.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Jun_Chip_seq_bin_size_25_pval1e-3_mkv_1_matrix_scan_sequence_names.tab"

# matrix.scan.file <- "/home/jaimicore/test/Clustered_vs_single_PSSMs_JUN_FOS/Clustered_vs_single_PSSMs_JUN_FOS_matrix_scan_results_PARSED.tab"
# prefix <- "/home/jaimicore/test/Clustered_vs_single_PSSMs_JUN_FOS/Clustered_vs_single_PSSMs_JUN_FOS"
# ID.to.names.correspondence.tab <- "/home/jaimicore/test/Clustered_vs_single_PSSMs_JUN_FOS/Clustered_vs_single_PSSMs_JUN_FOS_TF_ID_name_correspondence.tab"
# setwd("/home/jaimicore/test/Clustered_vs_single_PSSMs_JUN_FOS/")
# sequence.names.file <- "/home/jaimicore/test/Clustered_vs_single_PSSMs_JUN_FOS/Clustered_vs_single_PSSMs_JUN_FOS_matrix_scan_sequence_names.tab"

# matrix.scan.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Epromoters/K562_bin_size_25_pval1e-3_matrix_scan_results_PARSED.tab"
# prefix <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Epromoters/K562_bin_size_25_pval1e-3"
# ID.to.names.correspondence.tab <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Epromoters/K562_bin_size_25_pval1e-3_TF_ID_name_correspondence.tab"
# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Epromoters/")

##############################
## Read matrix-scan table 1
verbose(paste("Reading matrix-scan results table"), 1)
matrix.scan.results <- read.csv(file = matrix.scan.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results) <- c("seq_id", "ft_name", "bspos", "Pval")

##############################
## Read sequence names table
verbose(paste("Reading sequence names table"), 1)
sequence.names.tab <- read.csv(file = sequence.names.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(sequence.names.tab) <- c("seq_id")
total.scanned.sequences <- length(as.vector(sequence.names.tab$seq_id))
scanned.sequences <- unique(as.vector(sequence.names.tab$seq_id))

######################################
## Create the column -log10(pvalue)
## Assign a class to each p-value
matrix.scan.results$Pval.minlog10 <- -log10(matrix.scan.results$Pval)
matrix.scan.results$Pval.class <- ceiling(matrix.scan.results$Pval.minlog10*2)/2

classes.pval <- sort(unique(matrix.scan.results$Pval.class))
classes.pval.letters <- LETTERS[1:length(classes.pval)]

matrix.scan.results$Pval.class.letter <- sapply(matrix.scan.results$Pval.class, function(x){
  p.class <- which(classes.pval == x)
  classes.pval.letters[p.class ]
})

min.pval.minus.log10 <- min(matrix.scan.results$Pval.minlog10)
max.pval.minus.log10 <- max(matrix.scan.results$Pval.minlog10)

#################
## Set p-value
p.val <- as.numeric(p.val)

#############################
## Get the sequences + motif name's
seq.id <- unique(as.vector(matrix.scan.results$seq_id))
matrix.names <- unique(as.vector(matrix.scan.results$ft_name))

###################################
## Calculate the sequence limits
seq.length <- as.numeric(seq.length)

###################################
## Divide the sequences in bins 
# bin <- 25
verbose(paste("Setting bins of size", bin), 1)
bin <- as.numeric(bin)

## If the sequences are centered in the peak summit
if(off.set == 0){
  
  limits <- seq.length/2
  limit.dw <- limits
  limit.up <- -limits
  
  windows <- IRanges(start = seq(from = limit.up, to = limit.dw - bin + 1, by = bin), width = bin)
  
  ## Calculate the shift in order to set the ranges relative to the center of the sequence
  ## or to the user specified position
  windows.min.value <- min(abs(min(windows)))
  windows <- shift(windows, windows.min.value)
  
} else {
  
  if (off.set.type == "start"){
    
    windows <- IRanges(start = seq(from = seq.length, to = 0 - bin + 1, by = bin), width = bin)
    windows <- shift(windows, seq.length - off.set)
    
    limit.dw <- seq.length - off.set
    limit.up <- -off.set
    
  } else if (off.set.type == "end"){
    
    windows <- IRanges(start = seq(from = seq.length, to = 0 - bin + 1, by = bin), width = bin)
    windows <- shift(windows, off.set)
    
    limit.dw <- seq.length - off.set
    limit.up <- off.set
  }
}

## Adapt the original BS position realtive to the limits 
## calculated in the step before
matrix.scan.results$bspos <- matrix.scan.results$bspos + limits

ID.names.tab <- ID.to.names.correspondence.tab
ID.names <- read.table(ID.names.tab, sep = "\t")

feature.attributes <- vector("list", length(matrix.names)) 
profiles <- NULL
windows.labels <- NULL

setwd(results.folder)

print.formats <- c("pdf", "jpeg")


##########################################################
## Plot the distribution of TFBSs at different p-values ##
##########################################################

## Assign a color to each p-value class

## The sequencial color palette has a maximum of 9 colors
nb.color.classes <- length(classes.pval.letters)
if(length(classes.pval.letters) > 9){
  nb.color.classes <- 9
}
pval.class.colors <- colorRampPalette(brewer.pal(nb.color.classes, "YlGnBu"), space="Lab")(length(classes.pval.letters))

classes.to.colors <- list()
for(x in 1:length(classes.pval.letters)){
  classes.to.colors[[classes.pval.letters[x]]] <- pval.class.colors[x]
}

## Create directory with the TFBSs distribution
dir.create(paste(basename(prefix), "_TFBSs_pval_distribution/", sep = ""), showWarnings = FALSE, recursive = TRUE)
dir.create(paste(basename(prefix), "_TFBSs_per_seq/", sep = ""), showWarnings = FALSE, recursive = TRUE)
verbose(paste("Creating plots with distribution of TFBSs at different p-values"), 1)
covered.sequences.per.motif <- list()

thr <- sapply(1:length(matrix.names), function(m){
  
  ## Get the matrix name
  matrix.query <- matrix.names[m]
  matrix.query.name <- as.vector(ID.names[,2][which(ID.names[,1] == matrix.query)][1])
  
  # print(matrix.query)

  ## Get the sub-table with the hits of the query matrix
  matrix.query.selection <- matrix.scan.results[matrix.scan.results$ft_name == matrix.query,]
  matrix.query.classes <- sort(unique(matrix.query.selection$Pval.class.letter))
  
  ##
  nb.hits.per.sequence <- table(as.vector(matrix.query.selection$seq_id))
  nb.hits.per.sequence.range <- range(nb.hits.per.sequence)
  min.nb.hits <- nb.hits.per.sequence.range[1]
  max.nb.hits <- nb.hits.per.sequence.range[2]
  
  nb.hits.per.sequence <- table(nb.hits.per.sequence)
  no.hit.nb <- total.scanned.sequences - sum(nb.hits.per.sequence)
  names(no.hit.nb) <- 0
  nb.hits.per.sequence <- append(nb.hits.per.sequence, no.hit.nb, after = 0)
  
  for(f in print.formats){
    
    if(f == "pdf"){
      TFBSs.per.seq.file <- paste(basename(prefix), "_TFBSs_per_seq/", matrix.query, "_TFBSs_per_seq.pdf", sep = "")
      pdf(TFBSs.per.seq.file)
    } else {
      TFBSs.per.seq.file <- paste(basename(prefix), "_TFBSs_per_seq/", matrix.query, "_TFBSs_per_seq.jpeg", sep = "")
      jpeg(TFBSs.per.seq.file)
    }
    
    plot(y = as.vector(nb.hits.per.sequence),
         x = as.numeric(names(nb.hits.per.sequence)),
         type = "l",
         col = "darkgreen",
         lty = 1, 
         lwd = 3,
         main = "Number of predicted TFBS per sequence",
         xlab = "Number of predicted TFBS",
         ylab = "Number of sequences"
         )
    legend("topright",
           legend= paste("Putative TFBSs: ", as.numeric(names(nb.hits.per.sequence)), " - Sequences: ", as.vector(nb.hits.per.sequence), sep = ""),
           bg="white",
           cex = 0.65
    )
    dev.off()
  }
  
  ## Get the number of putative TFBSs and the number of sequences with 
  ## at least one match of the query matrix
  nb.TFBSs <- dim(matrix.query.selection)[1]
  nb.seq <- length(as.vector(unique(matrix.query.selection$seq_id)))
  covered.sequences.per.motif[[matrix.query]] <<- as.vector(unique(matrix.query.selection$seq_id))
# })
  
  for(f in print.formats){
    
    if(f == "pdf"){
      TFBSs.pval.distribution.file <- paste(basename(prefix), "_TFBSs_pval_distribution/", matrix.query, "_TFBSs_pval_classes.pdf", sep = "")
      pdf(TFBSs.pval.distribution.file)
    } else {
      TFBSs.pval.distribution.file <- paste(basename(prefix), "_TFBSs_pval_distribution/", matrix.query, "_TFBSs_pval_classes.jpeg", sep = "")
      jpeg(TFBSs.pval.distribution.file)
    }
    
    class.counter <- 0
    
    ## Iterate in the p-val classes
    sapply(matrix.query.classes, function(pclass){
      
      ## Count the number of p-val classes per query matrix
      class.counter <<- class.counter + 1 
      
      ## Select the hits with the current pval class for the query matrix
      matrix.query.classes.selection <- matrix.query.selection[matrix.query.selection$Pval.class.letter == pclass,]
      
      ## X-Y Plot ( TFBS position vs -log10(pval) )
      if(class.counter == 1){
        plot(x = matrix.query.classes.selection$bspos,
             y = matrix.query.classes.selection$Pval.minlog10,
             # ylim = c( min(matrix.query.selection$Pval.minlog10, na.rm = TRUE), max(matrix.query.selection$Pval.minlog10, na.rm = TRUE)+0.5),
             ylim = c(min.pval.minus.log10, round(max.pval.minus.log10)),
             xlim = c(-limits, limits),
             main = paste("Distribution of TFBSs of ", matrix.query, sep = ""),
             ylab = "-log10(pval) TFBSs",
             xlab = "position (nt)",
             col = classes.to.colors[[pclass]],
             pch = "o",
             cex = 1.5,
             panel.first=grid(col = "grey", lty = "solid") 
        )
      } else {
        lines(x = matrix.query.classes.selection$bspos,
              y = matrix.query.classes.selection$Pval.minlog10,
              col = classes.to.colors[[pclass]],
              type = "p",
              pch = "o",
              cex = 1.5
              )
      }
    })
    
    ## Insert legend
    legend("topleft", legend = paste(c("Nb of putative TFBSs: ", "Nb of sequences: "), c(nb.TFBSs, nb.seq), sep = ""), bg="white")
    
    ## Insert logo
    logo.file <- paste(logo.folder, matrix.query, "_logo.jpeg", sep = "")
    logo <- readJPEG(logo.file)
    rasterImage(logo, 
                xleft = limits - (limits/3),
                xright = limits - 5, 
                ybottom = max.pval.minus.log10 - 1,
                ytop = max.pval.minus.log10 - 0.25
                )
    trash <- dev.off()
  }
})
rm(thr)
# dev.off()
# verbose(paste("Distribution of TFBSs at different p-values: ", TFBSs.pval.distribution.file), 1)

##########################
## Co-ocurrence heatmap ##
##########################
covered.seq.percentage <- vector()
thrash <- sapply(covered.sequences.per.motif, function(m1){
  sapply(covered.sequences.per.motif, function(m2){
    
    intersected.seq <- intersect(m1, m2)
    intersected.seq.nb <- length( intersected.seq)
    intersected.seq.per <- intersected.seq.nb/total.scanned.sequences
    covered.seq.percentage <<- append(covered.seq.percentage, intersected.seq.per)
  })
})
covered.seq.percentage <- round(covered.seq.percentage, digits = 4) * 100
covered.seq.percentage <- matrix(covered.seq.percentage, ncol = length(covered.sequences.per.motif))
colnames(covered.seq.percentage) <- matrix.names
rownames(covered.seq.percentage) <- matrix.names

## Set the colors
coocurrence.palette <- colorRampPalette(brewer.pal(9, "YlGnBu"), space="Lab")(20)

## Draw the co-ocurrence heatmap
verbose(paste("Creating co-ocurrence heatmap"), 1)
out.format <- c("pdf", "jpg")
heatmap.profiles <- NULL
comp.order.list.rows <- list()
comp.order.list.cols <- list()
domain <- seq(from = 0, to = 100, by = 5)

## Set the heatmap file name
co.ocurrence.heatmap.file <- paste(basename, "_co_ocurrence_heatmap.pdf", sep = "")
  
## Create the image
pdf(co.ocurrence.heatmap.file)

## Create the heatmap using 4 agglomeration rules
th <- sapply(c("average", "complete", "single", "ward"), function(m){
    
  if(m == "ward"){
    temp <- m
    m <- "ward.D"
  }
    
  ## Compute the heatmap
  hm.coocurrences <- heatmap.2(covered.seq.percentage,
                
    ## Dendrogram control
    dendrogram = "row",
              
    main = paste("co-ocurrence in sequences\n", length(covered.sequences.per.motif), " motifs - ", total.scanned.sequences, " sequences", sep = ""),
    xlab = "",
    ylab = "",
                
    hclustfun = function(d){hclust(d, method="complete")},
    distfun = function(x) Dist(x,method = 'pearson'),
                
    ## Color
    col = coocurrence.palette,
    breaks = domain,
                
    ## Trace
    trace = "none",
                
    ## Key control
    key = TRUE,
    keysize = 1.5,
    density.info = "none",
    key.xlab = "% sequences with the two motifs",
    key.ylab = "",
    key.title = "",
    # cexRow = 0.25
    offsetCol = 0.25
  )
  thrash <- dev.off()      
    
  if(m == "ward.D"){
     m <- "ward"
  }
    
  comp.order.list.rows[[m]] <<- paste(rev(hm.coocurrences[[1]]), collapse = ",")
  comp.order.list.cols[[m]] <<- paste(rev(hm.coocurrences[[2]]), collapse = ",")
})
comp.order.list.rows[["default"]] <- paste(1:length(matrix.names), collapse = ",")
comp.order.list.cols[["default"]] <- paste(1:length(matrix.names), collapse = ",")
  
## Convert the coverage table to the format required in D3 heatmap
coocurrence.table.d3 <- paste(basename(prefix), "_coocurrence_table.tsv", sep = "")
y <- NULL
for(j in 1:dim(covered.seq.percentage)[1]){
  for(i in 1:dim(covered.seq.percentage)[2]){
    y <<- rbind(y, matrix(c(j,i, as.numeric(covered.seq.percentage[j,i])), nrow = 1))
  }
}
colnames(y) <- c("Row", "Col", "Value")
y <- as.data.frame(y)
verbose(paste("Exporting data with co-ocurreence percentage for D3 dynamic heatmap", coverage.table.d3), 2)
write.table(y, file = coocurrence.table.d3, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

left <- (max(as.vector(sapply(matrix.names, nchar))) + 2) * 10

#########################################################################
## Create count table from matrix-scan results (if not exist in input) ##
#########################################################################

input.count.table <- 0
seq.count.per.motif <- list()

if(input.count.table == 0){
  
  verbose(paste("Creating counts and frequencies tables"), 1)
  counts.per.bin <-  sapply(1:length(matrix.names), function(m){
    
    ## Select the matches of the query motif
    matrix.query <- matrix.names[m]
    
    matrix.query.selection <- matrix.scan.results[matrix.scan.results$ft_name == matrix.query,]
    
    nb.seq <- length(unique(as.vector(matrix.scan.results[matrix.scan.results$ft_name == matrix.query,]$seq_id)))
    seq.count.per.motif[[matrix.query]] <<- nb.seq
    
    ## As the reference point in matrix-scan was the end of the sequence and as we are working with peaks
    ## we add 300 to the position to have -/+ position around the summit    
    if(off.set != 0){
      
      if (off.set.type == "start"){
        matrix.query.selection$bspos <- matrix.query.selection$bspos + seq.lengt - off.set
      } else if (off.set.type == "end"){  
        matrix.query.selection$bspos <- matrix.query.selection$bspos + off.set
      }
      
    } else {
      matrix.query.selection$bspos <- matrix.query.selection$bspos# + limit.dw
    }
    
    ## Convert the BSs in Ranges
    selection.IR <- IRanges(start = matrix.query.selection$bspos, end = matrix.query.selection$bspos)
    
    ## Count the overlap of BS in the bins
    counts.per.bin <- countOverlaps(windows, selection.IR)
    
    return(counts.per.bin)
  })
  
}
counts.per.bin.table <- t(counts.per.bin)
#counts.per.bin.table <- round(counts.per.bin.table, digits = 3)
rm(counts.per.bin)
rownames(counts.per.bin.table) <- matrix.names

xlab <- data.frame(windows)$start
xlab <- ifelse(xlab >= 0, xlab + bin, xlab)
colnames(counts.per.bin.table) <- as.character(xlab)

###################################
## Calculate the Frecuency table
frequency.per.bin.table <- apply(counts.per.bin.table, 1, function(r){
  
  r.sum <- sum(r)
  r/r.sum
})
frequency.per.bin.table <- t(frequency.per.bin.table)
frequency.per.bin.table  <- round(frequency.per.bin.table , digits = 3)
max.y <- max(frequency.per.bin.table, na.rm = TRUE)

##########################################
## Export Counts and Frequencies tables
density.tab.file <- paste(basename, "_density_per_bin_profiles.tab", sep = "") 
write.table(frequency.per.bin.table, file = density.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

counts.tab.file <- paste(basename, "_counts_per_bin_profiles.tab", sep = "") 
write.table(counts.per.bin.table, file = counts.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

###############################################################
## Chi-squared calculation section                           ##
## If an input table with precalculated counts is not given, ##
## it takes the table calculated in the previous section     ##
###############################################################

# if(input.count.table == 0){
#   counts.per.bin.table <- read.csv("")
# }

#####################################
## Initialize (pre-allocated) list
feature.attributes <- vector("list", dim(counts.per.bin.table)[1])
feature.log2.ratio <- vector("list", dim(counts.per.bin.table)[1])

#####################################################
## Calculate X2, p-value, e-value and significance
verbose(paste("Calculating P-value. Null Hypothesis: homogenous distribution of the TFBSs"), 1)
thrash <- sapply(1:dim(counts.per.bin.table)[1], function(m){
  
  #   print(m)
  counts.per.bin <- counts.per.bin.table[m,]
  
  # plot(x = 1:12, y = counts.per.bin.table[1,], type = "l", ylim = c(0,100))
  # plot(x = 1:12, y = counts.per.bin.table[2,], type = "l", ylim = c(0,100))
  
  # case1 <- round(counts.per.bin.table[2,]/sum(counts.per.bin.table[2,]), digits = 2)
  # mean(abs(counts.per.bin.table[2,] - min(counts.per.bin.table[2,])) )
  
  ## Select the matches of the query feature
  feature.query <- rownames(counts.per.bin.table)[m]
  feature.attributes[[m]][["feature_id"]] <<- feature.query
  
  ## Goodnes-of-fit X2 test
  ## We don't require a probability vector because we assume the TFBSs (counts) 
  ## are distributed homogenously along the sequences
  chi <- chisq.test(counts.per.bin, correct = TRUE)
  
  ## The expected values are calculated in the next way:
  ## (2 * P-val) * (Sequence_length - Motif_length + 1 )
  # motif.name <- rownames(counts.per.bin.table)[m]
  # nb.seq <- seq.count.per.motif[[motif.name]]
  # expected <- (p.val * 2 * seq.length * nb.seq)
  # nb.bins <- dim(counts.per.bin.table)[2]
  # expected <- round(expected/nb.bins)
  # expected <- rep(expected, times= nb.bins)
  # feature.log2.ratio[[m]][["feature_id"]] <<- as.vector(log2(chi[[6]]/expected))
  
  ## The expected values are calculated in the next way:
  ## (sum(nb.sites) /  Nb.seq/Nb.bin)
  # motif.name <- rownames(counts.per.bin.table)[m]
  # nb.hits <- sum(counts.per.bin)
  # nb.seq <- seq.count.per.motif[[motif.name]]
  # nb.bins <- dim(counts.per.bin.table)[2]
  # # expected <- round(nb.hits / (nb.seq/nb.bins) )
  # 
  # expected <- round((nb.hits/nb.seq)*nb.bins)
  # 
  # expected <- rep(expected, times= nb.bins)
  # feature.log2.ratio[[m]][["feature_id"]] <<- as.vector(log2(chi[[6]]/expected))
  
  ## The expected values are calculated from the Observed values
  feature.log2.ratio[[m]][["feature_id"]] <<- as.vector(log2(chi[[6]]/(chi[[7]])))
  
  # nb.bins <- dim(counts.per.bin.table)[2]
  # tfbd.med <- median(chi[[6]]) + 1
  # expected <- rep(tfbd.med, times = nb.bins)
  
  # feature.log2.ratio[[m]][["feature_id"]] <<- as.vector(log2(chi[[6]]/(chi[[7]])))
  
  # feature.log2.ratio[[m]][["feature_id"]] <<- as.vector(round(log2(counts.per.bin/median(counts.per.bin)), digits = 2))
  
  ## Chi-squared
  cs.val <- round(chi[[1]], digits = 3)
  feature.attributes[[m]][["chi_squared"]] <<- cs.val
  
  ## Degrees of freedom
  df <- chi[[2]]
  feature.attributes[[m]][["degrees"]] <<- df
  
  ## Get p-value
  chi.pval <- round(as.numeric(chi[[3]]), digits = 100000000000000000)
  
  ## Calculate e-value
  chi.eval <- chi.pval * length(matrix.names)
  
  ## Calculate significance
  sig <- round(-log10(chi.eval), digits = 3)
  feature.attributes[[m]][["significance"]] <<- sig
  
  chi.pval <- prettyNum(chi.pval, scientific=TRUE, digits = 2)
  chi.eval <- prettyNum(chi.eval, scientific=TRUE, digits = 2)
  
  feature.attributes[[m]][["pval"]] <<- chi.pval
  feature.attributes[[m]][["eval"]] <<- chi.eval
  
})
rm(thrash)
names(feature.attributes) <- matrix.names
names(feature.log2.ratio) <- matrix.names

## Convert the list in a data frame
feature.attributes <- data.frame(t(
  matrix(as.vector(unlist(feature.attributes)), 
         ncol = length(feature.attributes))))
colnames(feature.attributes) <- c("Feature", "Chi_squared", "Degrees", "Sig", "P_val", "E_val")

feature.log2.ratio <- data.frame(t(
  matrix(as.vector(unlist(feature.log2.ratio)), 
         ncol = length(feature.log2.ratio))))
rownames(feature.log2.ratio) <- matrix.names
colnames(feature.log2.ratio) <- as.character(data.frame(windows)$start)

####################################
## Calculate the profile clusters ##
####################################
verbose(paste("Calculating the motif profile clusters"),1)

#############################
## Define profile clusters
tree.profiles <- hclust(Dist(feature.log2.ratio, method = "correlation"), method = "ward.D2")
# ordered.names.tree.profiles <- tree.profiles[[4]][tree.profiles[[3]]]
clusters.tree.profiles <- cutreeDynamic(tree.profiles, minClusterSize = 1, method = "tree")
names(clusters.tree.profiles) <- tree.profiles[[4]]

## Generate a color palette
nb.profile.clusters <- length(unique(clusters.tree.profiles))
cluster.tree.profiles.palette <- colorRampPalette(brewer.pal(9, "Set1"), space="Lab")(nb.profile.clusters)

## Assign a different color to each cluster
color.clusters.tree.profiles <- as.vector(sapply(clusters.tree.profiles, function(color){
    cluster.tree.profiles.palette[color]
}))

###############################################
## Fill the Profile_cluster attribute
profile.clusters.names <- as.vector(sapply(matrix.names, function(m){
  profile.cluster <- as.vector(clusters.tree.profiles[m])
  paste("Profile_cluster_", profile.cluster, sep = "")
}))
profile.clusters.names.unique <- unique(profile.clusters.names)
feature.attributes$Profile_cluster <- profile.clusters.names

## Remove the special characters -there characters broke the CSS variables-
## Only for the CSS section
motifs.names.parsed <- as.vector(sapply(matrix.names, function(m){
  m.temp <- gsub("-", "", m)
  m.temp <- gsub("\\.", "_", m.temp)
  m.temp <- gsub(":", "", m.temp)
  m.temp <- gsub("\\s+", "", m.temp, perl = TRUE)
  return(m.temp)
}))

## Get the member motif IDs of each cluster
cluster.profiles.motifs <- list()
cluster.profiles.motif.names <- list()
cluster.profiles.counter <- 0
thrash <- sapply(1:nb.profile.clusters, function(cl){
  
  cluster.profiles.counter <<- cluster.profiles.counter + 1
  cluster.profiles.motifs[[cluster.profiles.counter]] <<- names(which(clusters.tree.profiles == cluster.profiles.counter))
  
  ## Get the motif name
  cluster.profiles.motif.names[[cluster.profiles.counter]] <<- as.vector(
    sapply(cluster.profiles.motifs[[cluster.profiles.counter]], function(n){
      ID.names[which(ID.names[,2] == n),1]
    })
  )
})
rm(thrash)

####################################################################################
## Draw Profiles heatmap showing the frequencies of hits per bin for each feature ##
####################################################################################
verbose(paste("Drawing Heatmap profiles"),1)
## Profile Heatmap Color palette (user-defined)
rgb.palette <- rev(colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")(heatmap.color.classes))

log2.tab <- as.matrix(feature.log2.ratio)
log2.tab[is.infinite(log2.tab)] <- 0

## Print the heatmap
out.format <- c("pdf", "jpg")
heatmap.profiles <- NULL
for(format in out.format){
  
  profiles.heatmap.file <- paste(basename, "_profiles_heatmap.", format, sep = "") 
  
  if(format == "pdf"){
    pdf(profiles.heatmap.file)
  } else if (format == "jpg"){
    jpeg(profiles.heatmap.file)
  }
  
  ## Heatmap
  heatmap.profiles <<- heatmap.2(log2.tab,
                   
                   main = "Profile Heatmap",
                   xlab = "Position (bp)",
                   ylab = "Motifs",
                   
                   ## The order of the values is set according these dendrograms
                   Rowv = as.dendrogram(tree.profiles),
                   Colv = FALSE,
                   dendrogram = "row",
                   
                   ## Color
                   col = rgb.palette,
                   
                   ## Trace
                   trace = "none",
                   
                   ## Side colors
                   RowSideColors = color.clusters.tree.profiles,
                   
                   ## Key control
                   key = TRUE,
                   keysize = 1,
                   density.info = "none",
                   key.xlab = "Log2 Ratio",
                   key.ylab = "",
                   key.title = "",
                   # cexRow = 0.25
                   offsetCol = 0.25
  )
  trash <- dev.off()
}
heatmap.row.order <- rev(heatmap.profiles[[1]])
heatmap.row.order.names <- rownames(feature.log2.ratio)[heatmap.row.order]
heatmap.rows <- data.frame(row = heatmap.row.order, names = heatmap.row.order.names)
heatmap.rows.file <- paste(basename, "_heatmap_row_order.tab", sep = "")
write.table(heatmap.rows, file = heatmap.rows.file, quote = FALSE, sep = "\t", row.names = FALSE)

## Calculate q-values
## This step is executed once all the p-values were calculated
## The variable with class 'qvalue' is stored to its further exportation
pp <- as.numeric(as.vector(feature.attributes$P_val))
features.qvalues <- p.adjust(pp, method = "BH")
feature.attributes$Q_val <- prettyNum(features.qvalues, scientific=TRUE, digits = 2)

############################################################
## Additional columns                                     ##
## For matrix-scan results, but maybe not for other input ##
## (i.e. oligo-analysis table)                            ##
############################################################

## Add number of hits column
feature.attributes[,"Nb_hits"] <- as.vector(apply(counts.per.bin.table, 1, sum))

## Get additional data
## 1) Nb of sequences with one hit
## 2) Max p-val (among hits)
## 3) Min p-val (among hits)
additional.data <- vector("list", dim(counts.per.bin.table)[1]) 
thrash <- sapply(1:dim(counts.per.bin.table)[1], function(f){
  
  feature.query <- rownames(counts.per.bin.table)[f]
  
  ## Select the matches of the query motif
  matrix.query.selection <- matrix.scan.results[matrix.scan.results$ft_name == feature.query,]
  
  ## Get the number of sequences with at least one hit
  nb.seq.with.hits <- length(unique(as.vector(matrix.query.selection$seq_id)))
 
  ## Export the covered/non_covered sequences names 
  covered.seq <- unique(as.vector(matrix.query.selection$seq_id))
  not.covered.seq <- setdiff(scanned.sequences, covered.seq)

  covered.sequences.table <- as.data.frame(covered.seq)
  not.covered.sequences.table <- as.data.frame(not.covered.seq)
  
  overed.sequences.file <- file.path(covered.tables.dir, paste(feature.query, "_covered_sequences_IDs.tab", sep = ""))
  not.covered.sequences.file <- file.path(covered.tables.dir, paste(feature.query, "_not_covered_sequences_IDs.tab", sep = ""))

  # print("Aqui1")
  # # write.table(covered.sequences.table, file = covered.sequences.file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  # # write.table(not.covered.sequences.table, file = not.covered.sequences.file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  # print("Aqui2")
  
  ## Calculate the coverture rate
  coverture <- round(nb.seq.with.hits/total.scanned.sequences, digits = 4)*100
  additional.data[[f]][["Coverture"]] <<- paste(coverture, "%", sep = "")
  
  additional.data[[f]][["sequences"]] <<- nb.seq.with.hits
  
  ## Get the max p-value among the p-values of the matches
  max.pval.hit <- max(matrix.query.selection$Pval)
  additional.data[[f]][["max_pval"]] <<- max.pval.hit
  
  ## Get the min p-value among the p-values of the matches
  min.pval.hit <- min(matrix.query.selection$Pval)
  additional.data[[f]][["min_pval"]] <<- min.pval.hit
})
rm(thrash)
names(additional.data) <- matrix.names

## Convert the list in a data frame
additional.data <- data.frame(t(
  matrix(as.vector(unlist(additional.data)), 
         ncol = length(additional.data))))
colnames(additional.data) <- c("Coverture", "Nb_sequences", "Max_pval", "Min_pval")


#################################################################
## Merge the dataframes (additional data + feature attributes) ##
#################################################################
feature.attributes <- cbind(feature.attributes, additional.data)
feature.attributes  <- feature.attributes[,c(1,7,5,6,4,8,2,3,9,11,10)]

feature.attributes.file <- paste(basename, "_attributes.tab", sep = "")
# write.table(feature.attributes, file = feature.attributes.file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# rm(additional.data)

###############################################################
## Compute the XY-plot for Profile significance vs Coverture ##
###############################################################

verbose(paste("Drawing Significance vs Coverture plot"),1)

## Calculate X-Y coordinates
x.sig <- as.numeric(as.vector(feature.attributes$Sig))
x.sig[x.sig == Inf] <- 350
names(x.sig) <- as.vector(feature.attributes$Feature)
y.cov <- as.numeric(gsub("%", "", feature.attributes$Coverture))
names(y.cov) <- as.vector(feature.attributes$Feature)

## Print the plot
sig.coverture.file <- paste(basename, "_significance_vs_coverture.pdf", sep = "") 
pdf(sig.coverture.file)

## X-Y plot
plot(x.sig,
     y.cov,
     ylim = c(0,100),
     xlab = "Significance -log10(Corrected p-val)",
     ylab = "Sequence coverture",
     main = "Profile Significance vs Sequence Coverture",
     col = ifelse((x.sig >= 20 & y.cov >= 66), "darkgreen", "gray"),
     panel.first=grid(col = "grey", lty = "solid"),
     pch = "o",
     cex = 1.5
     )

## Mark the TFBMs sattisfying the threshold
selected.TFBMs <- feature.attributes[which(as.vector(feature.attributes$Sig) >= 20 & as.numeric(gsub("%", "", feature.attributes$Coverture)) >= 66), "Feature"]
selected.TFBMs <- as.vector(selected.TFBMs)

if(length(selected.TFBMs) > 0){
  
  ## Add the text to the selected TFBMs
  text(x = x.sig[c(selected.TFBMs)],
       y = y.cov[c(selected.TFBMs)],
       labels = selected.TFBMs,
       cex = 0.6, 
       pos = 3, 
       col="red")
}
thrash <- dev.off()

## Convert the X-Y values to the format required for C3 plot
xx.sig <- paste("['x',", paste(round(as.vector(x.sig)), collapse = ","), "],", sep = "")
yy.cov <- paste("['y',", paste(round(as.vector(y.cov)), collapse = ","), "],", sep = "")
x.y.coverture <- paste(xx.sig, yy.cov, collapse = "\n")

##############################################################
## Plot each profile individually (if it is user-specified) ##
## Require the icon of the feature (e.g. PSSM logo)         ##
##############################################################
# individual.plots <- 0
if(individual.plots == 1){
  
  ## Create folder for individual profile plots
  dir.create(paste(basename(prefix), "_TFBSs_positional_profiles/", sep = ""), recursive = TRUE, showWarnings = FALSE )
  
  verbose(paste("Printing all the profiles in a PDF file"), 1)
  
  # pdf.file.name <- paste(basename, "_positional_profiles.pdf", sep = "")
  # pdf(pdf.file.name)
  
  thrash <- sapply(1:dim(frequency.per.bin.table)[1], function(f){
    
    feature.query <- rownames(frequency.per.bin.table)[f]
    # matrix.query.ID <- ID.names[,1][which(ID.names[,2] == feature.query)][1]
    
    for(pf in print.formats){
      
      if(pf == "pdf"){
        pdf.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", feature.query, "_positional_profile.pdf", sep = "")
        pdf(pdf.file.name)
      } else {
        jpeg.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", feature.query, "_positional_profile.jpeg", sep = "")
        jpeg(jpeg.file.name)
      }
      
      y.val <- frequency.per.bin.table[f,]
      x.val <- as.numeric(colnames(frequency.per.bin.table))
      
      ## Draw the lines for the active promoters 
      # lines(x = x.val, y = y.val, )
      
      plot(x = x.val,
           y = y.val,
           type = "l",
           col = "#00BFC4",
           lty = 1, 
           lwd = 3,
           ylim = c(0, max(max.y)),

           ## Labels
           main = paste("Motif:", feature.query),
           xlab = "Distance to center",
           ylab = "TFBSs fraction",
           ## Hide x-axis
           xaxt='n'
      ) 
      
      ## Draw the grid
      abline(v=(x.val), col="lightgray", lty="dotted")
      abline(h=(seq(from = 0, to = 1, by = 0.01)), col="lightgray", lty="dotted")
      
      ## Draw the TSS (+1) position
      abline(v = 0, col="#045a8d", lwd = 2, lty = 2)
      
      ## Set x-axis values 
      axis(side = c(1,2,3,4), at = as.character(x.val), labels = as.character(x.val))
      
      # ## Draw the lines for the active promoters 
      # lines(x = x.val, y = y.val, type = "l", col = "#00BFC4", lty = 1, lwd = 3)
      
      ## Draw the legend
      legend("topleft", legend = c(paste(feature.query , "profile"), "Center"), fill = c("#00BFC4", "#045a8d"), bty="o", bg="white")
      
      logo.file <- paste(logo.folder, feature.query, "_logo.jpeg", sep = "")
      logo <- readJPEG(logo.file)
      rasterImage(logo, 
                  xleft = limits - (bin*3),
                  xright = limits, 
                  ybottom = max.y - 0.1,
                  ytop = max.y - 0.05)
      trash <- dev.off()
    }
  })
  # dev.off()    
}


#################################################################
## Create the dynamic report in HTML                           ##
## This is a general procedure (not only matrix-scan specific) ##
## A HTML template is modified with the current data stored in ##
## the dataframe frequency.per.bin.table complemented with the ##
## statistics calculated in the dataframe feature.attributes   ##
#################################################################

verbose(paste("Creating HTML dynamic report"), 1)

## Order the TF.names (and other variables) according the Significance (-log10(E-value))
order.by.eval <- order(as.numeric(as.vector(feature.attributes$E_val)))
feature.attributes <- feature.attributes[order.by.eval,]

## Set the motif names and IDs
TF.IDs <- as.vector(feature.attributes[,1])

## Get the IDs of the TF
TF.names <- as.vector(sapply(TF.IDs, function(x){
  as.vector(ID.names[which(ID.names[,1] == x),2])
}))

## Set colors
set.colors <- colorRampPalette(brewer.pal(10,"Paired"))(length(TF.IDs))

counter <- 0
x.correspondence <- NULL
x.y <- NULL
x.y.coverture <- NULL
x.y.coverture.names <- NULL
plot.names <- NULL
plot.names.cover <- NULL
area <- NULL
all.motifs <- NULL
all.motifs.cover <- NULL
all.motif.names <- NULL
hash.motif.IDs <- list()


## Each column of the variable profiles correspond to the counts per bin of each motif
thrash <- apply(frequency.per.bin.table[order.by.eval,], 1, function(values){
  
  counter <<- counter + 1
  
  ## Here we create a unique ID without CSS special characters
  ## Only to manipulate the objects in the HTML form
  motif <- paste(counter, "_", TF.IDs[counter], "_", counter, sep = "")  
  motif.cover <- paste(counter, counter, "_", TF.IDs[counter], "_", counter, counter, sep = "")
  
  motif <- gsub("_", "", motif)
  motif <- gsub("-", "", motif)
  motif <- gsub("\\.", "", motif)
  motif <- gsub(":", "", motif)
  motif <- gsub("\\s+", "", motif, perl = TRUE)
  
  hash.motif.IDs[[TF.IDs[counter]]] <<- motif
  
  all.motifs <<- append(all.motifs, motif)
  all.motif.names <<- append(all.motif.names, TF.names[counter])
  
  ## Create the name's data for the HTML file
  ## The name correspond to the motif name. 
  ## Note that two motifs for the same TF will have the same name and this name will appears twice 
  ## in the report. However their IDs are unique.
  plot.names <<- append(plot.names, paste("'", motif, "' : '",  TF.names[counter],"',", sep = ""))
  
  ## Create the area's data for the HTML file (optional)
  area <<- append(area, paste("'", motif, "' : 'area',", sep = ""))
  
  ## Create the X data
  ## For each TF we will have a row containing the ID plus the counts per bin
  ## The X-axis data for the plot
  y <- paste("['", 
             motif, 
             "',",
             paste(values, 
                   collapse = ","),
             "],",
             sep = "")
  x.y <<- rbind(x.y, y) 
  
  ## Convert the X-Y values to the format required for C3 plot coverture
  xx.sig <- paste("['", motif.cover, "_x',", round(x.sig[order.by.eval])[counter], "],", sep = "")
  yy.cov <- paste("['", motif.cover, "',", round(y.cov[order.by.eval])[counter], "],", sep = "")
  x.y.coverture <<- rbind(x.y.coverture, xx.sig)
  x.y.coverture <<- rbind(x.y.coverture, yy.cov)
  
  ## Add the motifs IDs sentences for the coverture plot
  plot.names.cover <<- append(plot.names.cover, paste("'", motif.cover, "' : '",  TF.names[counter],"',", sep = ""))
  
  ## Append all the motifs IDs for the coverture plot
  all.motifs.cover <<- append(all.motifs.cover, motif.cover)
    
  ## Add the motif names for the coverture plot
  name.cov <- paste("'", motif.cover, "' : '", motif.cover, "_x',", sep = "")
  x.y.coverture.names <<- rbind(x.y.coverture.names, name.cov)

})
all.motifs.cover.hash <- data.frame(all.motifs.cover, all.motif.names)


## Set the line width according the significance -log10(E-value)
## Higher significance means a wider line
significance <- as.numeric(as.vector(feature.attributes$Sig))
line.w <- sapply(significance, function(s){
  if(s <= 0){
    w <- 1
  } else if (s > 0 & s <= 10){
    w <- 2
  } else if (s > 10 & s <= 50){
    w <- 3
  }  else if (s > 50 & s <= 100){
    w <- 4
  } else if (s > 100){
    w <- 5
  }
})


##############################################
## Insert the cluster functions and buttons
all.cluster.functions <- vector()
all.cluster.buttons <- vector()
all.cluster.functions.cov <- vector()
all.cluster.buttons.cov <- vector()

## JS code to show the motifs corresponding to one cluster
## We require one function for each cluster
show.profile.cluster.function <- 'function --profile_cluster_show--() {
chart.hide([--all_motifs_function--]);
chart.show([--names_function--]);
}'
show.profile.cluster.button <- "<button class='small button_chart' onclick='--function_name--();'>--cluster_name--</button>"

## JS code for coverage plot
show.profile.cluster.function.cover <- 'function --profile_cluster_show--() {
coverchart.hide([--all_motifs_function--]);
coverchart.show([--names_function--]);
}'
show.profile.cluster.button.cover <- "<button class='small button_chart' onclick='--function_name--();'>--cluster_name--</button>"

all.motifs.function <- paste(paste("'", all.motifs, "'", sep = ""), collapse = ",")
all.motifs.function.cov <- paste(paste("'", all.motifs.cover, "'", sep = ""), collapse = ",")

## Generate the JS function to show the clusters
thrash <- sapply(1:length(cluster.profiles.motif.names), function(cl){
  
  ## Define the profile cluster name
  profile.cluster.name <- paste("Profile_cluster_", cl, sep = "")
  
  cluster.function <- show.profile.cluster.function
  cluster.function.cov <- show.profile.cluster.function.cover
  
  ##############################################################
  ## Change the function's name for the corresponding cluster
  
  ## Button name
  cluster.function.names <- paste(profile.cluster.name, "_show", sep = "")
  cluster.function.names.cov <- paste(profile.cluster.name, "_show_cov", sep = "")
  
  ## Cluster member names
  cluster.member.names <- as.vector(unlist(sapply(as.vector(unlist(cluster.profiles.motif.names[[cl]])), function(m){
    as.vector(hash.motif.IDs[[m]])
  })))
  cluster.member.names <- paste(paste("'", cluster.member.names, "'", sep = ""), collapse = ",")
  
  ## Cluster member names (cover plot)
  cluster.member.names.cov <- as.vector(unlist(sapply(as.vector(unlist(cluster.profiles.motif.names[[cl]])), function(m){
    as.vector(all.motifs.cover.hash[which(all.motifs.cover.hash[,2] == m),1])
  })))
  cluster.member.names.cov <- paste(paste("'", cluster.member.names.cov, "'", sep = ""), collapse = ",")
  
  ## Substitution
  cluster.function <- gsub("--profile_cluster_show--", cluster.function.names, cluster.function)
  cluster.function <- gsub("--all_motifs_function--", all.motifs.function, cluster.function)
  cluster.function <- gsub("--names_function--", cluster.member.names, cluster.function)
  
  ## Substitution (cov)
  cluster.function.cov <- gsub("--profile_cluster_show--", cluster.function.names.cov, cluster.function.cov)
  cluster.function.cov <- gsub("--all_motifs_function--", all.motifs.function.cov, cluster.function.cov)
  cluster.function.cov <- gsub("--names_function--", cluster.member.names.cov, cluster.function.cov)
  
  ## Concat all the functions
  all.cluster.functions <<- append(all.cluster.functions, cluster.function)
  all.cluster.functions.cov <<- append(all.cluster.functions.cov, cluster.function.cov)
  
  ##########################################################
  ## Create the button to show all the motifs in a cluster
  cluster.button <- show.profile.cluster.button
  cluster.button <- gsub("--cluster_name--", profile.cluster.name, cluster.button)
  cluster.button <- gsub("--function_name--", cluster.function.names, cluster.button)
  
  ## Create the button to show all the motifs in a cluster (cov)
  cluster.button.cov <- show.profile.cluster.button
  cluster.button.cov <- gsub("--cluster_name--", profile.cluster.name, cluster.button.cov)
  cluster.button.cov <- gsub("--function_name--", cluster.function.names.cov, cluster.button.cov)
  
  ## Concat all the functions
  all.cluster.buttons <<- append(all.cluster.buttons, cluster.button)
  all.cluster.buttons.cov <<- append(all.cluster.buttons.cov, cluster.button.cov)

})
all.cluster.functions <- paste(all.cluster.functions, collapse = "\n")
all.cluster.buttons <- paste(all.cluster.buttons, collapse = "\n")
all.cluster.functions.cov <- paste(all.cluster.functions.cov, collapse = "\n")
all.cluster.buttons.cov <- paste(all.cluster.buttons.cov, collapse = "\n")

## Write the logo's path
logos.F <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo.jpeg", sep = "")
})

logos.R <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo_rc.jpeg", sep = "")
})

## Write the Profile and TFBSs plots path
profiles.plots <- sapply(TF.IDs, function(i) {
  paste(basename(prefix), "_TFBSs_positional_profiles/", i, "_positional_profile.jpeg", sep = "")
})

tfbss.plots <- sapply(TF.IDs, function(i) {
  paste(basename(prefix), "_TFBSs_pval_distribution/", i, "_TFBSs_pval_classes.jpeg", sep = "")
})

tfbss.per.seq.plots <- sapply(TF.IDs, function(i) {
  paste(basename(prefix), "_TFBSs_per_seq/", i, "_TFBSs_per_seq.jpeg", sep = "")
})

## Write the path to the covered/non_covered sequences tables
covered.files <- sapply(TF.IDs, function(i) {
  paste(covered.tables.dir, i, "_covered_sequences_IDs.tab", sep = "")
})
not.covered.files <- sapply(TF.IDs, function(i) {
  paste(covered.tables.dir, i, "_not_covered_sequences_IDs.tab", sep = "")
})

## Create a Dataframe containing the information of all motifs
## This table will be exported and displayed as a dynamic table in the report
all.pval.match <- rep(p.val, times = length(TF.names))
datatable.info.tab <- feature.attributes
datatable.info.tab$P_val_threshold <- all.pval.match
datatable.info.tab$Names <- TF.names
datatable.info.tab$Profiles <- profiles.plots
datatable.info.tab$TFBS <- tfbss.plots
datatable.info.tab$TFBS_per_seq <- tfbss.per.seq.plots
datatable.info.tab$Logo <- logos.F
datatable.info.tab$Logo_RC <- logos.R
datatable.info.tab$covered_files <- covered.files
datatable.info.tab$not_covered_files <- not.covered.files
all.motifs <- all.motifs

############################
## Fill the HTML template
## Substitute the words marked in the template by the data
html.report <- readLines(html.template.file)
# [1] "Feature"         "Profile_cluster"           "P_val"           "E_val"          
# [5] "Sig"             "Q_val"           "Chi_squared"     "Degrees"        
# [9] "Nb_hits"         "Nb_sequences"    "Coverture"       "P_val_threshold"
# [13] "IDs"             "Profiles"        "TFBS"            "TFBS_per_site"            
# [17] "Logo"           "Logo_RC"
profile.data.tab.html <- create.html.tab(datatable.info.tab[,c(1,13,3:6,9:11,7,2,14:20)], img = c(15,16), plot = c(12,13,14), link.text.covered = 17, link.text.not.covered = 18)

profile.data.tab.html <- gsub("Inf", "&infin;", profile.data.tab.html)

profile.data.tab.html <- paste(profile.data.tab.html, collapse = "\n")
html.report <- gsub("--tab--", profile.data.tab.html, html.report)

## Define the x-axis categories
x.axis.categories <- paste(paste("'", colnames(frequency.per.bin.table), "'", sep = ""), collapse = ",")

## CSS section to set the line width
## Note: the width is proportional to the significance
line.w <- paste("#chart .c3-line-", all.motifs, "{ stroke-width: ", line.w, "px; }", sep = "")
line.w <- paste(line.w, collapse = "\n")
html.report <- gsub("--lines_w--", line.w, html.report)

## Add the TF_names data
TF.names <- paste("TF_names['", all.motif.names, "'] = '", all.motifs, "';", sep = "")
TF.names <- paste(TF.names, collapse = "\n")
html.report <- gsub("--TF_names--", TF.names, html.report)

## Add the TF_names data
tfs <- paste(paste("'", all.motif.names, "'", sep = ""), collapse = ",")
html.report <- gsub("--tfs--", tfs, html.report)

## Add the e-values data
## They are inserted in the JS section
evalues <- paste("evalues['", all.motifs, "'] = '", as.vector(datatable.info.tab$E_val), "';", sep = "")
evalues <- paste(evalues, collapse = "\n")
html.report <- gsub("--evalues--", evalues, html.report)

## Add the p-values (to display in the tooltip)
## They are inserted in the JS section
pvalues <- paste("pvalues['", all.motifs, "'] = '", as.vector(datatable.info.tab$P_val), "';", sep = "")
pvalues <- paste(pvalues, collapse = "\n")
html.report <- gsub("--pvalues--", pvalues, html.report)

## Add the profile clusters (to display in the tooltip)
## They are inserted in the JS section
profile.clusters.array <- paste(" profile_clusters['", all.motifs, "'] = '", as.vector(datatable.info.tab$Profile_cluster), "';", sep = "")
profile.clusters.array <- paste(profile.clusters.array, collapse = "\n")
html.report <- gsub("--profile_clusters_array--",  profile.clusters.array, html.report)

## Add the q-values (to display in the tooltip)
## They are inserted in the JS section
qvalues <- paste("qvalues['", all.motifs, "'] = '", as.vector(datatable.info.tab$Q_val), "';", sep = "")
qvalues <- paste(qvalues, collapse = "\n")
html.report <- gsub("--qvalues--", qvalues, html.report)

## Add the real motif IDs (to display in the tooltip)
## They are inserted in the JS section
## I called them 'real' because are those found on the original motif file
IDs <- paste("IDs['", all.motifs, "'] = '", TF.IDs, "';", sep = "")
IDs <- paste(IDs, collapse = "\n")
html.report <- gsub("--IDs--", IDs, html.report)

## Add the real motif logo path (to display in the tooltip)
## They are inserted in the JS section
logos <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo.jpeg", sep = "")
})
logos <- paste("pics['", all.motifs, "'] = '", as.vector(datatable.info.tab$Logo), "';", sep = "")
logos <- paste(logos, collapse = "\n")
html.report <- gsub("--pics--", logos, html.report)

logos.co <- paste("pics_m1['", all.motifs, "'] = '", as.vector(datatable.info.tab$Logo), "';", sep = "")
logos.co <- paste(logos.co, collapse = "\n")
html.report <- gsub("--pics_m1--", logos.co, html.report)

## Logos in Reverse complement
logos.rc <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo_rc.jpeg", sep = "")
})
logos.rc <- paste("pics_rc['", all.motifs, "'] = '", as.vector(datatable.info.tab$Logo_RC), "';", sep = "")
logos.rc <- paste(logos.rc, collapse = "\n")
html.report <- gsub("--pics_rc--", logos.rc, html.report)

## Add the signficance (to display in the tooltip)
## They are inserted in the JS section
sig <- paste("significances['", all.motifs, "'] = ", as.vector(datatable.info.tab$Sig), ";", sep = "")
sig <- paste(sig, collapse = "\n")
html.report <- gsub("--significances--", sig, html.report)

## Add the covertures (to display in the tooltip)
## They are inserted in the JS section
cc <- as.numeric(gsub("%", "", feature.attributes$Coverture))
coverture <- paste("TF_coverture['", all.motifs, "'] = ", as.vector(cc), ";", sep = "")
coverture <- paste(coverture, collapse = "\n")
html.report <- gsub("--TF_covertures--", coverture, html.report)

## The plot heigth depends in the number of motifs
motif.total <- length(all.motifs)
chart.heigth <- 500
if(motif.total >= 200){
  chart.heigth <- 800
} else if(motif.total >= 300){
  chart.heigth <- 1000
} else if(motif.total >= 600){
  chart.heigth <- 1400
}
html.report <- gsub("--chart_h--", chart.heigth, html.report)

## Add x values (one row per motif)
## They are inserted in the C3 section
xx <- paste(x.y, collapse = "\n")
html.report <- gsub("--x_y--", xx, html.report)

## Add the X-axis categories
html.report <- gsub("--categories--", x.axis.categories, html.report)

## Add the color code (one color per motif)
## They are inserted in the C3 section
set.colors <- paste(paste("'", sample(set.colors), "'", sep = ""), collapse = ",")
# set.colors <- paste(paste("'", rev(set.colors), "'", sep = ""), collapse = ",")
html.report <- gsub("--color_pattern--", set.colors, html.report)

## Insert the motif names
## They are inserted in the C3 section
plot.names <- paste(plot.names, collapse = "\n")
html.report <- gsub("--names--", plot.names, html.report)

## Option to draw the AUC of each motif profile
## It is inserted in the C3 section
## Only if specified by the user
if(draw.area == 1){
  area <- paste(area, collapse = "\n")
  html.report <- gsub("--area--", area, html.report)
}

## Insert the Y axis limits
## They are inserted in the C3section
max.y <- max(frequency.per.bin.table) + 0.02
html.report <- gsub("--y_axis--", max.y, html.report)

## Fill the parameters table
html.report <- gsub("--bin_l--", bin, html.report)
html.report <- gsub("--bin_nb--", ncol(counts.per.bin.table), html.report)
html.report <- gsub("--seq_nb--", total.scanned.sequences, html.report)     ## Don't forget length(seq.id)
html.report <- gsub("--motif_nb--", length(matrix.names), html.report)
html.report <- gsub("--p--", prettyNum(p.val), html.report)

## Fill the heatmap section
html.report <- gsub("--heatmap_png--", paste(basename, "_profiles_heatmap.jpg", sep = ""), html.report)
html.report <- gsub("--heatmap_pdf--", paste(basename, "_profiles_heatmap.pdf", sep = ""), html.report)

## Insert the Hit counts per bin table
html.report <- gsub("--hit_counts_table--", counts.tab.file, html.report)

## Insert the density of hits per bin table
html.report <- gsub("--density_table--", density.tab.file, html.report)

## Insert the Full motif description table
html.report <- gsub("--full_motif_table--", feature.attributes.file, html.report)

## Insert the Full motif description table
html.report <- gsub("Inf;", "Infinity;", html.report)

## Insert the Full motif description table
min.sig <- min(as.numeric(as.vector(datatable.info.tab$Sig)))
html.report <- gsub("--sig_min--", min.sig, html.report)

## Insert the JavaScript libraries path
html.report <- gsub("--c3_css--", c3.css.base, html.report)
html.report <- gsub("--c3--", c3.base, html.report)
html.report <- gsub("--d3--", d3.base, html.report)
html.report <- gsub("--jquery--", jquery.base, html.report)
html.report <- gsub("--datatable--", datatable.base, html.report)
html.report <- gsub("--datatable_css--", datatable.css.base, html.report)

## Insert the X-Y scatterplot values
x.y.coverture <- paste(x.y.coverture, collapse = "\n")
html.report <- gsub("--x_y_coverture--", x.y.coverture, html.report)

## Insert the X-Y scatterplot xs
x.y.coverture.names <- paste(x.y.coverture.names, collapse = "\n")
html.report <- gsub("--xs_coverture--", x.y.coverture.names, html.report)

## Insert the names in the coverture XY-plot
plot.names.cover <- paste(plot.names.cover, collapse = "\n")
html.report <- gsub("--names_cov--", plot.names.cover, html.report)

## Add the real motif IDs (to display in the tooltip)
## They are inserted in the JS section
IDs.cov <- paste("cov_IDs['", all.motifs.cover, "'] = '", TF.IDs, "';", sep = "")
IDs.cov <- paste(IDs.cov, collapse = "\n")
html.report <- gsub("--IDs_cov--", IDs.cov, html.report)

## Add the profile cluster (to display in the tooltip)
## They are inserted in the JS section
cov.profile.clusters <- paste("cov_profile_clusters['", all.motifs.cover, "'] = '", as.vector(datatable.info.tab$Profile_cluster), "';", sep = "")
cov.profile.clusters <- paste(cov.profile.clusters, collapse = "\n")
html.report <- gsub("--profile_clusters_array_cov--", cov.profile.clusters, html.report)

## Add the real motif logo path (to display in the tooltip)
## They are inserted in the JS section
logos.cov <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo.jpeg", sep = "")
})
logos.cov <- paste("cov_pics['", all.motifs.cover, "'] = '", as.vector(datatable.info.tab$Logo), "';", sep = "")
logos.cov <- paste(logos.cov, collapse = "\n")
html.report <- gsub("--pics_cov--", logos.cov, html.report)

## Logos in Reverse complement
logos.rc.cov <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo_rc.jpeg", sep = "")
})
logos.rc.cov <- paste("cov_pics_rc['", all.motifs.cover, "'] = '", as.vector(datatable.info.tab$Logo_RC), "';", sep = "")
logos.rc.cov <- paste(logos.rc.cov, collapse = "\n")
html.report <- gsub("--pics_rc_cov--", logos.rc.cov, html.report)

## Add the signficance (to display in the tooltip)
## They are inserted in the JS section
ss <- as.numeric(gsub("%", "", feature.attributes$Sig))
ss[ss == Inf] <- 350
sig.cov <- paste("cov_significances['", all.motifs.cover, "'] = ", as.vector(ss), ";", sep = "")
sig.cov <- paste(sig.cov, collapse = "\n")
html.report <- gsub("--significances_cov--", sig.cov, html.report)

## Add the covertures (to display in the tooltip)
## They are inserted in the JS section
cc <- as.numeric(gsub("%", "", feature.attributes$Coverture))
coverture.cov <- paste("cov_TF_coverture['", all.motifs.cover, "'] = ", as.vector(cc), ";", sep = "")
coverture.cov <- paste(coverture.cov, collapse = "\n")
html.report <- gsub("--TF_covertures_cov--", coverture.cov, html.report)

## Add the covertures (to display in the tooltip)
## They are inserted in the JS section
all.profiles.pics <- paste("'", as.vector(datatable.info.tab$Profiles), "'", sep = "")
profiles.pics.cov <- paste("cov_pics_profile['", all.motifs.cover, "'] = ", all.profiles.pics, ";", sep = "")
profiles.pics.cov <- paste(profiles.pics.cov, collapse = "\n")
html.report <- gsub("--profile_pics_cov--", profiles.pics.cov, html.report)

all.profiles.pics.co <- paste("'", as.vector(datatable.info.tab$Profiles), "'", sep = "")
profiles.pics.cov.co <- paste("profiles_m1['", all.motifs, "'] = ", all.profiles.pics.co, ";", sep = "")
profiles.pics.cov.co <- paste(profiles.pics.cov.co, collapse = "\n")
html.report <- gsub("--profiles_m1--", profiles.pics.cov.co, html.report)

## Insert the motif names (to hide/show all) in coverture plot
## They are inserted in the JQuery section
all.motifs.cover <- paste(paste("'", all.motifs.cover, "'", sep = ""), collapse = ",")
html.report <- gsub("--all_cover--", all.motifs.cover, html.report)

###################################################
## Fill the data for the co-ocurrence d3 heatmap
html.report <- gsub("--default_r_number--", comp.order.list.rows[["default"]], html.report)
html.report <- gsub("--default_c_number--", comp.order.list.cols[["default"]], html.report)

html.report <- gsub("--average_r_number--", comp.order.list.rows[["average"]], html.report)
html.report <- gsub("--average_c_number--", comp.order.list.cols[["average"]], html.report)

html.report <- gsub("--complete_r_number--", comp.order.list.rows[["complete"]], html.report)
html.report <- gsub("--complete_c_number--", comp.order.list.cols[["complete"]], html.report)

html.report <- gsub("--single_r_number--", comp.order.list.rows[["single"]], html.report)
html.report <- gsub("--single_c_number--", comp.order.list.cols[["single"]], html.report)

html.report <- gsub("--ward_r_number--", comp.order.list.rows[["ward"]], html.report)
html.report <- gsub("--ward_c_number--", comp.order.list.cols[["ward"]], html.report)

html.report <- gsub("--c_numb--", length(matrix.names), html.report)
html.report <- gsub("--r_numb--", length(matrix.names), html.report)

cell.size <- 20
html.report <- gsub("--cell_size--", cell.size , html.report)

html.report <- gsub("--left--", left, html.report)

html.report <- gsub("--heatmap_coocurrence_pdf--", co.ocurrence.heatmap.file, html.report)

coocurrence.table.d3 <- paste(basename(prefix), "_coocurrence_table.tsv", sep = "")
html.report <- gsub("--file--", coocurrence.table.d3, html.report)

domain <- paste(domain, collapse=",")
html.report <- gsub("--domain--", domain, html.report)

coocurrence.palette <- c("#FFFFFF", coocurrence.palette) 
gradient <- paste("[", paste(paste("'", coocurrence.palette, "'", sep=""), collapse=","), "];", sep = "")
html.report <- gsub("--gradient--", gradient, html.report)

legend <- seq(from = 0, to = 100, by = 5)
legend <- paste(legend, collapse=",")
html.report <- gsub("--data_legend--", legend, html.report)

col.labels <- paste(paste("'", matrix.names, "'", sep = ""), collapse = ",")
row.labels <- paste(paste("'", matrix.names, "'", sep = ""), collapse = ",")
html.report <- gsub("--col_label--", col.labels, html.report)
html.report <- gsub("--row_label--", row.labels, html.report)

html.body.size <- 200 + left + (length(matrix.names)*20) + 30
html.report <- gsub("--body--", html.body.size, html.report)

html.report <- gsub("--function_profile_clusters--", all.cluster.functions, html.report)
html.report <- gsub("--button_clusters--", all.cluster.buttons, html.report)
html.report <- gsub("--function_profile_clusters_points--", all.cluster.functions.cov, html.report)
html.report <- gsub("--cluster_points--", all.cluster.buttons.cov, html.report)

## Insert the motif names (to hide/show all)
## They are inserted in the JQuery section
all.motifs <- paste(paste("'", all.motifs, "'", sep = ""), collapse = ",")
html.report <- gsub("--all--", all.motifs, html.report)

## Export the report
html.report.file <- paste(basename, "_scan_profile_report.html", sep = "")
write(html.report, file = html.report.file)

# grep -v '^;' matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments.tab | awk -F '\t'  ' $2!="limit" && ($11 >= 4) {print $1"\t"$3"\t"($6+$5)/2"\t"$9} ' > matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments_PARSED.tab
# cat /home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/enrichment_by_scan/Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " matrix.scan.active = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; matrix.scan.inactive = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; p.val = '1e-4'; bin = '50'; pdf.file = './test.pdf'"