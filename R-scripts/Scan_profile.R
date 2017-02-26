## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}
dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")

## Load some libraries
source(file.path(dir.rsat, 'R-scripts/config.R'))

#############################
## Step 1
## Load required libraries
#############################
required.packages <- c("RColorBrewer",
                       "gplots",
                       "png",
                       "amap",
                       "ggplot2",
                       "grid",
                       "zoo",
                       "reshape2",
                       "plyr",
                       "gridExtra",
                       "cowplot",
                       "dynamicTreeCut")

# sapply(required.packages, function(x){library(x, character.only = TRUE)})

## List of RSAT-specific packages to be compiled on the server
for (pkg in c(required.packages)) { #required.packages.bioconductor
  suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE, lib.loc=c(dir.rsat.rlib, .libPaths())))
}

#################################################################################################
## Functions
create.html.tab <- function(tab, img = 0, plot = 0, link.text.covered = 0, link.text.not.covered = 0){
  
  full.tab <- NULL
  head.tab <- "<div id='individual_motif_tab' style='width:1500px;display:none' class='tab div_chart_sp'><p style='font-size:12px;padding:0px;border:0px'><b>Individual Motif View</b></p><table id='Motif_tab' class='hover compact stripe' cellspacing='0' width='1190px' style='padding:15px;align:center;'><thead><tr><th class=\"tab_col\"> Motif_ID </th><th class=\"tab_col\"> Motif_name </th> <th class=\"tab_col\"> P-value </th> <th class=\"tab_col\"> E-value </th> <th class=\"tab_col\"> Significance </th> <th class=\"tab_col\"> FDR </th> <th class=\"tab_col\"> Nb of hits </th><th class=\"tab_col\"> Nb of sequences </th><th class=\"tab_col\">Fraction of sequences</th><th class=\"tab_col\">Scan pvalue</th><th class=\"tab_col\"> Chi-squared</th><th class=\"tab_col\"> Chi applicability</th><th class=\"tab_col\">Profile cluster</th> <th class=\"tab_col\"> Profile </th> <th class=\"tab_col\"> TFBSs </th><th class=\"tab_col\"> TFBSs per seq </th> <th class=\"tab_col\"> Logo </th> <th class=\"tab_col\"> Logo (RC) </th> <th class=\"tab_col\"> Covered sequences </th> <th class=\"tab_col\"> Not Covered sequences </th> </tr></thead><tbody>"
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

###################################################
## Step 2: Read arguments from the command line.
##
## Arguments passed on the command line
## will over-write the default arguments
## specified above.
############################################
args <- commandArgs(trailingOnly=TRUE)
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print.formats <- c("pdf", "jpeg")
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
logo.folder <- paste(basename(prefix), "_logos", sep = "")

## Heatmap dendogram position
if (heatmap.dendo == "show"){
  heatmap.dendo <- "row"
} else if(heatmap.dendo == "hide"){
  heatmap.dendo <- "none"
}
print(heatmap.dendo)

## Create a file to store the resulting tables
covered.tables.dir <- paste(basename(prefix), "_covered_sequences_info", sep = "")
dir.create(covered.tables.dir, showWarnings = FALSE)

## Induced
# matrix.scan.file <- "/home/jaime/Desktop/induced/position_scan_CapStarrSeq_K562_IFN_induced_matrix_scan_results_PARSED.tab"
# sequence.names.file <- "/home/jaime/Desktop/induced/position_scan_CapStarrSeq_K562_IFN_induced_matrix_scan_sequence_names.tab"
# ID.to.names.correspondence.tab <- "/home/jaime/Desktop/induced/position_scan_CapStarrSeq_K562_IFN_induced_TF_ID_name_correspondence.tab"
# bin <- 50
# seq.length <- 2000

## SPPS
# matrix.scan.file <- "/home/jaime/Desktop/position_scan_diff/SPPS_discovered_motifs_diff_mode_matrix_scan_results_PARSED.tab"
# sequence.names.file <- "/home/jaime/Desktop/position_scan_diff/SPPS_discovered_motifs_diff_mode_matrix_scan_sequence_names.tab"
# ID.to.names.correspondence.tab <- "/home/jaime/Desktop/position_scan_diff/SPPS_discovered_motifs_diff_mode_TF_ID_name_correspondence.tab"
# bin <- 50
# seq.length <- 600


## SPPS
# matrix.scan.file <- "/home/jcastro/Desktop/SPPS_discovered_motifs_matrix_scan_results_PARSED.tab"
# sequence.names.file <- "/home/jcastro/Desktop/SPPS_discovered_motifs_matrix_scan_sequence_names.tab"
# ID.to.names.correspondence.tab <-"/home/jcastro/Desktop/SPPS_discovered_motifs_TF_ID_name_correspondence.tab"
# bin <- 50
# seq.length <- 600

####################################
## Step 3: Read matrix-scan table
verbose(paste("Reading matrix-scan results table"), 1)
matrix.scan.results <- read.csv(file = matrix.scan.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results) <- c("seq_id", "ft_name", "bspos", "Pval", "Weight")

## Remove duplicated rows
## One row can be duplicated because some motifs are palindromic and they could have
## an identical match in the same position for the F and R strand
matrix.scan.results <- matrix.scan.results[!duplicated(matrix.scan.results), ]

## Filter by p-value
p.val <- as.numeric(p.val)
matrix.scan.results <- subset(matrix.scan.results, Pval <= p.val)

# matrix.scan.results$ft_name <- gsub(":", "_", matrix.scan.results$ft_name)

#######################################
## Step 4: Read sequence names table
verbose(paste("Reading sequence names table"), 1)
sequence.names.tab <- read.csv(file = sequence.names.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(sequence.names.tab) <- c("seq_id")
total.scanned.sequences <- length(as.vector(sequence.names.tab$seq_id))
scanned.sequences <- unique(as.vector(sequence.names.tab$seq_id))

###################################
## Calculate the sequence limits
seq.length <- as.numeric(seq.length)
limits <- seq.length/2

## Define the bins size to be tested
## Min:5 - Max:100
## All the bins that will be a divisor of the sequences length between this range will be considered
bins <- seq(5, seq.length/2, 1)
bins <- bins[(seq.length/2) %% bins == 0]
bins <- as.numeric(bin)

#################################################################
## Step 5: Create the column -log10(pvalue)                    ##
## Assign a class to each p-value, and one color to each class ##
## To later plot the qualitative distribution of TFBSs         ##
#################################################################
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

## Define the positional classes
matrix.scan.results$position_class<- ceiling(((matrix.scan.results$bspos + limits)/bins))*bins

######################################
## Get the sequences + motif name's
seq.id <- unique(as.vector(matrix.scan.results$seq_id))
matrix.names <- unique(as.vector(matrix.scan.results$ft_name))
nb.motifs <- length(matrix.names)

################################
## Step 6: Load the motif IDs
ID.names.tab <- ID.to.names.correspondence.tab
ID.names <- read.table(ID.names.tab, sep = "\t")
ID.names$V1 <- gsub(":", "_", ID.names$V1)
ID.names$V2 <- gsub(":", "_", ID.names$V2)

setwd(results.folder)

##################################################################
## Step 7: Plot the distribution of TFBSs at different p-values ##
##################################################################
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

# thr <- sapply(1:nb.motifs, function(m){
# 
#   ## Get the matrix name
#   matrix.query <- matrix.names[m]
#   matrix.query.name <- as.vector(ID.names[,2][which(ID.names[,1] == matrix.query)][1])
# 
#   ## Get the sub-table with the hits of the query matrix
#   matrix.query.selection <- matrix.scan.results[matrix.scan.results$ft_name == matrix.query,]
#   matrix.query.selection$bspos <- matrix.query.selection$bspos + limits
#   matrix.query.classes <- sort(unique(matrix.query.selection$Pval.class.letter))
# 
#   ##
#   nb.hits.per.sequence <- table(as.vector(matrix.query.selection$seq_id))
#   nb.hits.per.sequence.range <- range(nb.hits.per.sequence)
#   min.nb.hits <- nb.hits.per.sequence.range[1]
#   max.nb.hits <- nb.hits.per.sequence.range[2]
# 
#   nb.hits.per.sequence <- table(nb.hits.per.sequence)
#   no.hit.nb <- total.scanned.sequences - sum(nb.hits.per.sequence)
#   names(no.hit.nb) <- 0
#   nb.hits.per.sequence <- append(nb.hits.per.sequence, no.hit.nb, after = 0)
# 
#   ## Generate GGplot for number of TFBSs per sequence
#   TFBSs.per.seq.file <- paste(basename(prefix), "_TFBSs_per_seq/", matrix.query, "_TFBSs_per_seq", sep = "")
# 
#   aa <- data.frame(nb.seq = nb.hits.per.sequence, nb.hits = as.numeric(names(nb.hits.per.sequence)))
#   ggplot(data=aa, aes(y = nb.seq, x=nb.hits)) +
#     geom_bar(aes(fill=nb.seq), stat="identity", position=position_dodge()) +
#     labs(title="Number of hits", y = "Number of sequences", x="Number of TFBSs") +
#     geom_text(aes(label=nb.seq), vjust=-0.15, color="black", size=3)
# 
#   suppressMessages(ggsave(paste(TFBSs.per.seq.file, ".jpeg", sep = "")))
#   suppressMessages(ggsave(paste(TFBSs.per.seq.file, ".pdf", sep = "")))
# 
#   ## Get the number of putative TFBSs and the number of sequences with
#   ## at least one match of the query matrix
#   nb.TFBSs <- dim(matrix.query.selection)[1]
#   nb.seq <- length(as.vector(unique(matrix.query.selection$seq_id)))
# 
#   TFBSs.pval.distribution.file <- paste(basename(prefix), "_TFBSs_pval_distribution/", matrix.query, "_TFBSs_pval_classes", sep = "")
# 
#   ################################
#   #Create a custom color scale
#   myColors <- colorRampPalette(brewer.pal(nb.color.classes, "YlGnBu"), space="Lab")(length(classes.pval.letters))
#   names(myColors) <- classes.pval.letters
# 
#   ## Range of p-values for the query motif
#   pval.class.matrix.query <- sort(unique(matrix.query.selection$Pval.class))
# 
#   ## Insert logo
#   logo.file <- paste(logo.folder, "/", matrix.query, "_logo.png", sep = "")
#   logo <- readPNG(logo.file)
#   logo.roster <- rasterGrob(logo, interpolate = TRUE)
# 
#   ## X position of plot annotations
#   text.xmax <- min(matrix.query.selection$bspos) + max(matrix.query.selection$bspos) / 4
#   text.center <- (min(matrix.query.selection$bspos) - text.xmax)*2
# 
#   ggplot(matrix.query.selection, aes(x=bspos, y=Pval.minlog10)) +
#     # ylim(c(min(matrix.scan.results$Pval.minlog10), max(matrix.scan.results$Pval.minlog10))) +
#     geom_point(aes(colour = Pval.class.letter), shape = 18, size = 4, stroke = 0.5) +
#     # geom_rug(position='jitter') +
#     labs(title=paste("Qualitative distribution of ", matrix.query, " TFBSs", sep = ""), y = "-log10(P-value)", x = "Position") +
#     scale_colour_manual(name = "-log10(P-value)",values = myColors, labels = paste(">", pval.class.matrix.query, sep = "")) +
#     theme_minimal() +
#     annotate("text", x = -limits + ((limits*2)/10), y = max.pval.minus.log10 - 0.25, label = paste("Nb of TFBSs: ", nb.TFBSs, sep = ""), size = 4, hjust = 0) +
#     annotate("text", x = -limits + ((limits*2)/10), y = max.pval.minus.log10 - 0.55, label = paste("Nb of sequences: ", nb.seq, sep = ""), size = 4, hjust = 0) +
#     scale_y_continuous(breaks=seq(round(min(matrix.scan.results$Pval.minlog10)),round(max(matrix.scan.results$Pval.minlog10)),1)) +
#     annotation_custom(logo.roster, xmax = limits - (limits/3), xmin = limits - 5, ymin = max.pval.minus.log10 - 1, ymax = max.pval.minus.log10 - 0.05)
# 
#   suppressMessages(ggsave(paste(TFBSs.pval.distribution.file, ".pdf", sep = "")))
#   suppressMessages(ggsave(paste(TFBSs.pval.distribution.file, ".jpeg", sep = "")))
# })
# rm(thr)


#################################################################################
## Step :  ##
#################################################################################
# verbose(paste("Creating boxplot with the distribution of TFBSs Weights at each bin"), 1)
# 
# boxplot.weights.dir <- paste(basename(prefix), "_Boxplot_Weight", sep = "")
# dir.create(boxplot.weights.dir, showWarnings = FALSE, recursive = TRUE)
# 
# matrix.scan.results <- ddply(matrix.scan.results, .(ft_name, position_class), mutate, median.w = round(median(Weight), digits = 1))
# matrix.scan.results$median.w[matrix.scan.results$median.w < 0] <- 0
# max.W <- max(matrix.scan.results$Weight)
# # aa <- matrix.scan.results
# # matrix.scan.results <- aa
# 
# median.W.values <- seq(0, max(as.vector(unique(round(sort(matrix.scan.results$median.w))))), by = 1)
# median.W.classes <- length(median.W.values)
# 
# ## Calculate a color palette for the Weight classes
# myColors <- colorRampPalette(brewer.pal(9, "YlGnBu"), space="Lab")(median.W.classes)
# 
# matrix.scan.results$median.w <- as.factor(matrix.scan.results$median.w)
# matrix.scan.results$position_class <- as.factor(matrix.scan.results$position_class)
# 
# # matrix.scan.results$class.color <- as.vector(myColors[matrix.scan.results$mean.w])
# 
# thr <- sapply(1:nb.motifs, function(m){
# 
#   print(matrix.query)
#   
#   ## Get the matrix name
#   matrix.query <- matrix.names[m]
#   matrix.query.name <- as.vector(ID.names[,2][which(ID.names[,1] == matrix.query)][1])
# 
#   ## Get the sub-table with the hits of the query matrix
#   matrix.query.selection <- subset(matrix.scan.results, ft_name == matrix.query)
# 
#   ## Calculate the mean W per bin
#   median.w.per.bin <- ddply(matrix.query.selection, "position_class", summarize, median.w = median(Weight))[,2]
#   median.w <- median(matrix.query.selection$Weight)
# 
#   ## Get the possitional classes
#   pos.class <- unique(as.vector(matrix.query.selection$position_class))
#   pos.class.sep.list <- split(matrix.query.selection, f = as.factor(matrix.query.selection$position_class))
# 
#   ## Calculate a T-test of the W ditribution on each bin vs the overal W distribution
# 
#   co <- 0
#   student.pvalues <- sapply(pos.class.sep.list, function(l){
# 
#     co <<- co + 1
# 
#     if( dim(l)[1] < 2){
#       NA
#     } else {
#       student <- t.test(x = l$Weight,
#                         y = matrix.query.selection$Weight)
#       student[["p.value"]]
#     }
# 
#   })
#   student.pvalues.corrected <- p.adjust(student.pvalues, method = "bonferroni")
# 
#   ## Assign colors
#   median.w <- round(median.w.per.bin)
#   median.w[median.w < 0] <- 0
#   median.w.to.color <- as.vector(myColors[(median.w+1)])
# 
#   ## Create the P-val_t_test column
#   list.counter <- 0
#   student.pvalues.corrected.vector <- as.vector(unlist(sapply(pos.class.sep.list, function(l){
#     list.counter <<- list.counter + 1
#     rep(student.pvalues.corrected[list.counter], times = nrow(l))
#   })))
# 
#   student.pvalues.corrected.pretty <- prettyNum(student.pvalues.corrected.vector, scientific=TRUE, digits = 1)
#   student.pvalues.corrected.pretty <- as.character(student.pvalues.corrected.pretty)
#   matrix.query.selection$P_val_t_test <- student.pvalues.corrected.pretty
# 
#   weight.boxplot.bin.file <- paste(basename(prefix), "_Boxplot_Weight/", matrix.query, "_Weight_distribution_per_bin", sep = "")
# 
#   boxplot.bins <- ggplot(matrix.query.selection, aes(x=position_class, y=Weight, group = position_class, fill = position_class)) +
#     geom_boxplot(notch = FALSE) +
#     guides(fill=FALSE) +
#     labs(title=paste("Distribution of TFBSs Weigths (separated by bins) for ", matrix.query, sep = ""), y = "", x = "Position") +
#     ylim(c(-1, max.W)) +
#     scale_fill_manual(values=(median.w.to.color)) +
#     geom_text(data = matrix.query.selection, aes(x = position_class, y = 0, label = median.w), size = 3.5, colour = "#424242") +
#     geom_text(data = matrix.query.selection, aes(x = (length(median.w.to.color)/2), y = 0.55, label = "Median Weight"), size = 5, colour = "#424242") +
#     geom_text(data = matrix.query.selection, aes(x = position_class, y = -1, label = P_val_t_test), size = 3.5, colour = "#424242") +
#     geom_text(data = matrix.query.selection, aes(x = (length(median.w.to.color)/2), y = -0.65, label = "Corrected P-value (t-test)"), size = 5, colour = "#424242") +
#     theme_minimal()
# 
#   boxplot.overall <- ggplot(matrix.query.selection, aes(x=ft_name, y=Weight, group = ft_name, fill = ft_name)) +
#     geom_boxplot( notch = FALSE) +
#     ylim(c(-1, max.W)) +
#     geom_text(data = matrix.query.selection, aes(x = 1, y = 0, label = median(matrix.query.selection$Weight)), size = 3.5, colour = "#424242") +
#     geom_text(data = matrix.query.selection, aes(x = 1, y = 1, label = "Median\nWeight"), size = 5, colour = "#424242") +
#     guides(fill=FALSE) +
#     labs(title=paste("Distribution of TFBSs Weigths\n", matrix.query, sep = ""), y = "Weight", x = "") +
#     theme_minimal()
# 
#   xx <- plot_grid(boxplot.overall, boxplot.bins, labels=c("", ""), ncol = 2, nrow = 1, rel_widths = c(0.2, 1.75))
# 
#   save_plot(filename = paste(weight.boxplot.bin.file, ".pdf", sep = ""), plot = xx, ncol = 2, base_aspect_ratio = 5, base_width = 15, base_height = 10)
#   save_plot(filename = paste(weight.boxplot.bin.file, ".jpeg", sep = ""), plot = xx, ncol = 2, base_aspect_ratio = 5, base_width = 15, base_height = 10)
# 
# })
# rm(thr)

#################################################################################
## Step 8: Calculate the raw counts of TFBSs on each position of the sequences ##
#################################################################################
verbose(paste("Calculating the raw counts of TFBSs in the sequences"),1)
matrix.scan.results$bspos.left <- matrix.scan.results$bspos + seq.length
raw.counts.all.motifs <- list()
p.counter <- 0
tested.pvalues <- unique(sort(matrix.scan.results$Pval.class))[1:4] -0.5
for(p in tested.pvalues){
  
  p.counter <- p.counter + 1
  
  raw.counts.all.motifs[[p.counter]] <- sapply(matrix.names, function(m){
    
    matrix.scan.results.pval.filter <- subset(matrix.scan.results, Pval.class >= p)
    
    sites.m <- as.vector(matrix.scan.results[which(matrix.scan.results.pval.filter$ft_name == m), "bspos.left"])
    raw.counts.m <- tabulate(sites.m, nbins = seq.length)
    # return(raw.counts.m)
    raw.counts.m
  })
  raw.counts.all.motifs[[p.counter]] <- as.data.frame(raw.counts.all.motifs[[p.counter]])
  
}

#######################################################
## Step 9: Calculate the raw counts of TFBSs per bin ##
#######################################################

## Define the bin 'names'
xlab <- vector("list", length(bins))
xlab <- sapply(bins, function(b){
  x.neg <- seq(-(seq.length/2), -(b) , b)
  append(x.neg, rev(-x.neg))
})
if(length(bins) == 1){
  xlab <- list(as.vector(xlab))
}

## For each bin, sums the number of TFBSs (function rollaply)
## The results are stored in a list where each element correspond to the counts of a given bin size
counts.per.bin <- list()
b <- 1 
p.counter <- 0
for(p in tested.pvalues){
  
  p.counter <- p.counter + 1
  print(p.counter)
  
  # verbose(paste("Calculating the raw counts of TFBSs in the sequences - Bin of size:", bins[b]),1)
  counts.per.bin.query <- rollapply(raw.counts.all.motifs[[p.counter]], bins[b], sum, na.rm = TRUE, by = bins[b], partial = TRUE, align="left")
  colnames(counts.per.bin.query) <- matrix.names
  rownames(counts.per.bin.query) <- xlab[[b]]
  counts.per.bin[[p.counter]] <- counts.per.bin.query
}


##############################################
## Step 10: Chi-squared calculation section ##
##############################################
counts.per.bin <- lapply(counts.per.bin, t)

pval.sig.hash <- list("3" = 1e-4, "3.5" = 5e-4, "4" = 1e-4, "4.5" = 5e-5, "5" = 1e-5, "5.5" = 5e-6, "6" = 1e-6)
all.chi.results <- data.frame()
## Iterate over the motifs
thrash <- sapply(matrix.names, function(m){
  
  print(m)
  p.counter <- 0
  
  ## Iterate over the count per bin tables at different p-values
  motif.chi.statistics <- sapply(counts.per.bin, function(cpb){
    
    p.counter <<- p.counter + 1
    print(p.counter)
    
    if(sum(cpb[m,]) > 0){
      
      motif.chi <- chisq.test(cpb[m,], correct = TRUE)
      
      chi.applicability.condition <- NULL
      if(unique(motif.chi[[7]]) >= 5){
        chi.applicability.condition <- "Yes"
      } else {
        chi.applicability.condition <- "No"
      }
      
      selected.fields <- c(m ,motif.chi[[1]], motif.chi[[2]], motif.chi[[3]], bins[1], chi.applicability.condition, tested.pvalues[p.counter])
      names(selected.fields) <- NULL
      selected.fields
      
    } else {
      selected.fields <- c(NA, NA, NA, NA, NA, NA, NA)
      names(selected.fields) <- NULL
      selected.fields
    }
  })
  motif.chi.statistics <- data.frame(t(motif.chi.statistics))
  colnames(motif.chi.statistics) <- c("Motif", "Chi", "DF", "Pvalue", "Bin_size", "Chi_app", "Scan_pval")
  
  ## Convert the columns to vectors
  motif.chi.statistics$Chi_app <- as.vector(motif.chi.statistics$Chi_app)
  motif.chi.statistics$Scan_pval <- as.vector(motif.chi.statistics$Scan_pval)
  motif.chi.statistics$Pvalue <- as.vector(motif.chi.statistics$Pvalue)
  
  ## Select the best scan P-value: the program minimizes the p-value when 
  ## the Chi squared applicability condition is satisfied.
  ## If it's nos satisfied, selects the higher scan Pvalue.
  
  ## Rows with Chi applicablity criterion
  chip.app.nrow <- nrow(subset(motif.chi.statistics, Chi_app == "Yes"))
  if(chip.app.nrow > 0){
    
    motif.chi.statistics <- subset(motif.chi.statistics, Chi_app == "Yes")
    row.min.pval <- which(motif.chi.statistics$Pvalue == min(motif.chi.statistics$Pvalue))[1]
    
    motif.chi.statistics <- motif.chi.statistics[row.min.pval,]
    highest.pval.sig <- as.vector(motif.chi.statistics[1,7])
    highest.pval <- pval.sig.hash[[highest.pval.sig]]
    motif.chi.statistics[1,7] <- highest.pval
    
  } else {
    highest.pval.sig <- as.vector(motif.chi.statistics[1,7])
    highest.pval <- pval.sig.hash[[highest.pval.sig]]
    motif.chi.statistics[1,7] <- as.vector(highest.pval)
    
    motif.chi.statistics <- motif.chi.statistics[1,]
  }
  
  all.chi.results <<- rbind(all.chi.results, motif.chi.statistics)
  
})
rownames(all.chi.results) <- NULL
rownames(all.chi.results) <- as.vector(all.chi.results$Motif)

##################################
## Enrichment in central region ##
##################################
# center.region.up <- 50
# center.region.dw <- 50 
# center.raw.counts.up <- (seq.length/2) - center.region.up
# center.raw.counts.dw <- (seq.length/2) + center.region.dw
# 
# mean.counts <- floor(colMeans(raw.counts.all.motifs))
# counts.in.center <- t(raw.counts.all.motifs[center.raw.counts.up:center.raw.counts.dw, ])
# 
# motif <- "cluster_11"
# motif <- "cluster_10"
# motif <- "cluster_1"
# 
# sapply(rownames(counts.in.center), function(motif){
#   
#   n.pos <- ncol(counts.in.center)
#   expected.counts <- as.vector(n.pos * (mean.counts[motif]))
#   sum.counts <- sum(counts.in.center[motif,])
#   
#   print(paste(sum.counts, expected.counts, sep = " - "))
#   
#   if(sum.counts >= expected.counts){
#     ppois(q = sum.counts, lambda = expected.counts, lower.tail = TRUE, log.p = TRUE)
#   } else {
#     ppois(q = sum.counts, lambda = expected.counts, lower.tail = FALSE, log.p = TRUE)
#   }
#   
# 
#   ppois(q = sum.counts, lambda = expected.counts, lower.tail = TRUE)
#   
# })
# 
# counts.per.bin.query <- rollapply(raw.counts.all.motifs, bins[b], sum, na.rm = TRUE, by = bins[b], partial = TRUE, align="left")
# 

#############################################################
## Step 11: Calculate number of hits and matched sequences ##
## Export files with matched and unmatched sequences       ##
#############################################################

## Count the number of matched sequences
verbose(paste("Counting the number of hits per motif"),1)
matched.sequences <- sapply(matrix.names, function(m){
  
  motif.scan.pval <- as.numeric(all.chi.results[m,"Scan_pval"])
  
  matrix.scan.results.subset <- subset(matrix.scan.results, ft_name == m & Pval <= motif.scan.pval)
  length(unique(as.vector(matrix.scan.results.subset$seq_id)))
})


## Count matches for each motif
verbose(paste("Counting the number of matched sequences per motif"),1)
matches.per.motif <- sapply(matrix.names, function(m){
  
  motif.scan.pval <- as.numeric(all.chi.results[m,"Scan_pval"])
  
  matrix.scan.results.subset <- subset(matrix.scan.results, ft_name == m & Pval <= motif.scan.pval)
  nrow(matrix.scan.results.subset)
})


## Matched and unmatched sequences
verbose(paste("Exporting names of matched and not-matched sequences per motif"),1)
covered.seq <- rep(0, nb.motifs)
not.covered.seq <- rep(0, nb.motifs)
names(covered.seq) <- matrix.names
names(not.covered.seq) <- matrix.names
list.counter <- 0
thrash <- sapply(matrix.names, function(m){
  
  list.counter <<- list.counter + 1
  
  motif.scan.pval <- as.numeric(all.chi.results[m,"Scan_pval"])
  matrix.scan.results.subset <- subset(matrix.scan.results, ft_name == m & Pval <= motif.scan.pval)
  
  covered.seq.query <- unique(as.vector(matrix.scan.results.subset$seq_id))
  not.covered.seq.query <- setdiff(scanned.sequences, covered.seq)
  
  covered.seq[list.counter] <<- file.path(covered.tables.dir, paste(m, "_covered_sequences_IDs.tab", sep = ""))
  not.covered.seq[list.counter] <<- file.path(covered.tables.dir, paste(m, "_not_covered_sequences_IDs.tab", sep = ""))

  covered.sequences.table <- as.data.frame(covered.seq.query)
  not.covered.sequences.table <- as.data.frame(not.covered.seq.query)
  
  # write.table(covered.sequences.table, file = covered.seq[list.counter], sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  # write.table(not.covered.sequences.table, file = not.covered.seq[list.counter], sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
})
rm(thrash)

#############################################
## Step 12: Calculate the Frequency tables ##
#############################################
dir.create(paste(basename(prefix), "_tables/", sep = ""), recursive = TRUE, showWarnings = FALSE )
freq.per.bin <- NULL
list.counter <- 0
counts.tab.file <- NULL
density.tab.file <- NULL
counts.per.bin <- sapply(matrix.names, function(m){
  
  pval.scan.index <- as.numeric(which(pval.sig.hash == as.numeric(all.chi.results[m,"Scan_pval"])))[1]
  cpb <- counts.per.bin[[pval.scan.index]][m,]
  cpb
})  

## Count the total of TFBSs
total.hits <- colSums(counts.per.bin)
  
## Divide the sum at each bin by the total
freq.per.bin <- as.data.frame(t(counts.per.bin)/total.hits)
  
##########################################
## Export Counts and Frequencies tables
density.tab.file <<- paste(basename(prefix), "_tables/density_per_bin_profiles_bin", bins[1], ".tab", sep = "")
write.table(freq.per.bin, file = density.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  
counts.tab.file <<- paste(basename(prefix), "_tables/counts_per_bin_profiles_bin", bins[1], ".tab", sep = "")
write.table(t(counts.per.bin), file = counts.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  
freq.log2.ratio <- list()
freq.log2.ratio[[1]] <- -log2(t(counts.per.bin)/rowMeans(t(counts.per.bin)))
freq.log2.ratio[[1]][is.infinite(freq.log2.ratio[[1]])] <- 0


########################################################################
## Step 13: calculate Q-value, E-value, and significance -log(pvalue) ##
## Re-order the table                                                 ##
########################################################################

all.chi.results$Chi <- round(as.numeric(as.vector(all.chi.results$Chi)), digits = 1)

## Calculate E-value: the correction for multitesting factor is the number of motifs
## that satisfied the Chi-applicability condition
chi.app.nb <- sum(all.chi.results$Chi_app == "Yes")
all.chi.results$Evalue <- as.numeric(as.vector(all.chi.results$Pvalue)) * chi.app.nb

## Calculate Significance
all.chi.results$Sig <- round(-log10(all.chi.results$Evalue), digits = 1)
  
## Calculate q-values
## This step is executed once all the p-values were calculated
## The variable with class 'qvalue' is stored to its further exportation
pp <- as.numeric(as.vector(all.chi.results$Pvalue))
features.qvalues <- p.adjust(pp, method = "BH")
all.chi.results$Qvalue <- prettyNum(features.qvalues, scientific=TRUE, digits = 2)

## P-val -> Pretty number
all.chi.results$Pvalue <- round(as.numeric(as.vector(all.chi.results$Pvalue)), digits = 100000000000000000)
all.chi.results$Pvalue <- prettyNum(all.chi.results$Pvalue, scientific=TRUE, digits = 2)

## E-val -> Pretty number
all.chi.results$Evalue <- prettyNum(all.chi.results$Evalue, scientific=TRUE, digits = 2)

## Add the matrix Id as feature ID
all.chi.results$feature <- rownames(all.chi.results)

## Add the coverage column
all.chi.results$coverage <- round(matched.sequences/total.scanned.sequences, digits = 2)
all.chi.results$sequences <- matched.sequences
  
## Add the number of matches per motif
all.chi.results$hits <- matches.per.motif

all.chi.results$cov <- covered.seq
all.chi.results$not_cov <- not.covered.seq

features.table <- all.chi.results[, c(1,2,3,6,9,4,8,10,12,13,14,7,15,16)]
colnames(features.table) <- c("Feature", "Chi_squared", "Degrees", "Chi_app", "Sig", "P_val" ,"E_val", "Q_val", "Coverage", "Sequences", "Nb_hits", "Pval_scan", "Covered_seq", "Not_covered_seq")
  
feature.attributes.file <- paste(basename, "_attributes.tab", sep = "")
write.table(features.table, file = feature.attributes.file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


###################################
## Step 14: Calculate BG matches ##
###################################
verbose(paste("Calculating the BG number of matches per motif"),1)

list.counter <- 1
counts.norm <- list()
counts.norm[[list.counter]] <- -log2(t(counts.per.bin)/rowMeans(t(counts.per.bin)))

# thrash <- sapply(bins, function(b){
# 
#   list.counter <<- list.counter + 1
# 
#   ## Mean of TFBSs per bin
#   mean.counts.all.motifs.bins <- rollapply(raw.counts.all.motifs, b, mean, na.rm = TRUE, by = b, partial = TRUE, align="left")
# 
#   ## Mean of TFBSs per bin
#   sum.counts.all.motifs.bins <- rollapply(raw.counts.all.motifs, b, sum, na.rm = TRUE, by = b, partial = TRUE, align="left")
# 
#   df_by_bins_norm <- apply(sum.counts.all.motifs.bins, 1, function(r){
#     r/mean(r)
#   })
#   df_by_bins_norm <- -log2(df_by_bins_norm)
#   df_by_bins_norm[is.infinite(df_by_bins_norm)] <- 0
# 
#   # ## Sum the means per bin
#   # all.columns.trimmed.means <- colSums(mean.counts.all.motifs.bins, na.rm = TRUE)
#   #
#   # ## The median of the means will be considered as the "basal" number of hits
#   # for (m in 1:ncol(mean.counts.all.motifs.bins)){
#   #   # row.mean <- mean(mean.counts.all.motifs.bins[,m],trim=0.3)
#   #   row.mean <- median(mean.counts.all.motifs.bins[,m])
#   #   all.columns.trimmed.means[m] <- row.mean
#   # }
#   #
#   # df_by_bins_norm = t(t(mean.counts.all.motifs.bins)-all.columns.trimmed.means)
# 
#   counts.norm[[list.counter]] <<- df_by_bins_norm
# })


# counts.norm[[1]] <- -log2(t(counts.per.bin[[1]])/rowMeans(t(counts.per.bin[[1]])))

## To check
# z <- t(counts.per.bin[[1]])[1,]
# zz <- -log2(z/mean(z))
# plot(zz,
#      type = "l")

######################################
## Step 15: Define profile clusters ##
######################################
verbose(paste("Cluster the positional profiles"),1)
list.counter <- 1
cluster.tree.profiles.palette <- vector("list", length(counts.norm))
color.clusters.tree.profiles <- vector("list", length(counts.norm))
profile.clusters.names <- vector("list", length(counts.norm))
tree.profiles <- vector("list", length(counts.norm))
cluster.profiles.motifs <- NULL
cluster.profiles.motif.names <- NULL

freq.per.bin.norm <- list()
freq.per.bin.norm[[1]] <- freq.per.bin/total.scanned.sequences
## Calculate the highest frequency of TFBSs per each list
max.y <- max(unlist(freq.per.bin.norm[[1]]), na.rm = TRUE)
# max.y <- lapply(freq.per.bin.norm[[1]], max, na.rm = TRUE)

cn <- freq.per.bin
  
tree.profiles[[list.counter]] <- hclust(Dist(cn, method = "correlation"), method = "complete")
clusters.tree.profiles <- cutreeDynamic(tree.profiles[[list.counter]], minClusterSize = 1, method = "tree")
names(clusters.tree.profiles) <- matrix.names

## Generate a color palette
nb.profile.clusters <- length(unique(clusters.tree.profiles))
cluster.tree.profiles.palette[[list.counter]] <- colorRampPalette(brewer.pal(9, "Set1"), space="Lab")(nb.profile.clusters)

## Assign a different color to each cluster
color.clusters.tree.profiles[[list.counter]] <- as.vector(sapply(clusters.tree.profiles, function(color){
  cluster.tree.profiles.palette[[list.counter]][color]
}))
  
############################################
## Fill the Profile_cluster attribute 
## Add a new column to all.chi.results df
profile.clusters.names[[list.counter]] <- paste("Profile_cluster", clusters.tree.profiles, sep = "_")

## Get the member motif IDs of each cluster
cluster.profiles.counter <- 0
thrash <- sapply(1:nb.profile.clusters, function(cl){
  
  cluster.profiles.counter <<- cluster.profiles.counter + 1
  cluster.profiles.motifs[[list.counter]][[cluster.profiles.counter]] <<- names(which(clusters.tree.profiles == cluster.profiles.counter))
  
  ## Get the motif name
  cluster.profiles.motif.names[[list.counter]][[cluster.profiles.counter]] <<- as.vector(
    sapply(cluster.profiles.motifs[[list.counter]][[cluster.profiles.counter]], function(n){
      ID.names[which(ID.names[,1] == n),2]
    })
  )
})


##########################################################
## Step 16: Draw Profiles heatmap                       ##
## Shows the frequencies of hits per bin for each motif ##
##########################################################
verbose(paste("Drawing Heatmap profiles"),1)
list.counter <- 0
thrash <- lapply(freq.log2.ratio, function(cn){
  
  list.counter <<- list.counter + 1
  
  ## Profile Heatmap Color palette (user-defined)
  rgb.palette <- rev(colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")(heatmap.color.classes))
  
  ## Remove noisy extremities
  df.limits <- round(range(cn))
  df.limits.abs <- abs(df.limits)
  
  # if (df.limits.abs[1] > df.limits.abs[2]){
  #   cn[cn < (-df.limits[2])] <- (-df.limits[2])
  #   cn[cn > df.limits[2]] <- df.limits[2]
  # } else {
  #   cn[cn > (-df.limits[1])] <- (-df.limits[1])
  #   cn[cn < df.limits[1]] <- df.limits[1]
  # }
  
  # cn.t <- t(cn)
  cn.t <- cn
  colnames(cn.t) <- xlab[[list.counter]]
  
  ## Print the heatmap
  heatmap.profiles <- NULL
  for(format in print.formats){
    
    profiles.heatmap.file <- paste(basename(prefix), "_profiles_heatmap.", format, sep = "")
    
    if(format == "pdf"){
      pdf(profiles.heatmap.file)
    } else if (format == "jpg"){
      jpeg(profiles.heatmap.file)
    }
    
    ## Heatmap
    heatmap.profiles <<- heatmap.2(cn.t,
                                   
                                   main = "Profile Heatmap",
                                   xlab = "Position (bp)",
                                   ylab = "Motifs",
                                   
                                   ## The order of the values is set according these dendrograms
                                   Rowv = as.dendrogram(tree.profiles[[list.counter]]),
                                   Colv = FALSE,
                                   dendrogram = "row",
                                   
                                   ## Color
                                   col = rgb.palette,
                                   
                                   ## Trace
                                   trace = "none",
                                   
                                   ## Side colors
                                   RowSideColors = color.clusters.tree.profiles[[list.counter]],
                                   
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
    thr <- dev.off()
    
    ## Export the heatmap row order in a tab file
    heatmap.row.order <- rev(heatmap.profiles[[1]])
    heatmap.row.order.names <- matrix.names[heatmap.row.order]
    heatmap.rows <- data.frame(row = heatmap.row.order, names = heatmap.row.order.names)
    heatmap.rows.file <- paste(basename, "_heatmap_row_order_bin_", bins[list.counter],"tab", sep = "")
    write.table(heatmap.rows, file = heatmap.rows.file, quote = FALSE, sep = "\t", row.names = FALSE)
    
  }
  rm(thr)
})
rm(thrash)

#######################################################################
## Step 17: Compute the XY-plot for Profile significance vs Coverage ##
#######################################################################
verbose(paste("Drawing Significance vs Coverage plot"),1)
x.y.coverage <- vector("list", length(bins))
list.counter <- 0
x.sig <- vector("list", length(bins))
y.cov <- vector("list", length(bins))
th <- sapply(bins, function(b){
  
  list.counter <<- list.counter + 1
  
  ## Calculate X-Y coordinates
  ## X = Significance
  x.sig[[list.counter]] <<- as.numeric(as.vector(features.table$Sig))
  inf.val <- which(x.sig[[list.counter]] == Inf)
  x.sig[[list.counter]][inf.val] <<- 350
  names(x.sig[[list.counter]]) <- as.vector(features.table$Feature)
  y.cov[[list.counter]] <<- as.numeric(gsub("%", "", features.table$Coverage)) *100
  names(y.cov[[list.counter]]) <- as.vector(features.table$Feature)
  
  ## X-Y plot
  plot(x.sig[[list.counter]],
       y.cov[[list.counter]],
       ylim = c(0,100),
       xlab = "Significance -log10(Corrected E-val)",
       ylab = "Sequence coverage",
       main = "Profile Significance vs Sequence coverage",
       col = ifelse((x.sig[[list.counter]] >= 20 & y.cov[[list.counter]] >= 66), "darkgreen", "gray"),
       panel.first=grid(col = "grey", lty = "solid"),
       pch = "o",
       cex = 1.5
  )
  
  ## Mark the TFBMs sattisfying the threshold
  selected.TFBMs <- features.table[which(as.vector(features.table$Sig) >= 20 & as.numeric(gsub("%", "", features.table$coverage)) >= 66), "Feature"]
  selected.TFBMs <- as.vector(selected.TFBMs)
  
  if(length(selected.TFBMs) > 0){
    
    ## Add the text to the selected TFBMs
    text(x = x.sig[[list.counter]][c(selected.TFBMs)],
         y = y.cov[[list.counter]][c(selected.TFBMs)],
         labels = selected.TFBMs,
         cex = 0.6,
         pos = 3,
         col="red")
  }
  thrash <- dev.off()
  
  ## Convert the X-Y values to the format required for C3 plot
  xx.sig <- paste("['x',", paste(round(as.vector(x.sig[[list.counter]])), collapse = ","), "],", sep = "")
  yy.cov <- paste("['y',", paste(round(as.vector(y.cov[[list.counter]])), collapse = ","), "],", sep = "")
  x.y.coverage[[list.counter]] <<- paste(xx.sig, yy.cov, collapse = "\n")
})
rm(th)

##############################################################
## Step 18: Plot each profile individually                  ##
## Require the icon of the feature (e.g. PSSM logo)         ##
##############################################################
verbose(paste("Printing all the profiles in a PDF file"), 1)
## Create folder for individual profile plots
dir.create(paste(basename(prefix), "_TFBSs_positional_profiles/", sep = ""), recursive = TRUE, showWarnings = FALSE )
list.counter <- 0
p.counter <- 0
# thr <- sapply(bins, function(b){
# 
#   list.counter <<- list.counter + 1
# 
#   thrash <- sapply(matrix.names, function(feature.query){
# 
#     p.counter <<- p.counter + 1
#     # print(p.counter)
# 
#     ## Set Output file name
#     file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", feature.query, "_positional_profile_bins_", b, sep = "")
# 
#     ## Get the coordinates
#     y.val <- as.numeric(freq.per.bin.norm[[list.counter]][feature.query,])
#     x.val <- as.numeric(xlab[[list.counter]])
# 
#     ## Load the logo
#     # logo.file <- paste(logo.folder, "/", feature.query, "_logo.png", sep = "")
#     # logo <- readPNG(logo.file)
#     # logo.roster <- rasterGrob(logo, interpolate = TRUE)
# 
#     ## Plot the profile (using ggplot2)
#     xy.df <- data.frame(x = x.val, y = y.val)
#     ggplot(xy.df, aes(x=x, y=y)) +
#       geom_line(colour = "#00BFC4", size = 3) +
#       ylim(0, max.y) +
#       labs(title=paste(feature.query, " binding profile", sep = ""), y = "Frequency of TFBSs", x = "Position") +
#       # geom_rug(position='jitter', sides="l") +
#       geom_area(fill = "#00BFC4", alpha=0.35) #+
#       #annotation_custom(logo.roster, xmax = limits, xmin = limits - sum(abs(range(xy.df$x)))/5, ymin = max.y[[list.counter]] - 0.01, ymax = max.y[[list.counter]] - 0.075)
# 
#     ## Export the file
#     suppressMessages(ggsave(paste(file.name, ".jpeg", sep = ""), plot = last_plot()))
#     suppressMessages(ggsave(paste(file.name, ".pdf", sep = ""), plot = last_plot()))
#   })
# })
# rm(thr)

##################################################################
## Step 19: Get the name of the matched sequences of each motif ##
##################################################################
covered.sequences.per.motif <- ddply(matrix.scan.results, "ft_name", summarize, x = paste(unique(seq_id), collapse = ","))
covered.sequences.per.motif <- lapply(covered.sequences.per.motif$x, function(x){
  unlist(strsplit(x, ","))
})
names(covered.sequences.per.motif) <- matrix.names

###########################################################
## Step 20: Draw Co-ocurrence heatmap (to display in D3) ##
###########################################################
# covered.seq.percentage <- NULL
# thrash <- sapply(covered.sequences.per.motif, function(m1){
#   sapply(covered.sequences.per.motif, function(m2){
#     
#     intersected.seq <- intersect(m1, m2)
#     intersected.seq.nb <- length( intersected.seq)
#     intersected.seq.per <- intersected.seq.nb/total.scanned.sequences
#     covered.seq.percentage <<- append(covered.seq.percentage, intersected.seq.per)
#   })
# })
# covered.seq.percentage <- round(covered.seq.percentage, digits = 4) * 100
# covered.seq.percentage <- matrix(covered.seq.percentage, ncol = length(covered.sequences.per.motif))
# colnames(covered.seq.percentage) <- matrix.names
# rownames(covered.seq.percentage) <- matrix.names
# 
# ## Set the colors
# coocurrence.palette <- colorRampPalette(brewer.pal(9, "YlGnBu"), space="Lab")(20)
# 
# ## Draw the co-ocurrence heatmap
# verbose(paste("Creating co-ocurrence heatmap"), 1)
# 
# heatmap.profiles <- NULL
# comp.order.list.rows <- vector("list", 4)
# comp.order.list.cols <- vector("list", 4)
# domain <- seq(from = 0, to = 100, by = 5)
# 
# ## Calculate the distance table between the co-ocurrences
# dist.coocurrence <- Dist(covered.seq.percentage ,method = 'pearson')
# 
# ## Create the heatmap using 4 agglomeration rules
# th <- sapply(c("average", "complete", "single", "ward"), function(m){
#   
#   if(m == "ward"){
#     temp <- m
#     m <- "ward.D"
#   }
#   
#   ## Calculate the hierarchical tree
#   col.order <- hclust(dist.coocurrence, method = m)[[3]]
#   
#   if(m == "ward.D"){
#     m <- "ward"
#   }
#   
#   ## The col and rows have the same order
#   comp.order.list.rows[[m]] <<- paste(col.order, collapse = ",")
#   comp.order.list.cols[[m]] <<- paste(col.order, collapse = ",")
#   
# })
# comp.order.list.rows[["default"]] <- paste(1:length(matrix.names), collapse = ",")
# comp.order.list.cols[["default"]] <- paste(1:length(matrix.names), collapse = ",")
# 
# ## Convert the coverage table to the format required in D3 heatmap
# coocurrence.table.d3 <- paste(basename(prefix), "_coocurrence_table.tsv", sep = "")
# y <- NULL
# for(j in 1:dim(covered.seq.percentage)[1]){
#   for(i in 1:dim(covered.seq.percentage)[2]){
#     y <<- rbind(y, matrix(c(j,i, as.numeric(covered.seq.percentage[j,i])), nrow = 1))
#   }
# }
# colnames(y) <- c("Row", "Col", "Value")
# y <- as.data.frame(y)
# verbose(paste("Exporting data with co-ocurreence percentage for D3 dynamic heatmap", coverage.table.d3), 2)
# write.table(y, file = coocurrence.table.d3, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# left <- (max(as.vector(sapply(matrix.names, nchar))) + 2) * 10

###################################################################
## Step 21: Add the columns Profile cluster to the feature table ##
###################################################################
list.counter <- 1
features.table$Profile_cluster <- profile.clusters.names[[list.counter]]


#################################################################
## Step 22: Create the dynamic report in HTML                  ##
## This is a general procedure (not only matrix-scan specific) ##
## A HTML template is modified with the current data stored in ##
## the dataframe freq.per.bin complemented with the statistics ##
## calculated in the dataframe feature.attributes              ##
#################################################################
verbose(paste("Creating HTML dynamic report"), 1)
#############
## JS code ##
#############

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


binthrash <- sapply(1:length(bins), function(list.counter){
  
  ## Order the TF.names (and other variables) according the Significance (-log10(E-value))
  order.by.eval <- order(as.numeric(as.vector(features.table$E_val)))
  features.table <- features.table[order.by.eval,]
  
  ## Set the motif names and IDs
  TF.IDs <- matrix.names[order.by.eval]
  
  ## Get the IDs of the TF
  TF.names <- as.vector(sapply(TF.IDs, function(x){
    as.vector(ID.names[which(ID.names[,1] == x),2])
  }))
  
  ## Set colors
  set.colors <- colorRampPalette(brewer.pal(10,"Paired"))(nb.motifs)
  
  ## Initialize variables
  counter <- 0
  x.correspondence <- vector("list", length(bins))
  x.y <- vector("list", length(bins))
  x.y.coverage <- vector("list", length(bins))
  x.y.coverage.names <- vector("list", length(bins))
  plot.names <- vector("list", length(bins))
  plot.names.cover <- vector("list", length(bins))
  area <- vector("list", length(bins))
  all.motifs <- vector("list", length(bins))
  all.motifs.cover <- vector("list", length(bins))
  all.motif.names <- vector("list", length(bins))
  hash.motif.IDs <- list()
  
  
  ## Each column of the variable profiles correspond to the counts per bin of each motif
  verbose(paste("Drawing dynamic profile plot"), 1)
  # thrash <- apply(freq.per.bin[[list.counter]][order.by.eval,], 1, function(values){
    thrash <- apply(freq.per.bin[order.by.eval,], 1, function(values){    
    
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
    
    all.motifs[[list.counter]][counter] <<- motif
    all.motif.names[[list.counter]][counter] <<- TF.names[counter]
    
    ## Create the name's data for the HTML file
    ## The name correspond to the motif name. 
    ## Note that two motifs for the same TF will have the same name and this name will appears twice 
    ## in the report. However their IDs are unique.
    plot.names[[list.counter]] <<- append(plot.names[[list.counter]], paste("'", motif, "' : '",  TF.names[counter],"',", sep = ""))
    
    ## Create the area's data for the HTML file (optional)
    area[[list.counter]] <<- append(area[[list.counter]], paste("'", motif, "' : 'area',", sep = ""))
    
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
    x.y[[list.counter]] <<- rbind(x.y[[list.counter]], y) 
    
    ## Convert the X-Y values to the format required for C3 plot coverage
    xx.sig <- paste("['", motif.cover, "_x',", round(x.sig[[list.counter]][order.by.eval][counter]), "],", sep = "")
    yy.cov <- paste("['", motif.cover, "',", round(y.cov[[list.counter]][order.by.eval][counter], digits = 2), "],", sep = "")
    x.y.coverage[[list.counter]] <<- append(x.y.coverage[[list.counter]], xx.sig)
    x.y.coverage[[list.counter]] <<- append(x.y.coverage[[list.counter]], yy.cov)
    
    ## Add the motifs IDs sentences for the coverage plot
    plot.names.cover[[list.counter]] <<- append(plot.names.cover[[list.counter]], paste("'", motif.cover, "' : '",  TF.names[counter],"',", sep = ""))
    
    ## Append all the motifs IDs for the coverage plot
    all.motifs.cover[[list.counter]] <<- append(all.motifs.cover[[list.counter]], motif.cover)
    
    ## Add the motif names for the coverage plot
    name.cov <- paste("'", motif.cover, "' : '", motif.cover, "_x',", sep = "")
    x.y.coverage.names[[list.counter]][counter] <<- name.cov
    
  })
  all.motifs.cover.hash <- data.frame(all.motifs.cover[[list.counter]], all.motif.names[[list.counter]])
  
  ## Set the line width according the significance -log10(E-value)
  ## Higher significance means a wider line
  significance <- as.numeric(as.vector(features.table$Sig))
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
  
  all.motifs.function <- paste(paste("'", all.motifs[[list.counter]], "'", sep = ""), collapse = ",")
  all.motifs.function.cov <- paste(paste("'", all.motifs.cover[[list.counter]], "'", sep = ""), collapse = ",")
  
  ## Generate the JS function to show the clusters
  thrash <- sapply(1:length(cluster.profiles.motif.names[[list.counter]]), function(cl){
    
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
    cluster.member.names <- as.vector(unlist(sapply(as.vector(cluster.profiles.motif.names[[list.counter]][[cl]]), function(m){
        as.vector(hash.motif.IDs[[m]])
    })))
    cluster.member.names <- paste(paste("'", cluster.member.names, "'", sep = ""), collapse = ",")
    
    ## Cluster member names (cover plot)
    cluster.member.names.cov <- as.vector(unlist(sapply(as.vector(cluster.profiles.motif.names[[list.counter]][[cl]]), function(m){
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
    paste(logo.folder, "/", i, "_logo.jpeg", sep = "")
  })
  
  logos.R <- sapply(TF.IDs, function(i){
    paste(logo.folder, "/", i, "_logo_rc.jpeg", sep = "")
  })
  
  ## Write the Profile and TFBSs plots path
  profiles.plots <- sapply(TF.IDs, function(i) {
    paste(basename(prefix), "_TFBSs_positional_profiles/", i, "_positional_profile_bins_", bins[list.counter],".jpeg", sep = "")
  })
  
  tfbss.plots <- sapply(TF.IDs, function(i) {
    paste(basename(prefix), "_TFBSs_pval_distribution/", i, "_TFBSs_pval_classes.jpeg", sep = "")
  })
  
  tfbs.weight.boxplot <- sapply(TF.IDs, function(i) {
    paste(basename(prefix), "_Boxplot_Weight/", i, "_Weight_distribution_per_bin.jpeg", sep = "")
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
  all.pval.match <- rep(p.val, times = nb.motifs)
  datatable.info.tab <- features.table
  datatable.info.tab$P_val_threshold <- all.pval.match
  datatable.info.tab$ID <- TF.IDs
  datatable.info.tab$Names <- TF.names
  datatable.info.tab$Profiles <- profiles.plots
  datatable.info.tab$TFBS <- tfbss.plots
  datatable.info.tab$TFBS_per_seq <- tfbss.per.seq.plots
  datatable.info.tab$Logo <- logos.F
  datatable.info.tab$Logo_RC <- logos.R
  datatable.info.tab$covered_files <- covered.files
  datatable.info.tab$not_covered_files <- not.covered.files
  datatable.info.tab$w_boxplot <- tfbs.weight.boxplot
  all.motifs <- all.motifs[[list.counter]]
  all.motif.names <- all.motif.names[[list.counter]]
  
  #############################################################
  ## Fill the HTML template                                  ##
  ## Substitute the words marked in the template by the data ##
  #############################################################
  html.report <- readLines(html.template.file)
  verbose(paste("Fill html report"), 1)
  

  # [1] "Feature"           "Chi_squared"       "Degrees"          
  # [4] "Chi_app"           "Sig"               "P_val"            
  # [7] "E_val"             "Q_val"             "Coverage"         
  # [10] "Sequences"         "Nb_hits"           "Pval_scan"        
  # [13] "Covered_seq"       "Not_covered_seq"   "Profile_cluster"  
  # [16] "P_val_threshold"   "ID"                "Names"            
  # [19] "Profiles"          "TFBS"              "TFBS_per_seq"     
  # [22] "Logo"              "Logo_RC"           "covered_files"    
  # [25] "not_covered_files" "w_boxplot"        

  
  # Motif_ID	Motif_name	P-value	E-value	Significance	FDR	Nb of hits	Nb of sequences	Fraction of sequences	Chi-squared	Chi applicability	Profile cluster	Profile	TFBSs	TFBSs per seq	Logo	Logo (RC)	Covered sequences	Not Covered sequences
  
  profile.data.tab.html <- create.html.tab(datatable.info.tab[,c(1,18,6,7,5,8,11,10,9,12,2,4,15,19,20,21,22,23,24,25)], img = c(17,18), plot = c(14,15,16), link.text.covered = 19, link.text.not.covered = 20)
  profile.data.tab.html <- gsub("Inf", "&infin;", profile.data.tab.html)
  
  profile.data.tab.html <- paste(profile.data.tab.html, collapse = "\n")
  html.report <- gsub("--tab--", profile.data.tab.html, html.report)
  
  ## Define the x-axis categories
  x.axis.categories <- paste(paste("'", colnames(freq.per.bin), "'", sep = ""), collapse = ",")
  
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
    paste(logo.folder, "/", i, "_logo.jpeg", sep = "")
  })
  logos <- paste("pics['", all.motifs, "'] = '", as.vector(datatable.info.tab$Logo), "';", sep = "")
  logos <- paste(logos, collapse = "\n")
  html.report <- gsub("--pics--", logos, html.report)
  
  logos.co <- paste("pics_m1['", all.motifs, "'] = '", as.vector(datatable.info.tab$Logo), "';", sep = "")
  logos.co <- paste(logos.co, collapse = "\n")
  html.report <- gsub("--pics_m1--", logos.co, html.report)
  
  ## Logos in Reverse complement
  logos.rc <- sapply(TF.IDs, function(i){
    paste(logo.folder, "/", i, "_logo_rc.jpeg", sep = "")
  })
  logos.rc <- paste("pics_rc['", all.motifs, "'] = '", as.vector(datatable.info.tab$Logo_RC), "';", sep = "")
  logos.rc <- paste(logos.rc, collapse = "\n")
  html.report <- gsub("--pics_rc--", logos.rc, html.report)
  
  ## Add the signficance (to display in the tooltip)
  ## They are inserted in the JS section
  sig <- paste("significances['", all.motifs, "'] = ", as.vector(datatable.info.tab$Sig), ";", sep = "")
  sig <- paste(sig, collapse = "\n")
  html.report <- gsub("--significances--", sig, html.report)
  
  ## Add the coverages (to display in the tooltip)
  ## They are inserted in the JS section
  cc <- as.numeric(gsub("%", "", features.table$Coverage))
  coverage <- paste("TF_coverage['", all.motifs, "'] = ", as.vector(cc), ";", sep = "")
  coverage <- paste(coverage, collapse = "\n")
  html.report <- gsub("--TF_covertures--", coverage, html.report)
  
  ## The plot heigth depends in the number of motifs
  motif.total <- length(all.motifs)
  chart.heigth <- 500
  if(motif.total >= 200){
    chart.heigth <- 850
  } else if(motif.total >= 300){
    chart.heigth <- 1200
  } else if(motif.total >= 600){
    chart.heigth <- 1600
  }
  html.report <- gsub("--chart_h--", chart.heigth, html.report)
  
  ## Add x values (one row per motif)
  ## They are inserted in the C3 section
  xx <- paste(x.y[[list.counter]], collapse = "\n")
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
  plot.names <- paste(plot.names[[list.counter]], collapse = "\n")
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
  max.y <- max(freq.per.bin) + 0.02
  html.report <- gsub("--y_axis--", max.y, html.report)
  
  ## Fill the parameters table
  verbose(paste("Creating parameters table"), 1)
  html.report <- gsub("--bin_l--", bins[list.counter], html.report)
  html.report <- gsub("--bin_nb--", ncol(freq.per.bin), html.report)
  html.report <- gsub("--seq_nb--", total.scanned.sequences, html.report)     ## Don't forget length(seq.id)
  html.report <- gsub("--motif_nb--", nb.motifs, html.report)
  html.report <- gsub("--p--", prettyNum(p.val), html.report)
  
  ## Fill the heatmap section
  verbose(paste("Creating Heatmap section"), 1)
  html.report <- gsub("--heatmap_png--", paste(basename(prefix), "_profiles_heatmap.jpeg", sep = ""), html.report)
  html.report <- gsub("--heatmap_pdf--", paste(basename(prefix), "_profiles_heatmap.pdf", sep = ""), html.report)
  
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
  verbose(paste("Drawing Significance vs Coverage plot"), 1)
  x.y.coverage <- paste(x.y.coverage[[list.counter]], collapse = "\n")
  html.report <- gsub("--x_y_coverture--", x.y.coverage, html.report)
  
  ## Insert the X-Y scatterplot xs
  x.y.coverage.names <- paste(x.y.coverage.names[[list.counter]], collapse = "\n")
  html.report <- gsub("--xs_coverture--", x.y.coverage.names, html.report)
  
  ## Insert the names in the coverage XY-plot
  plot.names.cover <- paste(plot.names.cover[[list.counter]], collapse = "\n")
  html.report <- gsub("--names_cov--", plot.names.cover, html.report)
  
  ## Add the real motif IDs (to display in the tooltip)
  ## They are inserted in the JS section
  IDs.cov <- paste("cov_IDs['", all.motifs.cover[[list.counter]], "'] = '", TF.IDs, "';", sep = "")
  IDs.cov <- paste(IDs.cov, collapse = "\n")
  html.report <- gsub("--IDs_cov--", IDs.cov, html.report)
  
  ## Add the profile cluster (to display in the tooltip)
  ## They are inserted in the JS section
  cov.profile.clusters <- paste("cov_profile_clusters['", all.motifs.cover[[list.counter]], "'] = '", as.vector(datatable.info.tab$Profile_cluster), "';", sep = "")
  cov.profile.clusters <- paste(cov.profile.clusters, collapse = "\n")
  html.report <- gsub("--profile_clusters_array_cov--", cov.profile.clusters, html.report)
  
  ## Add the real motif logo path (to display in the tooltip)
  ## They are inserted in the JS section
  logos.cov <- sapply(TF.IDs, function(i){
    paste(logo.folder, i, "_logo.jpeg", sep = "")
  })
  logos.cov <- paste("cov_pics['", all.motifs.cover[[list.counter]], "'] = '", as.vector(datatable.info.tab$Logo), "';", sep = "")
  logos.cov <- paste(logos.cov, collapse = "\n")
  html.report <- gsub("--pics_cov--", logos.cov, html.report)
  
  ## Logos in Reverse complement
  logos.rc.cov <- sapply(TF.IDs, function(i){
    paste(logo.folder, i, "_logo_rc.jpeg", sep = "")
  })
  logos.rc.cov <- paste("cov_pics_rc['", all.motifs.cover[[list.counter]], "'] = '", as.vector(datatable.info.tab$Logo_RC), "';", sep = "")
  logos.rc.cov <- paste(logos.rc.cov, collapse = "\n")
  html.report <- gsub("--pics_rc_cov--", logos.rc.cov, html.report)
  
  ## Add the signficance (to display in the tooltip)
  ## They are inserted in the JS section
  ss <- as.numeric(gsub("%", "", features.table$Sig))
  ss[ss == Inf] <- 350
  sig.cov <- paste("cov_significances['", all.motifs.cover[[list.counter]], "'] = ", as.vector(ss), ";", sep = "")
  sig.cov <- paste(sig.cov, collapse = "\n")
  html.report <- gsub("--significances_cov--", sig.cov, html.report)
  
  ## Add the coverages (to display in the tooltip)
  ## They are inserted in the JS section
  cc <- as.numeric(gsub("%", "", features.table$Coverage))
  coverage.cov <- paste("cov_TF_coverage['", all.motifs.cover[[list.counter]], "'] = ", as.vector(cc), ";", sep = "")
  coverage.cov <- paste(coverage.cov, collapse = "\n")
  html.report <- gsub("--TF_covertures_cov--", coverage.cov, html.report)
  
  ## Add the coverages (to display in the tooltip)
  ## They are inserted in the JS section
  all.profiles.pics <- paste("'", as.vector(datatable.info.tab$Profiles), "'", sep = "")
  profiles.pics.cov <- paste("cov_pics_profile['", all.motifs.cover[[list.counter]], "'] = ", all.profiles.pics, ";", sep = "")
  profiles.pics.cov <- paste(profiles.pics.cov, collapse = "\n")
  html.report <- gsub("--profile_pics_cov--", profiles.pics.cov, html.report)
  
  all.profiles.pics.co <- paste("'", as.vector(datatable.info.tab$Profiles), "'", sep = "")
  profiles.pics.cov.co <- paste("profiles_m1['", all.motifs, "'] = ", all.profiles.pics.co, ";", sep = "")
  profiles.pics.cov.co <- paste(profiles.pics.cov.co, collapse = "\n")
  html.report <- gsub("--profiles_m1--", profiles.pics.cov.co, html.report)
  
  ## Insert the motif names (to hide/show all) in coverage plot
  ## They are inserted in the JQuery section
  all.motifs.cover <- paste(paste("'", all.motifs.cover[[list.counter]], "'", sep = ""), collapse = ",")
  html.report <- gsub("--all_cover--", all.motifs.cover, html.report)
  
  ###################################################
  ## Fill the data for the co-ocurrence d3 heatmap
  # html.report <- gsub("--default_r_number--", comp.order.list.rows[["default"]], html.report)
  # html.report <- gsub("--default_c_number--", comp.order.list.cols[["default"]], html.report)
  # 
  # html.report <- gsub("--average_r_number--", comp.order.list.rows[["average"]], html.report)
  # html.report <- gsub("--average_c_number--", comp.order.list.cols[["average"]], html.report)
  # 
  # html.report <- gsub("--complete_r_number--", comp.order.list.rows[["complete"]], html.report)
  # html.report <- gsub("--complete_c_number--", comp.order.list.cols[["complete"]], html.report)
  # 
  # html.report <- gsub("--single_r_number--", comp.order.list.rows[["single"]], html.report)
  # html.report <- gsub("--single_c_number--", comp.order.list.cols[["single"]], html.report)
  # 
  # html.report <- gsub("--ward_r_number--", comp.order.list.rows[["ward"]], html.report)
  # html.report <- gsub("--ward_c_number--", comp.order.list.cols[["ward"]], html.report)
  # 
  # html.report <- gsub("--c_numb--", length(matrix.names), html.report)
  # html.report <- gsub("--r_numb--", length(matrix.names), html.report)
  # 
  # cell.size <- 20
  # html.report <- gsub("--cell_size--", cell.size , html.report)
  # 
  # html.report <- gsub("--left--", left, html.report)
  # 
  # coocurrence.table.d3 <- paste(basename(prefix), "_coocurrence_table.tsv", sep = "")
  # html.report <- gsub("--file--", coocurrence.table.d3, html.report)
  # 
  # domain <- paste(domain, collapse=",")
  # html.report <- gsub("--domain--", domain, html.report)
  # 
  # coocurrence.palette <- c("#FFFFFF", coocurrence.palette) 
  # gradient <- paste("[", paste(paste("'", coocurrence.palette, "'", sep=""), collapse=","), "];", sep = "")
  # html.report <- gsub("--gradient--", gradient, html.report)
  # 
  # legend <- seq(from = 0, to = 100, by = 5)
  # legend <- paste(legend, collapse=",")
  # html.report <- gsub("--data_legend--", legend, html.report)
  # 
  # col.labels <- paste(paste("'", matrix.names, "'", sep = ""), collapse = ",")
  # row.labels <- paste(paste("'", matrix.names, "'", sep = ""), collapse = ",")
  # html.report <- gsub("--col_label--", col.labels, html.report)
  # html.report <- gsub("--row_label--", row.labels, html.report)
  # 
  # html.body.size <- 200 + left + (length(matrix.names)*20) + 30
  # html.report <- gsub("--body--", html.body.size, html.report)
  
  html.report <- gsub("--function_profile_clusters--", all.cluster.functions, html.report)
  html.report <- gsub("--button_clusters--", all.cluster.buttons, html.report)
  html.report <- gsub("--function_profile_clusters_points--", all.cluster.functions.cov, html.report)
  html.report <- gsub("--cluster_points--", all.cluster.buttons.cov, html.report)
  
  ## Insert the motif names (to hide/show all)
  ## They are inserted in the JQuery section
  all.motifs <- paste(paste("'", all.motifs, "'", sep = ""), collapse = ",")
  html.report <- gsub("--all--", all.motifs, html.report)
  
  ## Export the report
  verbose(paste("Exporting html report"), 1)
  html.report.file <- paste(basename(prefix), "_scan_profile_report.html", sep = "")
  write(html.report, file = html.report.file)
})
rm(binthrash)
# grep -v '^;' matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments.tab | awk -F '\t'  ' $2!="limit" && ($11 >= 4) {print $1"\t"$3"\t"($6+$5)/2"\t"$9} ' > matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments_PARSED.tab
# cat /home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/enrichment_by_scan/Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " matrix.scan.active = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; matrix.scan.inactive = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; p.val = '1e-4'; bin = '50'; pdf.file = './test.pdf'"


########################
## Centrimo algorithm ##
########################  MA0024_3

# ## Select the 'best' site on each sequence
# matrix.query.selection.best.score <- ddply(matrix.query.selection, "seq_id", mutate, best = max(Pval.minlog10, na.rm = TRUE))
# matrix.query.selection.best.score <- subset(matrix.query.selection.best.score , Pval.minlog10 == best, select = c("seq_id", "bspos", "Pval.minlog10"))
# 
# ## Randomize the rows
# matrix.query.selection.best.score <- matrix.query.selection.best.score[sample(1:nrow(matrix.query.selection.best.score)),]
# 
# ## Take the not-duplicated sequences (for this the rows are firstly randomized)
# matrix.query.selection.best.score <- matrix.query.selection.best.score[!duplicated(matrix.query.selection.best.score$seq_id), ]
# 
# ## Count the number of hits per bin
# selection.IR.best <- IRanges(start = matrix.query.selection.best.score$bspos, end = matrix.query.selection.best.score$bspos)    ## Count the overlap of BS in the bins
# counts.per.bin.best <- countOverlaps(windows, selection.IR.best)
# 
# return(list(counts = counts.per.bin, counts.best = counts.per.bin.best))
## Total number of best hits per bin
# best.counts.per.bin.table <- as.data.frame(matrix(unlist(counts.per.bin[2,]), nrow = length(counts.per.bin[2,])))
# rownames(best.counts.per.bin.table) <- matrix.names
# colnames(best.counts.per.bin.table) <- as.character(xlab)
