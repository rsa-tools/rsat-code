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

# matrix.scan.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/PHO_summit_coordinates_best_score_pm300_position_scan_top500peaks_mkv_1_pval_1e-3_bin_size_25_all_discovered_motifs_matrix_scan_results_PARSED.tab"
# prefix <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/PHO_summit_coordinates_best_score_pm300_position_scan_top500peaks_mkv_1_pval_1e-3_bin_size_25_all_discovered_motifs"
# ID.to.names.correspondence.tab <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/PHO_summit_coordinates_best_score_pm300_position_scan_top500peaks_mkv_1_pval_1e-3_bin_size_25_all_discovered_motifs_TF_ID_name_correspondence.tab"
# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo")
# sequence.names.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/PHO_summit_coordinates_best_score_pm300_position_scan_top500peaks_mkv_1_pval_1e-3_bin_size_25_all_discovered_motifs_matrix_scan_sequence_names.tab"

# matrix.scan.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/Epromoters/HELA_bin_size_50_pval1e-3_matrix_scan_results_PARSED.tab"
# prefix <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/Epromoters/HELA_bin_size_50_pval1e-3"
# ID.to.names.correspondence.tab <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/Epromoters/HELA_bin_size_50_pval1e-3_TF_ID_name_correspondence.tab"
# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/Epromoters")
# sequence.names.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Demo/Epromoters/HELA_bin_size_50_pval1e-3_matrix_scan_sequence_names.tab"
# seq.length <- 2000
# bin <- 50

####################################
## Step 3: Read matrix-scan table
verbose(paste("Reading matrix-scan results table"), 1)
matrix.scan.results <- read.csv(file = matrix.scan.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results) <- c("seq_id", "ft_name", "bspos", "Pval")

#######################################
## Step 4: Read sequence names table
verbose(paste("Reading sequence names table"), 1)
sequence.names.tab <- read.csv(file = sequence.names.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(sequence.names.tab) <- c("seq_id")
total.scanned.sequences <- length(as.vector(sequence.names.tab$seq_id))
scanned.sequences <- unique(as.vector(sequence.names.tab$seq_id))

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

#################
## Set p-value
p.val <- as.numeric(p.val)

######################################
## Get the sequences + motif name's
seq.id <- unique(as.vector(matrix.scan.results$seq_id))
matrix.names <- unique(as.vector(matrix.scan.results$ft_name))
nb.motifs <- length(matrix.names)

###################################
## Calculate the sequence limits
seq.length <- as.numeric(seq.length)
limits <- seq.length/2

################################
## Step 6: Load the motif IDs
ID.names.tab <- ID.to.names.correspondence.tab
ID.names <- read.table(ID.names.tab, sep = "\t")

##################################################################
## Step 7: Plot the distribution of TFBSs at different p-values ##
##################################################################
setwd(results.folder)
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

thr <- sapply(1:nb.motifs, function(m){
  
  ## Get the matrix name
  matrix.query <- matrix.names[m]
  matrix.query.name <- as.vector(ID.names[,2][which(ID.names[,1] == matrix.query)][1])
  
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
  
  ## Generate GGplot for number of TFBSs per sequence
  TFBSs.per.seq.file <- paste(basename(prefix), "_TFBSs_per_seq/", matrix.query, "_TFBSs_per_seq", sep = "")
  
  aa <- data.frame(nb.seq = nb.hits.per.sequence, nb.hits = as.numeric(names(nb.hits.per.sequence)))
  ggplot(data=aa, aes(y = nb.seq, x=nb.hits)) +
    geom_bar(aes(fill=nb.seq), stat="identity", position=position_dodge()) +
    labs(title="Number of hits", y = "Number of sequences", x="Number of TFBSs") +
    geom_text(aes(label=nb.seq), vjust=-0.15, color="black", size=3)
  
  suppressMessages(ggsave(paste(TFBSs.per.seq.file, ".jpeg", sep = "")))
  suppressMessages(ggsave(paste(TFBSs.per.seq.file, ".pdf", sep = "")))
  
  ## Get the number of putative TFBSs and the number of sequences with 
  ## at least one match of the query matrix
  nb.TFBSs <- dim(matrix.query.selection)[1]
  nb.seq <- length(as.vector(unique(matrix.query.selection$seq_id)))
  
  TFBSs.pval.distribution.file <- paste(basename(prefix), "_TFBSs_pval_distribution/", matrix.query, "_TFBSs_pval_classes", sep = "")
  
  ################################
  #Create a custom color scale
  myColors <- colorRampPalette(brewer.pal(nb.color.classes, "YlGnBu"), space="Lab")(length(classes.pval.letters))
  names(myColors) <- classes.pval.letters
  
  ## Range of p-values for the query motif
  pval.class.matrix.query <- sort(unique(matrix.query.selection$Pval.class))
  
  ## Insert logo
  logo.file <- paste(logo.folder, "/", matrix.query, "_logo.png", sep = "")
  logo <- readPNG(logo.file)
  logo.roster <- rasterGrob(logo, interpolate = TRUE)
  
  ## X position of plot annotations
  text.xmax <- min(matrix.query.selection$bspos) + max(matrix.query.selection$bspos) / 4
  text.center <- (min(matrix.query.selection$bspos) - text.xmax)*2
  
  ggplot(matrix.query.selection, aes(x=bspos, y=Pval.minlog10)) +
    ylim(c(min(matrix.scan.results$Pval.minlog10), max(matrix.scan.results$Pval.minlog10))) +
    geom_point(aes(colour = Pval.class.letter), shape = "O", size = 5, stroke = 3) +
    geom_rug(position='jitter') +
    labs(title=paste("Qualitative distribution of ", matrix.query, " TFBSs", sep = ""), y = "-log10(P-value)", x = "Position") +
    scale_colour_manual(name = "-log10(P-value)",values = myColors, labels = paste(">", pval.class.matrix.query, sep = "")) +
    theme_minimal() +
    annotate("text", x = -limits + ((limits*2)/10), y = max.pval.minus.log10 - 0.25, label = paste("Nb of TFBSs: ", nb.TFBSs, sep = ""), size = 4, hjust = 0) + 
    annotate("text", x = -limits + ((limits*2)/10), y = max.pval.minus.log10 - 0.55, label = paste("Nb of sequences: ", nb.seq, sep = ""), size = 4, hjust = 0) + 
    annotation_custom(logo.roster, xmax = limits - (limits/3), xmin = limits - 5, ymin = max.pval.minus.log10 - 1, ymax = max.pval.minus.log10 - 0.05)
  
  suppressMessages(ggsave(paste(TFBSs.pval.distribution.file, ".pdf", sep = "")))
  suppressMessages(ggsave(paste(TFBSs.pval.distribution.file, ".jpeg", sep = "")))
})
rm(thr)

#################################################################################
## Step 8: Calculate the raw counts of TFBSs on each position of the sequences ##
#################################################################################

verbose(paste("Calculating the raw counts of TFBSs in the sequences"),1)
matrix.scan.results$bspos.left <- matrix.scan.results$bspos + seq.length
raw.counts.all.motifs <- sapply(matrix.names, function(m){
  
  sites.m <- as.vector(matrix.scan.results[which(matrix.scan.results$ft_name == m),"bspos.left"])
  raw.counts.m <- tabulate(sites.m, nbins = seq.length)
  return(raw.counts.m)
})
raw.counts.all.motifs <- as.data.frame(raw.counts.all.motifs)

#######################################################
## Step 9: Calculate the raw counts of TFBSs per bin ##
#######################################################
## Define the bins size to be tested
## Min:5 - Max:100
## All the bins that will be a divisor of the sequences length between this range will be considered
bins <- seq(5, seq.length/2, 1)
bins <- bins[(seq.length/2) %% bins == 0]
bins <- as.numeric(bin)

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
counts.per.bin <- vector("list", length(bins))
for(b in 1:length(bins)){
  verbose(paste("Calculating the raw counts of TFBSs in the sequences - Bin of size:", bins[b]),1)
  counts.per.bin.query <- rollapply(raw.counts.all.motifs, bins[b], sum, na.rm = TRUE, by = bins[b], partial = TRUE, align="left")
  colnames(counts.per.bin.query) <- matrix.names
  rownames(counts.per.bin.query) <- xlab[[b]]
  counts.per.bin[[b]] <- counts.per.bin.query
}

#############################################################
## Step 10: Calculate number of hits and matched sequences ##
## Export files with matched and unmatched sequences       ##
#############################################################

## Count the number of matched sequences
verbose(paste("Counting the number of hits per motif"),1)
matched.sequences <- sapply(matrix.names, function(m){
  matrix.scan.results.subset <- subset(matrix.scan.results, ft_name == m)
  length(unique(as.vector(matrix.scan.results.subset$seq_id)))
})

## Count matches for each motifs
verbose(paste("Counting the number of matched sequences per motif"),1)
matches.per.motif <- ddply(matrix.scan.results, "ft_name", summarize, x = length(bspos))
matches.per.motif <- matches.per.motif$x
names(matches.per.motif) <- matrix.names

## Matched and unmatched sequences
verbose(paste("Exporting names of matched and not-matched sequences per motif"),1)
covered.seq <- rep(0, nb.motifs)
not.covered.seq <- rep(0, nb.motifs)
names(covered.seq) <- matrix.names
names(not.covered.seq) <- matrix.names
list.counter <- 0
thrash <- sapply(matrix.names, function(m){
  
  list.counter <<- list.counter + 1
  
  matrix.scan.results.subset <- subset(matrix.scan.results, ft_name == m)
  
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
## Step 11: Calculate the Frequency tables ##
#############################################
dir.create(paste(basename(prefix), "_tables/", sep = ""), recursive = TRUE, showWarnings = FALSE )
freq.per.bin <- NULL
list.counter <- 0
counts.tab.file <- NULL
density.tab.file <- NULL
freq.per.bin <- lapply(counts.per.bin, function(cpb){
  
  list.counter <<- list.counter + 1
  
  ## Count the total of TFBSs
  total.hits <- colSums(cpb)
  ## Divide the sum at each bin by the total
  freq.per.bin.query <- as.data.frame(t(cpb)/total.hits)
  
  ##########################################
  ## Export Counts and Frequencies tables
  density.tab.file <<- paste(basename(prefix), "_tables/density_per_bin_profiles_bin", bins[list.counter], ".tab", sep = "")
  write.table(freq.per.bin.query, file = density.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  
  counts.tab.file <<- paste(basename(prefix), "_tables/counts_per_bin_profiles_bin", bins[list.counter], ".tab", sep = "")
  write.table(cpb, file = counts.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  
  return(freq.per.bin.query)
})
## Calculate the highest frequency of TFBSs per each list
max.y <- lapply(freq.per.bin, max, na.rm = TRUE)

##############################################
## Step 12: Chi-squared calculation section ##
##############################################
counts.per.bin <- lapply(counts.per.bin, t)
list.counter <- 0
all.chi.results <- NULL
thrash <- lapply(counts.per.bin, function(cpb){
  
  list.counter <<- list.counter + 1
  
  ## Calculate the X2 test
  ## H0 = Homogeneous distribution of 
  chi.vector <- apply(cpb, 1, chisq.test, correct = TRUE)
  
  chi.fields.matrix <- sapply(chi.vector, function(x){
    selected.fields <- as.vector(c(x[[1]], x[[2]], x[[3]], bins[list.counter]))
    names(selected.fields) <- NULL
    return(selected.fields)
  })
  chi.fields.matrix <- as.data.frame(t(chi.fields.matrix))
  colnames(chi.fields.matrix) <- c("Chi", "DF", "Pvalue", "Bin_size")
  
  ## Concatenate all results
  all.chi.results[[list.counter]] <<- rbind(all.chi.results, chi.fields.matrix)
})
rm(thrash)

########################################################################
## Step 13: calculate Q-value, E-value, and significance -log(pvalue) ##
## Re-order the table                                                 ##
########################################################################
features.table <- lapply(all.chi.results, function(df){
  df$Chi <- round(df$Chi, digits = 3)
  
  ## Calculate E-value
  df$Evalue <- df$Pvalue * nb.motifs
  
  ## Calculate Significance
  df$Sig <- round(-log10(df$Evalue), digits = 3)
  
  ## Calculate q-values
  ## This step is executed once all the p-values were calculated
  ## The variable with class 'qvalue' is stored to its further exportation
  pp <- as.numeric(as.vector(df$Pvalue))
  features.qvalues <- p.adjust(pp, method = "BH")
  df$Qvalue <- prettyNum(features.qvalues, scientific=TRUE, digits = 2)
  
  ## P-val -> Pretty number
  df$Pvalue <- round(df$Pvalue, digits = 100000000000000000)
  df$Pvalue <- prettyNum(df$Pvalue, scientific=TRUE, digits = 2)
  
  ## E-val -> Pretty number
  df$Evalue <- prettyNum(df$Evalue, scientific=TRUE, digits = 2)
  
  ## Add the matrix Id as feature ID
  df$feature <- rownames(df)
  
  ## Add the coverage column
  df$coverage <- round(matched.sequences/total.scanned.sequences, digits = 2)
  
  df$sequences <- matched.sequences
  
  ## Add the number of matches per motif
  df$hits <- matches.per.motif
  
  ## Aqui
  df$cov <- covered.seq
  df$not_cov <- not.covered.seq
  
  df.sub <- df[, c(8,1,2,6,3,5,7,9,10,11,12,13)]
  colnames(df.sub) <- c("Feature", "Chi_squared", "Degrees", "Sig", "P_val", "E_val", "Q_val", "Coverage", "Sequences", "Nb_hits", "Covered_seq", "Not_covered_seq")
  
  return(df.sub)
})
rm(all.chi.results)

feature.attributes.file <- paste(basename, "_attributes.tab", sep = "")
write.table(features.table, file = feature.attributes.file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


###################################
## Step 14: Calculate BG matches ##
###################################
verbose(paste("Calculating the BG number of matches per motif"),1)

list.counter <- 0
counts.norm <- list()
thrash <- sapply(bins, function(b){
  
  list.counter <<- list.counter + 1
  
  ## Mean of TFBSs per bin
  mean.counts.all.motifs.bins <- rollapply(raw.counts.all.motifs, b, mean, na.rm = TRUE, by = b, partial = TRUE, align="left")
  
  ## Sum the means per bin
  all.columns.trimmed.means <- colSums(mean.counts.all.motifs.bins, na.rm = TRUE)
  
  ## The median of the means will be considered as the "basal" number of hits
  for (m in 1:ncol(mean.counts.all.motifs.bins)){
    # row.mean <- mean(mean.counts.all.motifs.bins[,m],trim=0.3)
    row.mean <- median(mean.counts.all.motifs.bins[,m])
    all.columns.trimmed.means[m] <- row.mean
  }
  
  df_by_bins_norm = t(t(mean.counts.all.motifs.bins)-all.columns.trimmed.means)
  
  counts.norm[[list.counter]] <<- df_by_bins_norm
})

######################################
## Step 15: Define profile clusters ##
######################################
list.counter <- 0
cluster.tree.profiles.palette <- vector("list", length(counts.norm))
color.clusters.tree.profiles <- vector("list", length(counts.norm))
profile.clusters.names <- vector("list", length(counts.norm))
tree.profiles <- vector("list", length(counts.norm))
cluster.profiles.motifs <- NULL
cluster.profiles.motif.names <- NULL
thr <- lapply(counts.norm, function(cn){
  
  list.counter <<- list.counter + 1
  
  tree.profiles[[list.counter]] <<- hclust(Dist(t(cn), method = "correlation"), method = "complete")
  clusters.tree.profiles <- cutreeDynamic(tree.profiles[[list.counter]], minClusterSize = 1, method = "tree")
  names(clusters.tree.profiles) <- tree.profiles[[list.counter]][[4]]
  
  ## Generate a color palette
  nb.profile.clusters <- length(unique(clusters.tree.profiles))
  cluster.tree.profiles.palette[[list.counter]] <<- colorRampPalette(brewer.pal(9, "Set1"), space="Lab")(nb.profile.clusters)
  
  ## Assign a different color to each cluster
  color.clusters.tree.profiles[[list.counter]] <<- as.vector(sapply(clusters.tree.profiles, function(color){
    cluster.tree.profiles.palette[[list.counter]][color]
  }))
  
  ############################################
  ## Fill the Profile_cluster attribute 
  ## Add a new column to all.chi.results df
  profile.clusters.names[[list.counter]] <<- as.vector(sapply(matrix.names, function(m){
    profile.cluster <- as.vector(clusters.tree.profiles[m])
    paste("Profile_cluster_", profile.cluster, sep = "")
  }))
  
  ## Get the member motif IDs of each cluster
  cluster.profiles.counter <- 0
  thrash <- sapply(1:nb.profile.clusters, function(cl){
    
    cluster.profiles.counter <<- cluster.profiles.counter + 1
    cluster.profiles.motifs[[list.counter]][[cluster.profiles.counter]] <<- names(which(clusters.tree.profiles == cluster.profiles.counter))
    
    ## Get the motif name
    cluster.profiles.motif.names[[list.counter]][[cluster.profiles.counter]] <<- as.vector(
      sapply(cluster.profiles.motifs[[list.counter]][[cluster.profiles.counter]], function(n){
        ID.names[which(ID.names[,2] == n),1]
      })
    )
  })
  rm(thrash)
})
rm(thr)


##########################################################
## Step 16: Draw Profiles heatmap                       ##
## Shows the frequencies of hits per bin for each motif ##
##########################################################
verbose(paste("Drawing Heatmap profiles"),1)
list.counter <- 0
thrash <- lapply(counts.norm, function(cn){
  
  list.counter <<- list.counter + 1
  
  ## Profile Heatmap Color palette (user-defined)
  rgb.palette <- rev(colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")(heatmap.color.classes))
  
  ## Remove noisy extremities
  df.limits <- round(range(cn))
  df.limits.abs <- abs(df.limits)
  
  if (df.limits.abs[1] > df.limits.abs[2]){
    cn[cn < (-df.limits[2])] <- (-df.limits[2])
    cn[cn > df.limits[2]] <- df.limits[2]
  } else {
    cn[cn > (-df.limits[1])] <- (-df.limits[1])
    cn[cn < df.limits[1]] <- df.limits[1]
  }
  
  cn.t <- t(cn)
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
  x.sig[[list.counter]] <<- as.numeric(as.vector(features.table[[list.counter]]$Sig))
  inf.val <- which(x.sig[[list.counter]] == Inf)
  x.sig[[list.counter]][inf.val] <<- 350
  names(x.sig[[list.counter]]) <- as.vector(features.table[[list.counter]]$Feature)
  y.cov[[list.counter]] <<- as.numeric(gsub("%", "", features.table[[list.counter]]$Coverage)) *100
  names(y.cov[[list.counter]]) <- as.vector(features.table[[list.counter]]$Feature)
  
  ## Print the Significance vs Coverage plot
  # sig.coverage.file <- paste(basename, "_significance_vs_coverage_bins_", b, ".pdf", sep = "")
  # pdf(sig.coverage.file)
  
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
  selected.TFBMs <- features.table[[list.counter]][which(as.vector(features.table[[list.counter]]$Sig) >= 20 & as.numeric(gsub("%", "", features.table[[list.counter]]$coverage)) >= 66), "Feature"]
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
thr <- sapply(bins, function(b){
  
  list.counter <<- list.counter + 1
  
  thrash <- sapply(matrix.names, function(feature.query){
    
    ## Set Output file name
    file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", feature.query, "_positional_profile_bins_", b, sep = "")
    
    ## Get the coordinates
    y.val <- as.numeric(freq.per.bin[[list.counter]][feature.query,])
    x.val <- as.numeric(xlab[[list.counter]])
    
    ## Load the logo
    logo.file <- paste(logo.folder, "/", feature.query, "_logo.png", sep = "")
    logo <- readPNG(logo.file)
    logo.roster <- rasterGrob(logo, interpolate = TRUE)
    
    ## Plot the profile (using ggplot2)
    xy.df <- data.frame(x = x.val, y = y.val)
    ggplot(xy.df, aes(x=x, y=y)) +
      geom_line(colour = "#00BFC4", size = 3) +
      ylim(0, max(max.y[[list.counter]])) +
      labs(title=paste(feature.query, " binding profile", sep = ""), y = "Frequency of TFBSs", x = "Position") +
      geom_rug(position='jitter', sides="l") +
      geom_area(fill = "#00BFC4", alpha=0.35) + 
      annotation_custom(logo.roster, xmax = limits, xmin = limits - sum(abs(range(xy.df$x)))/5, ymin = max.y[[list.counter]] - 0.01, ymax = max.y[[list.counter]] - 0.075)
    
    ## Export the file
    suppressMessages(ggsave(paste(file.name, ".jpeg", sep = "")))
    suppressMessages(ggsave(paste(file.name, ".pdf", sep = "")))
  })
})
rm(thr)

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
covered.seq.percentage <- NULL
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

heatmap.profiles <- NULL
comp.order.list.rows <- vector("list", 4)
comp.order.list.cols <- vector("list", 4)
domain <- seq(from = 0, to = 100, by = 5)

## Calculate the distance table between the co-ocurrences
dist.coocurrence <- Dist(covered.seq.percentage ,method = 'pearson')

## Create the heatmap using 4 agglomeration rules
th <- sapply(c("average", "complete", "single", "ward"), function(m){
  
  if(m == "ward"){
    temp <- m
    m <- "ward.D"
  }
  
  ## Calculate the hierarchical tree
  col.order <- hclust(dist.coocurrence, method = m)[[3]]
  
  if(m == "ward.D"){
    m <- "ward"
  }
  
  ## The col and rows have the same order
  comp.order.list.rows[[m]] <<- paste(col.order, collapse = ",")
  comp.order.list.cols[[m]] <<- paste(col.order, collapse = ",")
  
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

###################################################################
## Step 21: Add the columns Profile cluster to the feature table ##
###################################################################
list.counter <- 0
features.table <- lapply(features.table, function(ft){
  
  list.counter <<- list.counter + 1
  ft$Profile_cluster <- profile.clusters.names[[list.counter]]
  return(ft)
})

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
  order.by.eval <- order(as.numeric(as.vector(features.table[[list.counter]]$E_val)))
  features.table[[list.counter]] <- features.table[[list.counter]][order.by.eval,]
  
  ## Set the motif names and IDs
  TF.IDs <- matrix.names
  
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
  thrash <- apply(freq.per.bin[[list.counter]][order.by.eval,], 1, function(values){
    
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
  significance <- as.numeric(as.vector(features.table[[list.counter]]$Sig))
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
  datatable.info.tab <- features.table[[list.counter]]
  datatable.info.tab$P_val_threshold <- all.pval.match
  datatable.info.tab$ID <- TF.IDs[order.by.eval]
  datatable.info.tab$Names <- TF.names
  datatable.info.tab$Profiles <- profiles.plots
  datatable.info.tab$TFBS <- tfbss.plots
  datatable.info.tab$TFBS_per_seq <- tfbss.per.seq.plots
  datatable.info.tab$Logo <- logos.F
  datatable.info.tab$Logo_RC <- logos.R
  datatable.info.tab$covered_files <- covered.files
  datatable.info.tab$not_covered_files <- not.covered.files
  all.motifs <- all.motifs[[list.counter]]
  all.motif.names <- all.motif.names[[list.counter]]
  
  #############################################################
  ## Fill the HTML template                                  ##
  ## Substitute the words marked in the template by the data ##
  #############################################################
  html.report <- readLines(html.template.file)
  verbose(paste("Fill html report"), 1)
  
  # [1] "Feature"           "Chi_squared"       "Degrees"           "Sig"               "P_val"            
  # [6] "E_val"             "Q_val"             "Coverage"          "Sequences"         "Nb_hits"          
  # [11] "Covered_seq"       "Not_covered_seq"   "Profile_cluster"   "P_val_threshold"   "ID"               
  # [16] "Names"             "Profiles"          "TFBS"              "TFBS_per_seq"      "Logo"             
  # [21] "Logo_RC"           "covered_files"     "not_covered_files"
  # c("Feature", "ID", "P_val", "E_val", "Sig", "Q_val", "Nb_hits", "Sequences", "Coverage", "Chi_squared", "Profile_cluster", "Profiles", "TFBS", "TFBS_per_site", "Logo", "Logo_RC", "covered_files", "not_covered_files")
  # c(1,15,5,6,4,7,10,9,8,2,13,17,18,19,20,21,22,23)
  
  profile.data.tab.html <- create.html.tab(datatable.info.tab[,c(1,15,5,6,4,7,10,9,8,2,13,17,18,19,20,21,22,23)], img = c(15,16), plot = c(12,13,14), link.text.covered = 17, link.text.not.covered = 18)
  profile.data.tab.html <- gsub("Inf", "&infin;", profile.data.tab.html)
  
  profile.data.tab.html <- paste(profile.data.tab.html, collapse = "\n")
  html.report <- gsub("--tab--", profile.data.tab.html, html.report)
  
  ## Define the x-axis categories
  x.axis.categories <- paste(paste("'", colnames(freq.per.bin[[list.counter]]), "'", sep = ""), collapse = ",")
  
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
  cc <- as.numeric(gsub("%", "", features.table[[list.counter]]$Coverage))
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
  max.y <- max(freq.per.bin[[list.counter]]) + 0.02
  html.report <- gsub("--y_axis--", max.y, html.report)
  
  ## Fill the parameters table
  verbose(paste("Creating parameters table"), 1)
  html.report <- gsub("--bin_l--", bins[list.counter], html.report)
  html.report <- gsub("--bin_nb--", ncol(freq.per.bin[[list.counter]]), html.report)
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
  ss <- as.numeric(gsub("%", "", features.table[[list.counter]]$Sig))
  ss[ss == Inf] <- 350
  sig.cov <- paste("cov_significances['", all.motifs.cover[[list.counter]], "'] = ", as.vector(ss), ";", sep = "")
  sig.cov <- paste(sig.cov, collapse = "\n")
  html.report <- gsub("--significances_cov--", sig.cov, html.report)
  
  ## Add the coverages (to display in the tooltip)
  ## They are inserted in the JS section
  cc <- as.numeric(gsub("%", "", features.table[[list.counter]]$Coverage))
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
  verbose(paste("Exporting html report"), 1)
  html.report.file <- paste(basename(prefix), "_scan_profile_report.html", sep = "")
  write(html.report, file = html.report.file)
})
rm(binthrash)
# grep -v '^;' matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments.tab | awk -F '\t'  ' $2!="limit" && ($11 >= 4) {print $1"\t"$3"\t"($6+$5)/2"\t"$9} ' > matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments_PARSED.tab
# cat /home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/enrichment_by_scan/Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " matrix.scan.active = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; matrix.scan.inactive = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; p.val = '1e-4'; bin = '50'; pdf.file = './test.pdf'"


########################
## Centrimo algorithm ##
########################

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