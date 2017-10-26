#############################
## Step 1
## Load required libraries
#############################

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

##############################################
## Repeat N times each element of one array
repeat.n <- function(x = c(1,2,3), times = 1){
  
  x.n <- sapply(x, function(y){
    rep(y, times = times)
  })
  x.n <- as.vector(matrix(x.n, nrow = 1))
  return(x.n)
}

#######################################
## Create html table for the website
create.html.tab <- function(tab, img = 0, plot = 0){
  
  full.tab <- NULL
  head.tab <- "<div id='individual_motif_tab' style='width:2500px;display:none' class='tab div_chart_sp'><p style='font-size:12px;padding:0px;border:0px'><b>Individual Motif View</b></p><table id='Motif_tab' class='hover compact stripe' cellspacing='0' width='900px' style='padding:15px;align:center;'><thead><tr><th class=\"tab_col\"> Motif_name </th><th class=\"tab_col\"> Motif_ID </th> <th class=\"tab_col\"> Chi-squared </th> <th class=\"tab_col\"> Degrees </th> <th class=\"tab_col\"> Sig (Chi) </th> <th class=\"tab_col\"> Eval (Chi) </th> <th class=\"tab_col\"> Pval (Chi) </th><th class=\"tab_col\"> Qval (Chi) </th> <th class=\"tab_col\"> Coverage (query) </th><th class=\"tab_col\"> Nb sequences (query) </th><th class=\"tab_col\"> Coverage (control) </th><th class=\"tab_col\"> Nb sequences (control) </th><th class=\"tab_col\"> Nb hits (query) </th><th class=\"tab_col\"> Nb hits (control) </th> <th class=\"tab_col\"> Profile cluster </th><th class=\"tab_col\"> Profiles </th><th class=\"tab_col\"> Site distribution </th> <th class=\"tab_col\"> Logo </th> <th class=\"tab_col\"> Logo (RC) </th></tr></thead><tbody>"
  content.tab <- apply(tab, 1, function(row){
    
    row.length <- length(row)
    rows.nb <- 1:row.length
    
    ## Get the number of the columns with/without picture or plot
    ## This is done because the tab require different arguments
    rows.simple <- rows.nb[!(rows.nb %in% img)]
    rows.simple <- rows.simple[!(rows.simple %in% plot)]
    
    rows.pic <- rows.nb[rows.nb %in% img]
    rows.plot <- rows.nb[rows.nb %in% plot]
    
    ## Columns with simple text
    rows.text <- paste("<td>", row[rows.simple], "</td>", collapse = "")
    
    ## Columns with images
    rows.pic.text <- paste("<td><img class='logo_tab' src ='", as.character(row[rows.pic]), "'/></td>", collapse = "")
  
    ## Columns with plots and links
    rows.plot.pdf <- sapply(row[rows.plot], function(x){
      gsub("jpeg","pdf", x)
    })
    rows.plot.text <- paste("<td><a href='", rows.plot.pdf, "' target='_blank'><img class='plot_tab' src ='", row[rows.plot], "'/></a></td>", collapse = "")

    ## Head and tail tags
    row.head <- "<tr>"
    row.tail <- "</tr>"
    paste(row.head, rows.text, rows.plot.text, rows.pic.text, row.tail, sep = "")    

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
if (!exists("matrix.scan.file.query")) {
  stop("Missing mandatory argument (Query matrix-scan results table): matrix.scan.file.query ")
} else if (!exists("matrix.scan.file.control")) {
  stop("Missing mandatory argument (Control matrix-scan results table): matrix.scan.file.control ")
} else if (!exists("sequence.names.file.query")) {
  stop("Missing mandatory argument (Query sequence names table): sequence.names.file.query ")
} else if (!exists("sequence.names.file.control")) {
  stop("Missing mandatory argument (Control sequence names table): sequence.names.file.control ")
} else if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
} else if (!exists("ID.to.names.correspondence.tab")) {
  stop("Missing mandatory argument (Correspondence of the motif IDs and names): ID.to.names.correspondence.tab ")
} else if (!exists("seq.length")) {
  
  ## For the moment we assume all the input sequences have the same size
  stop("Missing mandatory argument (Sequence length): seq.length ")
}  else if (!exists("results.folder")){
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
  logo.folder <- paste(basename(prefix), "_logos", sep = "")
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
bin <- as.numeric(bin)
## Heatmap dendogram position
if (heatmap.dendo == "show"){
  heatmap.dendo <- "row"
} else if(heatmap.dendo == "hide"){
  heatmap.dendo <- "none"
}
print(heatmap.dendo)
bins <- bin

## Create folder for individual profile plots
dir.create(paste(prefix, "_TFBSs_positional_profiles/", sep = ""), recursive = TRUE, showWarnings = FALSE )

# matrix.scan.file.query <- "PRE_cluster_1_2_3_4_vs_PRE_cluster_6_matrix_scan_results_PARSED.tab"
# matrix.scan.file.control <- "PRE_cluster_1_2_3_4_vs_PRE_cluster_6_matrix_scan_results_PARSED_control.tab"
# sequence.names.file.query <- "PRE_cluster_1_2_3_4_vs_PRE_cluster_6_matrix_scan_sequence_names.tab"
# sequence.names.file.control <- "PRE_cluster_1_2_3_4_vs_PRE_cluster_6_matrix_scan_sequence_names_control.tab"
# ID.to.names.correspondence.tab <- "PRE_cluster_1_2_3_4_vs_PRE_cluster_6_TF_ID_name_correspondence.tab"
# bin <- 50
# seq.length <- 600

###########################################################
## Step 3: Read matrix-scan tables for query and control ##
###########################################################

## Read matrix-scan table (Query)
verbose(paste("Reading Query matrix-scan results table"), 1)
matrix.scan.results.query <- read.csv(file = matrix.scan.file.query, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results.query) <- c("seq_id", "ft_name", "bspos", "Pval")

## Read matrix-scan table (Control)
verbose(paste("Reading Control matrix-scan results table"), 1)
matrix.scan.results.control <- read.csv(file = matrix.scan.file.control, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results.control) <- c("seq_id", "ft_name", "bspos", "Pval")

## Remove duplicated rows
## One row can be duplicated because some motifs are palindromic and they could have
## an identical match in the same position for the F and R strand
matrix.scan.results.query <- matrix.scan.results.query[!duplicated(matrix.scan.results.query), ]
matrix.scan.results.control <- matrix.scan.results.control[!duplicated(matrix.scan.results.control), ]

#############################################################
## Step 4: Read sequence names table for query and control ##
#############################################################
## Read sequence names table (Query)
verbose(paste("Reading sequence names table"), 1)
sequence.names.tab.query <- read.csv(file = sequence.names.file.query, sep = "\t", header = TRUE, comment.char = ";")
colnames(sequence.names.tab.query) <- c("seq_id")
total.scanned.sequences.query <- length(as.vector(sequence.names.tab.query$seq_id))
scanned.sequences.query <- unique(as.vector(sequence.names.tab.query$seq_id))

#######################################
## Read sequence names table (Control)
verbose(paste("Reading sequence names table"), 1)
sequence.names.tab.control <- read.csv(file = sequence.names.file.control, sep = "\t", header = TRUE, comment.char = ";")
colnames(sequence.names.tab.control) <- c("seq_id")
total.scanned.sequences.control <- length(as.vector(sequence.names.tab.control$seq_id))
scanned.sequences.control <- unique(as.vector(sequence.names.tab.control$seq_id))

#################
## Set p-value
p.val <- as.numeric(p.val)

#################################################################
## Step 5: Create the column -log10(pvalue)                    ##
## Assign a class to each p-value, and one color to each class ##
## To later plot the qualitative distribution of TFBSs         ##
#################################################################
classes.pval <- list()
classes.pval.letters <- list()
for(i in 1:2){
  
  if(i == 1){
    ms.tab <- matrix.scan.results.query
  } else {
    ms.tab <- matrix.scan.results.control
  }
  
  ms.tab$Pval.minlog10 <- -log10(ms.tab$Pval)
  ms.tab$Pval.class <- ceiling(ms.tab$Pval.minlog10*2)/2
  
  classes.pval[[i]] <- sort(unique(ms.tab$Pval.class))
  classes.pval.letters[[i]] <- LETTERS[1:length(classes.pval[[i]])]
  
  ms.tab$Pval.class.letter <- sapply(ms.tab$Pval.class, function(x){
    p.class <- which(classes.pval[[i]] == x)
    classes.pval.letters[[i]][p.class ]
  })
  
  if(i == 1){
    matrix.scan.results.query <- ms.tab
  } else {
    matrix.scan.results.control <- ms.tab
  }
}
min.pval.minus.log10.query <- min(matrix.scan.results.query$Pval.minlog10)
max.pval.minus.log10.query <- max(matrix.scan.results.query$Pval.minlog10) 
min.pval.minus.log10.control <- min(matrix.scan.results.control$Pval.minlog10)
max.pval.minus.log10.control <- max(matrix.scan.results.control$Pval.minlog10)

class.length <- sapply(classes.pval, length)
if(class.length[1] >= class.length[2]){
  classes.pval <- classes.pval[[1]]
  classes.pval.letters <- classes.pval.letters[[1]]
} else {
  classes.pval <- classes.pval[[2]]
  classes.pval.letters <- classes.pval.letters[[2]]
}
 
#############################
## Get the sequences + motif name's
seq.id.query <- unique(as.vector(matrix.scan.results.query$seq_id))
seq.id.control <- unique(as.vector(matrix.scan.results.control$seq_id))

matrix.names.query <- unique(as.vector(matrix.scan.results.query$ft_name))
matrix.names.control <- unique(as.vector(matrix.scan.results.control$ft_name))
matrix.names <- intersect(matrix.names.query, matrix.names.control)
nb.motifs <- length(matrix.names.query)
verbose(paste(length(matrix.names), "motifs found in the intersection of Query vs Control sequences"), 1)

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
verbose(paste("Creating plots with distribution of TFBSs at different p-values"), 1)

max.pval <- max(c(matrix.scan.results.query$Pval.minlog10, matrix.scan.results.control$Pval.minlog10))

thr <- sapply(1:nb.motifs, function(m){

  ## Get the matrix name
  motif <- matrix.names[m]
  motif.selected.name <- as.vector(ID.names[,2][which(ID.names[,1] == motif)][1])

  ## Get the sub-table with the hits of the query matrix
  motif.selection.query <- matrix.scan.results.query[matrix.scan.results.query$ft_name == motif,]
  motif.selection.query$bspos <- motif.selection.query$bspos + limits
  motif.classes.query <- sort(unique(motif.selection.query$Pval.class.letter))

  motif.selection.control <- matrix.scan.results.control[matrix.scan.results.control$ft_name == motif,]
  motif.selection.control$bspos <- motif.selection.control$bspos + limits
  motif.classes.control <- sort(unique(motif.selection.control$Pval.class.letter))

  ################################
  #Create a custom color scale
  myColors <- colorRampPalette(brewer.pal(nb.color.classes, "YlGnBu"), space="Lab")(length(classes.pval.letters))
  names(myColors) <- classes.pval.letters

  ## Insert logo
  logo.file <- paste(logo.folder, "/", motif, "_logo.png", sep = "")
  logo <- readPNG(logo.file)
  logo.roster <- rasterGrob(logo, interpolate = TRUE)

  ggplot.list <- list()
  TFBSs.pval.distribution.file <- paste(basename(prefix), "_TFBSs_pval_distribution/", motif, "_TFBSs_pval_classes", sep = "")
  for(i in 1:2){

    ## Select the matrix-scan sub-table
    if(i == 1){
      seq.type <- "Query"
      motif.selection <- motif.selection.query
    } else {
      seq.type <- "Control"
      motif.selection <- motif.selection.control
    }
    ## Range of p-values for the query motif
    pval.class.motif <- sort(unique(motif.selection$Pval.class))

    nb.TFBSs <- length(as.vector(motif.selection$bspos))
    nb.seq <- unique(as.vector(motif.selection$seq_id))

    ## X position of plot annotations
    text.xmax <- min(motif.selection$bspos) + max(motif.selection$bspos) / 4
    text.center <- (min(motif.selection$bspos) - text.xmax)*2

    ggplot.list[[i]] <- ggplot(motif.selection, aes(x=bspos, y=Pval.minlog10)) +
      ylim(c(min(matrix.scan.results.query$Pval.minlog10), max.pval)) +
      geom_point(aes(colour = Pval.class.letter), shape = 18, size = 4, stroke = 0.5) +
      # geom_rug(position='jitter') +
      labs(title=paste("Qualitative distribution of ", motif, " TFBSs\nin ", seq.type, " sequences", sep = ""), y = "-log10(P-value)", x = "Position") +
      scale_colour_manual(name = "-log10(P-value)",values = myColors, labels = paste(">", pval.class.motif, sep = "")) +
      theme(
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey"),
        panel.ontop = FALSE
      )
      # annotate("text", x = -limits + ((limits*2)/10), y = max.pval - 0.25, label = paste("Nb of TFBSs: ", nb.TFBSs, sep = ""), size = 4, hjust = 0) +
      # annotate("text", x = -limits + ((limits*2)/10), y = max.pval - 0.55, label = paste("Nb of sequences: ", nb.seq, sep = ""), size = 4, hjust = 0) +
      # annotation_custom(logo.roster, xmax = limits - (limits/3), xmin = limits - 5, ymin = max.pval - 1, ymax = max.pval - 0.05)
  }

  xx <- plot_grid(ggplot.list[[1]], ggplot.list[[2]], labels=c("", ""), ncol = 2, nrow = 1, rel_widths = c(1,1))
  save_plot(filename = paste(TFBSs.pval.distribution.file, ".pdf", sep = ""), plot = xx, ncol = 2, base_aspect_ratio = 2, base_width = 6, base_height = 6.5)
  save_plot(filename = paste(TFBSs.pval.distribution.file, ".jpeg", sep = ""), plot = xx, ncol = 2, base_aspect_ratio = 2, base_width = 6, base_height = 6)
})
rm(thr)

#################################################################################
## Step 8: Calculate the raw counts of TFBSs on each position of the sequences ##
#################################################################################
verbose(paste("Getting the raw counts of TFBSs in the sequences"),1)

matrix.scan.results.query$bspos.left <- matrix.scan.results.query$bspos + seq.length
matrix.scan.results.control$bspos.left <- matrix.scan.results.control$bspos + seq.length


get.raw.counts <- function(m, df){
  sites.m <- as.vector(df[which(df$ft_name == m),"bspos.left"])
  raw.counts.m <- tabulate(sites.m, nbins = seq.length)
  return(raw.counts.m)
}

## Get a DataFrame with the raw counts of TFBSs per each nucleotide of the sequence
raw.counts.all.motifs.query <- sapply(matrix.names, get.raw.counts, matrix.scan.results.query)
raw.counts.all.motifs.query <- as.data.frame(raw.counts.all.motifs.query)

raw.counts.all.motifs.control <- sapply(matrix.names, get.raw.counts, matrix.scan.results.control)
raw.counts.all.motifs.control <- as.data.frame(raw.counts.all.motifs.control)
  
#######################################################
## Step 9: Calculate the raw counts of TFBSs per bin ##
#######################################################
## Define the bin 'names'
x.neg <- seq(-(seq.length/2), -(bins) , bins)
xlab <-append(x.neg, rev(-x.neg))

## For each bin, sums the number of TFBSs (function rollaply)
## The results are stored in a list where each element correspond to the counts of a given bin size
verbose(paste("Separating the raw counts of TFBSs in bins of size:", bin),1)
counts.per.bin.query <- rollapply(raw.counts.all.motifs.query, bins, sum, na.rm = TRUE, by = bins, partial = TRUE, align="left")
colnames(counts.per.bin.query) <- matrix.names
rownames(counts.per.bin.query) <- xlab

counts.per.bin.control <- rollapply(raw.counts.all.motifs.control, bins, sum, na.rm = TRUE, by = bins, partial = TRUE, align="left")
colnames(counts.per.bin.control) <- matrix.names
rownames(counts.per.bin.control) <- xlab

#####################################################
## Step 10: Chi-squared and KS calculation section ##
#####################################################
verbose(paste("Calculating Chi-Squared and KS statistics"),1)
counts.per.bin.query <- t(counts.per.bin.query)
counts.per.bin.control <- t(counts.per.bin.control)
all.chi.ks.results <- sapply(1:nrow(counts.per.bin.query), function(r){

  #####################################
  ## Calculate the X2 in two senses:
  ## Query vs Control
  ## Control vs Query
  query <- counts.per.bin.query[r,]
  query.freq <- query/sum(query)
  control <- counts.per.bin.control[r,]
  control.freq <- control/sum(control)
  
  chisq.query.vs.control <- chisq.test(query, p = control.freq, correct = TRUE)
  chisq.control.vs.query <- chisq.test(control, p = query.freq, correct = TRUE)
  
  # #####################################
  # ## Calculate the KS in two senses:
  # ## Query vs Control
  # ## Control vs Query
  # ## Actually both return the same result!
  # ks.query.vs.control <- ks.test(unique(query.freq), unique(control.freq), alternative = "two.sided")
  # # ks.control.vs.query <- ks.test(control, query, alternative = "two.sided")
  # 
  # return(c(chisq.query.vs.control[[1]], chisq.query.vs.control[[3]], chisq.control.vs.query[[1]], chisq.control.vs.query[[3]], chisq.control.vs.query[[2]], round(ks.query.vs.control[[1]], digits = 2), ks.query.vs.control[[2]], bins))
  return(c(chisq.query.vs.control[[1]], chisq.query.vs.control[[3]], chisq.control.vs.query[[1]], chisq.control.vs.query[[3]], chisq.control.vs.query[[2]], bins))
  
})
all.chi.ks.results <- t(all.chi.ks.results)
rownames(all.chi.ks.results) <- matrix.names
colnames(all.chi.ks.results) <- c("Chi_Q_vs_C", "Chi_Pvalue_Q_vs_C", "Chi_C_vs_Q", "Chi_Pvalue_C_vs_Q", "DF", "Bin_size")


#############################################
## Step 11: Calculate the Frequency tables ##
#############################################
verbose(paste("Exporting counts and frequency tables"),1)
dir.create(paste(basename(prefix), "_tables/", sep = ""), recursive = TRUE, showWarnings = FALSE )

## Number of Hits per motif (query and control)
counts.per.bin.query.nb.hits <- apply(counts.per.bin.query, 1, sum)
counts.per.bin.control.nb.hits <- apply(counts.per.bin.control, 1, sum)

## Frequency of TFBSs (query and control)
counts.per.bin.query.freq <- counts.per.bin.query/counts.per.bin.query.nb.hits
counts.per.bin.control.freq <- counts.per.bin.control/counts.per.bin.control.nb.hits


query.counts.tab.file <- paste(basename(prefix), "_tables/counts_per_bin_profiles_bin", bins, "_query.tab", sep = "")
control.counts.tab.file <- paste(basename(prefix), "_tables/counts_per_bin_profiles_bin", bins, "_control.tab", sep = "")
fraction.tab.file.query <- paste(basename(prefix), "_tables/density_per_bin_profiles_bin", bins, "_query.tab", sep = "")
fraction.tab.file.control <- paste(basename(prefix), "_tables/density_per_bin_profiles_bin", bins, "_control.tab", sep = "")
for(i in 1:2){
  
  if(i == 1){
    freq.tab.file <- fraction.tab.file.query
    counts.tab.file <- query.counts.tab.file
    fpb <- counts.per.bin.query.freq
    cpb <- counts.per.bin.query
  } else {
    freq.tab.file <- fraction.tab.file.control
    counts.tab.file <- control.counts.tab.file
    fpb <- counts.per.bin.control.freq
    cpb <- counts.per.bin.control
  }
  
  write.table(cpb, file = counts.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  write.table(fpb, file = freq.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
  
}

## Calculate the highest frequency of TFBSs
max.y <- max(as.vector(c(counts.per.bin.query.freq, counts.per.bin.control.freq)), na.rm = TRUE)


#################################
## Step 12: calculate coverage ##
#################################
verbose(paste("Calculating coverage and number of hits per motif"),1)
coverage.all.motifs <- sapply(matrix.names, function(m){
  
  matched.query.sequences <- as.vector(subset(matrix.scan.results.query, ft_name == m)$seq_id)
  unique.matched.query.sequences <- unique(matched.query.sequences)
  m1 <- round(length(unique.matched.query.sequences)/total.scanned.sequences.query, digits = 2)
  m2 <- length(unique.matched.query.sequences)
  m3 <- length(matched.query.sequences)

  matched.control.sequences <- as.vector(subset(matrix.scan.results.control, ft_name == m)$seq_id)
  unique.matched.control.sequences <- unique(matched.control.sequences)
  m4 <- round(length(unique.matched.control.sequences)/total.scanned.sequences.control, digits = 2)
  m5 <- length(unique.matched.control.sequences)
  m6 <- length(matched.control.sequences)
  
  return(list(m1,m2,m3,m4,m5,m6))

})
coverage.query <- unlist(coverage.all.motifs[1,])
seq.matched.query <- unlist(coverage.all.motifs[2,])
nb.hits.query <- unlist(coverage.all.motifs[3,])
coverage.control <- unlist(coverage.all.motifs[4,])
seq.matched.control <- unlist(coverage.all.motifs[5,])
nb.hits.control <- unlist(coverage.all.motifs[6,])
rm(coverage.all.motifs)

# plot(x = coverage.query,
#      y = coverage.control,
#      xlim = c(0,1),
#      ylim = c(0,1))
# lines(x = seq(0,1,0.001), y = seq(0,1,0.001))

########################################################################
## Step 13: calculate Q-value, E-value, and significance -log(pvalue) ##
## Re-order the table                                                 ##
########################################################################
verbose(paste("Calculating  E-value, Q-value and significance"),1)
df <- as.data.frame(all.chi.ks.results)

## Round Chi values
df$Chi_Q_vs_C <- round(df$Chi_Q_vs_C, digits = 2)
df$Chi_C_vs_Q <- round(df$Chi_C_vs_Q, digits = 2)

## Calculate E-values
df$Chi_Evalue_Q_vs_C <- df$Chi_Pvalue_Q_vs_C * nb.motifs
df$Chi_Evalue_C_vs_Q <- df$Chi_Pvalue_C_vs_Q * nb.motifs

## Calculate Significance
df$Sig_Chi_Q_vs_C <- round(-log10(df$Chi_Evalue_Q_vs_C), digits = 2)
df$Sig_Chi_C_vs_Q <- round(-log10(df$Chi_Evalue_C_vs_Q), digits = 2)

## Calculate Q-values
df$Chi_Qvalue_Q_vs_C <- p.adjust(as.vector(df$Chi_Pvalue_Q_vs_C), method = "BH")
df$Chi_Qvalue_C_vs_Q <- p.adjust(as.vector(df$Chi_Pvalue_C_vs_Q), method = "BH")

####################################################
## Transform Pval, Eval and Qval to PrettyNumbers ##
####################################################
df$Chi_Pvalue_Q_vs_C <- prettyNum(df$Chi_Pvalue_Q_vs_C, scientific=TRUE, digits = 2)
df$Chi_Pvalue_C_vs_Q <- prettyNum(df$Chi_Pvalue_C_vs_Q, scientific=TRUE, digits = 2)

df$Chi_Evalue_Q_vs_C <- prettyNum(df$Chi_Evalue_Q_vs_C, scientific=TRUE, digits = 2)
df$Chi_Evalue_C_vs_Q <- prettyNum(df$Chi_Evalue_C_vs_Q, scientific=TRUE, digits = 2)

df$Chi_Qvalue_Q_vs_C <- prettyNum(df$Chi_Qvalue_Q_vs_C, scientific=TRUE, digits = 2)
df$Chi_Qvalue_C_vs_Q <- prettyNum(df$Chi_Qvalue_C_vs_Q, scientific=TRUE, digits = 2)

## Add the matrix Id as feature ID
df$feature <- rownames(df)

## Add the coverage columns (for query and control sequences)
df$coverage_query <- coverage.query
df$covered_seq_query <- seq.matched.query
df$coverage_control <- coverage.control
df$covered_seq_control <- seq.matched.control
df$nb_hits_query <- nb.hits.query
df$nb_hits_control <-nb.hits.control

features.table <- df[, c("feature", "Chi_Q_vs_C", "DF","Sig_Chi_Q_vs_C", "Chi_Evalue_Q_vs_C", "Chi_Pvalue_Q_vs_C", "Chi_Qvalue_Q_vs_C", "coverage_query", "covered_seq_query", "coverage_control", "covered_seq_control", "nb_hits_query", "nb_hits_control")]
colnames(features.table) <- c("Feature", "Chi_squared", "Degrees", "Sig_Chi", "E_val_Chi", "P_val_Chi", "Q_val_Chi","Coverage_query", "Sequences_query", "Coverage_control", "Sequences_control","Nb_hits_query", "Nb_hits_control")

###############################################################################
## Step 14: calculate -log2 ratio between Query and Control TFBS frequencies ##
## Print the positional profiles comparing query and control sequences       ##
###############################################################################
verbose(paste("Drawing positional profiles"),1)
counts.per.bin.log2 <- sapply(matrix.names, function(m){
  
  ## Set Output file name
  file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", m, "_positional_profile", sep = "")

  ## Load the logo
  logo.file <- paste(logo.folder, "/", m, "_logo.png", sep = "")
  logo <- readPNG(logo.file)
  logo.roster <- rasterGrob(logo, interpolate = TRUE)

  ## Plot profiles
  control.points <- data.frame(y = counts.per.bin.control.freq[m,], x = xlab, Sequences = rep("control", times = length(counts.per.bin.control.freq[m,])))
  query.points <- data.frame(y = counts.per.bin.query.freq[m,], x = xlab, Sequences = rep("query", times = length(counts.per.bin.query.freq[m,])))
  df.freq <- rbind(control.points, query.points)
  df.freq$Sequences <- as.factor(df.freq$Sequences)

  ggplot(df.freq, aes(y=y, x=x, group = Sequences, colour=Sequences)) +
    geom_line(size = 2) +
    ylim(0, max.y) +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = "grey"),
      panel.ontop = FALSE
    ) +
    geom_rug(position='jitter', sides="l") +
    labs(title=paste(m, " binding profile (Query vs Control)", sep = ""), y = "Frequency of TFBSs", x = "Position")
    # annotation_custom(logo.roster, xmax = limits, xmin = limits - sum(abs(range(df.freq$x)))/5, ymin = max.y - 0.01, ymax = max.y - 0.075)

  ## Export the file
  suppressMessages(ggsave(paste(file.name, ".jpeg", sep = "")))
  suppressMessages(ggsave(paste(file.name, ".pdf", sep = "")))

  ## Calculate the -log2 ratio
  -log2(counts.per.bin.query.freq[m,]/counts.per.bin.control.freq[m,])
  
})


######################################
## Step 15: Define profile clusters ##
######################################
verbose(paste("Clustering of differential binding profiles"),1)
cluster.profiles.motifs <- NULL
cluster.profiles.motif.names <- NULL
data.t <- t(counts.per.bin.log2)
tree.profiles <- hclust(Dist(data.t, method = "correlation"), method = "ward.D2")
clusters.tree.profiles <- cutreeDynamic(tree.profiles, minClusterSize = 1, method = "tree")
names(clusters.tree.profiles) <- matrix.names

## Generate a color palette
nb.profile.clusters <- length(unique(clusters.tree.profiles))
cluster.tree.profiles.palette <- colorRampPalette(brewer.pal(9, "Set1"), space="Lab")(nb.profile.clusters)

## Assign a different color to each cluster
color.clusters.tree.profiles <- as.vector(sapply(clusters.tree.profiles+1, function(color){
  cluster.tree.profiles.palette[color]
}))

profile.clusters.names <- as.vector(sapply(matrix.names, function(m){
  profile.cluster <- as.vector(clusters.tree.profiles[m])
  paste("Profile_cluster_", profile.cluster, sep = "")
}))

## Get the member motif IDs of each cluster
cluster.profiles.counter <- 0
thrash <- sapply(1:nb.profile.clusters, function(cl){
  
  cluster.profiles.counter <<- cluster.profiles.counter + 1
  cluster.profiles.motifs[[cluster.profiles.counter]] <<- names(which(clusters.tree.profiles == cluster.profiles.counter))
  
  ## Get the motif name
  cluster.profiles.motif.names[[cluster.profiles.counter]] <<- as.vector(
    sapply(cluster.profiles.motifs[[cluster.profiles.counter]], function(n){

      
      as.vector(ID.names[which(ID.names[,1] == n),2])
    })
  )
})
rm(thrash)

########################################
## Step 16: Draw -log2 ratio heatmaps ##
########################################

# data.t[is.infinite(data.t)] <- 0
# 
# verbose(paste("Drawing differential binding heatmap"),1)
# ## Create color palette
# heatmap.color.classes <- 11
# heatmap.color.palette <- "RdBu"
# rgb.palette <- rev(colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")(heatmap.color.classes))
# 
# hm.main <- paste("Profile Heatmap\nQuery over Control")
# 
# for(pf in print.formats){
# 
#   if(pf == "pdf"){
#     pdf.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", "Query_over_control_positional_profile.pdf", sep = "")
#     pdf(pdf.file.name)
#   } else {
#     jpeg.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", "Query_over_control_positional_profile.jpeg", sep = "")
#     jpeg(jpeg.file.name)
#   }
#   data.t <- as.matrix(t(counts.per.bin.log2))
# 
#   heatmap.2(data.t,
# 
#             main = hm.main,
#             xlab = "Position (bp)",
#             ylab = "Motifs",
# 
#             ## The order of the values is set according these dendrograms
#             Rowv = as.dendrogram(tree.profiles),
#             Colv = FALSE,
#             dendrogram = "row",
# 
#             ## Color
#             col = rgb.palette,
# 
#             ## Trace
#             trace = "none",
# 
#             ## Side colors
#             RowSideColors = color.clusters.tree.profiles,
# 
#             ## Key control
#             key = TRUE,
#             keysize = 1,
#             density.info = "none",
#             key.xlab = "Log2 Ratio",
#             key.ylab = "",
#             key.title = "",
#             cexRow = 0.66,
#             offsetCol = 0.25
#   )
#     t <- dev.off()
# }

######################################################################
## Step 17: Add the columns Profile cluster + logos + profile plots ##
######################################################################
verbose(paste("Exporting attributes table"),1)
#####################################
## Write the logo and profile path
logos.F <- sapply(matrix.names, function(i){
  paste(logo.folder, "/", i, "_logo.png", sep = "")
})

logos.R <- sapply(matrix.names, function(i){
  paste(logo.folder, "/", i, "_logo_rc.png", sep = "")
})

## Write the Profile and TFBSs plots path
profiles.plots <- sapply(matrix.names, function(i) {
  paste(basename(prefix), "_TFBSs_positional_profiles/", i, "_positional_profile.jpeg", sep = "")
})

## Write the Profile and TFBSs plots path
pval.site.plots <- sapply(matrix.names, function(i) {
  paste(basename(prefix), "_TFBSs_pval_distribution/", i, "_TFBSs_pval_classes.jpeg", sep = "")
})

############################
## Get the IDs of the TFs
TF.IDs <- as.vector(sapply(matrix.names, function(x){
  as.vector(ID.names[which(ID.names[,1] == x),2])
}))

## As some IDs include a period ('.') in their text, we change it by
## an underscore ('_')
TF.IDs.cp <- gsub("\\.", "", TF.IDs)

features.table$ID <- as.vector(TF.IDs.cp)
features.table$Profile_cluster <- profile.clusters.names

feature.attributes.file <- paste(basename(prefix), "_attributes.tab", sep = "")
write.table(features.table, file = feature.attributes.file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

features.table$Profiles <- profiles.plots
features.table$Logo <- logos.F
features.table$Logo_RC <- logos.R
features.table$Site_distrib <- pval.site.plots

## Order the table according the Significance (-log10(E-value))
order.by.eval <- order(as.numeric(as.vector(features.table$Sig_Chi)))
features.table <- features.table[order.by.eval,]

print(dim(features.table))

features.table <- features.table[,c(1,14,2:13,15:19)]

################################
## Create dynamic html report ##
################################
verbose(paste("Creating HTML dynamic report"),1)
## Set colors
set.colors <- colorRampPalette(brewer.pal(10,"Paired"))(length(TF.IDs))
set.colors <- sapply(sample(set.colors), function(x){
  rep(x, times = 2)
})
set.colors <- as.vector(matrix(set.colors, nrow = 1))

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
TF.names <- matrix.names

thrash <- sapply(order.by.eval, function(o){
  
  counter <<- counter + 1
  
  ## Get the counts and the fraction of TFBSs per bin in query and control
  counts.query <- round(counts.per.bin.query.freq[o,], digits = 3)
  counts.control <- round(counts.per.bin.control.freq[o,], digits = 3)
  
  ## Here we create a unique ID without CSS special characters
  ## Only to manipulate the objects in the HTML form
  motif <- paste(counter, "_", TF.IDs.cp[counter], "_", counter, sep = "")
  motif.c <- paste(counter, "_", TF.IDs.cp[counter], "_control_", counter, sep = "") 
  
  motif <- gsub("_", "", motif)
  motif <- gsub("-", "", motif)
  motif <- gsub("\\.", "", motif)
  motif <- gsub(":", "", motif)
  motif <- gsub("\\s+", "", motif, perl = TRUE)
  motif.c <- gsub("_", "", motif.c)
  motif.c <- gsub("-", "", motif.c)
  motif.c <- gsub("\\.", "", motif.c)
  motif.c <- gsub(":", "", motif.c)
  motif.c <- gsub("\\s+", "", motif.c, perl = TRUE)
  
  hash.motif.IDs[[TF.IDs.cp[counter]]] <<- motif
  
  all.motifs <<- append(all.motifs, motif)
  all.motifs <<- append(all.motifs, motif.c)
  all.motif.names <<- append(all.motif.names, TF.names[counter])
  all.motif.names <<- append(all.motif.names, paste(TF.names[counter], "_control", sep = ""))
  
  ## Create the name's data for the HTML file
  ## The name correspond to the motif name. 
  ## Note that two motifs for the same TF will have the same name and this name will appears twice 
  ## in the report. However their IDs are unique.
  plot.names <<- append(plot.names, paste("'", motif, "' : '",  TF.names[counter],"',", sep = ""))
  plot.names <<- append(plot.names, paste("'", motif.c, "' : '",  TF.names[counter],"_control',", sep = ""))
  
  ## Create the X data
  ## For each TF we will have a row containing the ID plus the counts per bin
  ## The X-axis data for the plot
  y <- paste("['", 
             motif, 
             "',",
             paste(counts.query, 
                   collapse = ","),
             "],",
             sep = "")
  x.y <<- rbind(x.y, y) 
  
  y.c <- paste("['", 
             motif.c, 
             "',",
             paste(counts.control, 
                   collapse = ","),
             "],",
             sep = "")
  x.y <<- rbind(x.y, y.c) 
  
})

## Set the line width according the significance -log10(E-value)
## Higher significance means a wider line
significance <- as.numeric(as.vector(features.table$Sig_Chi))
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
  w <- rep(w, times = 2)
})
line.w <- as.vector(matrix(line.w, ncol = 1))


##############################################
## Insert the cluster functions and buttons
all.cluster.functions <- vector()
all.cluster.buttons <- vector()

## JS code to show the motifs corresponding to one cluster
## We require one function for each cluster
show.profile.cluster.function <- 'function --profile_cluster_show--() {
chart.hide([--all_motifs_function--]);
chart.show([--names_function--]);
}'
show.profile.cluster.button <- "<button class='small button_chart' onclick='--function_name--();'>--cluster_name--</button>"

## Get the member motif IDs of each cluster
cluster.profiles.motifs <- list()
cluster.profiles.motif.names <- list()
cluster.profiles.counter <- 0
thrash <- sapply(1:nb.profile.clusters, function(cl){

  cluster.profiles.counter <<- cluster.profiles.counter + 1
  cluster.profiles.motifs[[cluster.profiles.counter]] <<- names(which(clusters.tree.profiles+1 == cluster.profiles.counter))
  
  ## Get the motif name
  cluster.profiles.motif.names[[cluster.profiles.counter]] <<- as.vector(
    sapply(cluster.profiles.motifs[[cluster.profiles.counter]], function(n){
      
      ID.names[which(ID.names[,1] == n),2]
    })
  )
})
rm(thrash)

all.motifs.function <- paste(paste("'", all.motifs, "'", sep = ""), collapse = ",")

## Generate the JS function to show the clusters
thrash <- sapply(1:length(cluster.profiles.motif.names), function(cl){
  
  ## Define the profile cluster name
  profile.cluster.name <- paste("Profile_cluster_", cl, sep = "")
  
  cluster.function <- show.profile.cluster.function
  
  ##############################################################
  ## Change the function's name for the corresponding cluster
  
  ## Button name
  cluster.function.names <- paste(profile.cluster.name, "_show", sep = "")
  
  ## Cluster member names
  cluster.member.names <- as.vector(unlist(sapply(cluster.profiles.motif.names[[cl]], function(m){
    
    m <- gsub("_", "", m)
    m <- gsub("-", "", m)
    m <- gsub("\\.", "", m)
    m <- gsub(":", "", m)
    m <- gsub("\\s+", "", m, perl = TRUE)
    
    hash.motif.IDs[[m]]
  })))
  cluster.member.names <- paste(paste("'", cluster.member.names, "'", sep = ""), collapse = ",")
  
  ## Substitution
  cluster.function <- gsub("--profile_cluster_show--", cluster.function.names, cluster.function)
  cluster.function <- gsub("--all_motifs_function--", all.motifs.function, cluster.function)
  cluster.function <- gsub("--names_function--", cluster.member.names, cluster.function)
  
  ## Concat all the functions
  all.cluster.functions <<- append(all.cluster.functions, cluster.function)
  
  ##########################################################
  ## Create the button to show all the motifs in a cluster
  cluster.button <- show.profile.cluster.button
  cluster.button <- gsub("--cluster_name--", profile.cluster.name, cluster.button)
  cluster.button <- gsub("--function_name--", cluster.function.names, cluster.button)

  ## Concat all the functions
  all.cluster.buttons <<- append(all.cluster.buttons, cluster.button)
})
all.cluster.functions <- paste(all.cluster.functions, collapse = "\n")
all.cluster.buttons <- paste(all.cluster.buttons, collapse = "\n")


############################
## Fill the HTML template
## Substitute the words marked in the template by the data
html.report <- readLines(html.template.file)
print(colnames(features.table))
profile.data.tab.html <- create.html.tab(features.table, img = c(17,18), plot = c(16,19))
profile.data.tab.html <- gsub("Inf", "&infin;", profile.data.tab.html)
profile.data.tab.html <- paste(profile.data.tab.html, collapse = "\n")
html.report <- gsub("--tab--", profile.data.tab.html, html.report)

## Define the x-axis categories
x.axis.categories <- paste(paste("'", xlab, "'", sep = ""), collapse = ",")

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

## Create the sahed lines for the control
all.motifs.control <- all.motifs[grep("control", all.motifs)]
dashed <- paste("'", all.motifs.control, "': [{'end':'", xlab[length(xlab)], "'}],", sep = "")
dashed <- paste(dashed, collapse = "\n")
html.report <- gsub("--dashed--", dashed, html.report)

## Add the e-values data
## They are inserted in the JS section
eval.rep <- repeat.n(as.vector(features.table$E_val_Chi), times = 2)
evalues <- paste("evalues['", all.motifs, "'] = '", eval.rep, "';", sep = "")
evalues <- paste(evalues, collapse = "\n")
html.report <- gsub("--evalues--", evalues, html.report)

## Add the p-values (to display in the tooltip)
## They are inserted in the JS section
pval.rep <- repeat.n(as.vector(features.table$P_val_Chi), times = 2)
pvalues <- paste("pvalues['", all.motifs, "'] = '", pval.rep, "';", sep = "")
pvalues <- paste(pvalues, collapse = "\n")
html.report <- gsub("--pvalues--", pvalues, html.report)

## Add the profile clusters (to display in the tooltip)
## They are inserted in the JS section
cl.rep <- repeat.n(as.vector(features.table$Profile_cluster), times = 2)
profile.clusters.array <- paste(" profile_clusters['", all.motifs, "'] = '", cl.rep, "';", sep = "")
profile.clusters.array <- paste(profile.clusters.array, collapse = "\n")
html.report <- gsub("--profile_clusters_array--",  profile.clusters.array, html.report)

## Add the q-values (to display in the tooltip)
## They are inserted in the JS section
qval.rep <- repeat.n(as.vector(features.table$Q_val_Chi), times = 2)
qvalues <- paste("qvalues['", all.motifs, "'] = '", qval.rep, "';", sep = "")
qvalues <- paste(qvalues, collapse = "\n")
html.report <- gsub("--qvalues--", qvalues, html.report)

## Add the real motif IDs (to display in the tooltip)
## They are inserted in the JS section
## I called them 'real' because are those found on the original motif file
TF.IDs.rep <- repeat.n(TF.IDs, times = 2)
IDs <- paste("IDs['", all.motifs, "'] = '", TF.IDs.rep, "';", sep = "")
IDs <- paste(IDs, collapse = "\n")
html.report <- gsub("--IDs--", IDs, html.report)

## Add the real motif logo path (to display in the tooltip)
## They are inserted in the JS section
logos <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo.png", sep = "")
})
logos.rep <- repeat.n(as.vector(features.table$Logo), times = 2)
logos <- paste("pics['", all.motifs, "'] = '", logos.rep, "';", sep = "")
logos <- paste(logos, collapse = "\n")
html.report <- gsub("--pics--", logos, html.report)

## Logos in Reverse complement
logos.rc <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo_rc.png", sep = "")
})
logos.rc.rep <- repeat.n(as.vector(features.table$Logo_RC), times = 2)
logos.rc <- paste("pics_rc['", all.motifs, "'] = '", logos.rc.rep, "';", sep = "")
logos.rc <- paste(logos.rc, collapse = "\n")
html.report <- gsub("--pics_rc--", logos.rc, html.report)

## Add the signficance (to display in the tooltip)
## They are inserted in the JS section
sig.rep <- repeat.n(as.vector(features.table$Sig_Chi), times = 2)
sig <- paste("significances['", all.motifs, "'] = ", sig.rep, ";", sep = "")
sig <- paste(sig, collapse = "\n")
html.report <- gsub("--significances--", sig, html.report)

## Add the covertures (to display in the tooltip)
## They are inserted in the JS section
cc <- as.numeric(gsub("%", "", features.table$Coverage_query))
cc.rep <- repeat.n(as.vector(cc), times = 2)
coverture <- paste("TF_coverture['", all.motifs, "'] = ", cc.rep, ";", sep = "")
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
set.colors <- paste(paste("'", set.colors, "'", sep = ""), collapse = ",")
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
max.y <- max.y + 0.02
html.report <- gsub("--y_axis--", max.y, html.report)

## Fill the parameters table
html.report <- gsub("--seq_l--", seq.length, html.report)
html.report <- gsub("--bin_l--", bin, html.report)
html.report <- gsub("--bin_nb--", length(xlab), html.report)
html.report <- gsub("--seq_nb_q--", total.scanned.sequences.query, html.report)
html.report <- gsub("--seq_nb_c--", total.scanned.sequences.control, html.report)
html.report <- gsub("--motif_nb--", length(matrix.names), html.report)
html.report <- gsub("--p--", prettyNum(p.val), html.report)

## Fill the heatmap section
html.report <- gsub("--heatmap_png--", paste(basename, "/Query_over_control_positional_profile.jpeg", sep = ""), html.report)
html.report <- gsub("--heatmap_pdf--", paste(basename, "/Query_over_control_positional_profile.pdf", sep = ""), html.report)

## Insert the Hit counts per bin table
html.report <- gsub("--hit_counts_table_query--", query.counts.tab.file, html.report)
html.report <- gsub("--hit_counts_table_control--", control.counts.tab.file, html.report)

## Insert the density of hits per bin table
html.report <- gsub("--density_table_query--", fraction.tab.file.query, html.report)
html.report <- gsub("--density_table_control--", fraction.tab.file.control, html.report)

## Insert the Full motif description table
feature.attributes.file <- paste(basename, "_attributes.tab", sep = "")
html.report <- gsub("--full_motif_table--", feature.attributes.file, html.report)

## Insert the Full motif description table
html.report <- gsub("Inf;", "Infinity;", html.report)

## Insert the Full motif description table
min.sig <- min(as.numeric(as.vector(features.table$Sig_Chi)))
html.report <- gsub("--sig_min--", min.sig, html.report)

## Insert the JavaScript libraries path
html.report <- gsub("--c3_css--", c3.css.base, html.report)
html.report <- gsub("--c3--", c3.base, html.report)
html.report <- gsub("--d3--", d3.base, html.report)
html.report <- gsub("--jquery--", jquery.base, html.report)
html.report <- gsub("--datatable--", datatable.base, html.report)
html.report <- gsub("--datatable_css--", datatable.css.base, html.report)

html.report <- gsub("--function_profile_clusters--", all.cluster.functions, html.report)
html.report <- gsub("--button_clusters--", all.cluster.buttons, html.report)

all.motifs <- paste(paste("'", all.motifs, "'", sep = ""), collapse = ",")
html.report <- gsub("--all--", all.motifs, html.report)

## Export the report
html.report.file <- paste(basename, "_scan_profile_report.html", sep = "")
write(html.report, file = html.report.file)