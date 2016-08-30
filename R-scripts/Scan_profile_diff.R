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
  head.tab <- "<div id='individual_motif_tab' style='width:2500px;display:none' class='tab div_chart_sp'><p style='font-size:12px;padding:0px;border:0px'><b>Individual Motif View</b></p><table id='Motif_tab' class='hover compact stripe' cellspacing='0' width='1190px' style='padding:15px;align:center;'><thead><tr><th class=\"tab_col\"> Motif_name </th><th class=\"tab_col\"> Motif_ID </th> <th class=\"tab_col\"> P-value </th> <th class=\"tab_col\"> E-value </th> <th class=\"tab_col\"> FDR </th> <th class=\"tab_col\"> Significance </th> <th class=\"tab_col\"> Nb of hits </th><th class=\"tab_col\"> Nb of sequences </th><th class=\"tab_col\">Fraction of sequences</th><th class=\"tab_col\"> Chi-squared</th><th class=\"tab_col\"> P-value </th> <th class=\"tab_col\"> E-value </th> <th class=\"tab_col\"> FDR </th> <th class=\"tab_col\"> Significance </th> <th class=\"tab_col\"> Nb of hits </th><th class=\"tab_col\"> Nb of sequences </th><th class=\"tab_col\">Fraction of sequences</th><th class=\"tab_col\"> Chi-squared</th><th class=\"tab_col\">Profile cluster</th><th class=\"tab_col\">DF</th> <th class=\"tab_col\"> Profile </th><th class=\"tab_col\"> Logo </th> <th class=\"tab_col\"> Logo (RC) </th></tr></thead><tbody>"
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

# else if (!exists("html.template.file")){
#   stop("Missing mandatory argument (HTML template to draw the profiles): html.template.file ")
# }

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

## Heatmap dendogram position
if (heatmap.dendo == "show"){
  heatmap.dendo <- "row"
} else if(heatmap.dendo == "hide"){
  heatmap.dendo <- "none"
}
print(heatmap.dendo)

## Create folder for individual profile plots
dir.create(paste(prefix, "_TFBSs_positional_profiles/", sep = ""), recursive = TRUE, showWarnings = FALSE )

# Borrar
# ## Create a file to store the resulting tables
# covered.tables.dir <- paste(prefix, "_covered_sequences_info", sep = "")
# dir.create(covered.tables.dir, showWarnings = FALSE)


#/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS_logos

# matrix.scan.file.query <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS_matrix_scan_results_PARSED.tab"
# matrix.scan.file.control <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS_matrix_scan_results_PARSED_control.tab"
# prefix <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS"
# ID.to.names.correspondence.tab <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS_TF_ID_name_correspondence.tab"
# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff")
# sequence.names.file.query <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS_matrix_scan_sequence_names.tab"
# sequence.names.file.control <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/Template/Diff/Pho_vs_SSPS_matrix_scan_sequence_names_control.tab"

####################################
## Read matrix-scan table (Query)
verbose(paste("Reading Query matrix-scan results table"), 1)
matrix.scan.results.query <- read.csv(file = matrix.scan.file.query, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results.query) <- c("seq_id", "ft_name", "bspos", "Pval")

######################################
## Read matrix-scan table (Control)
verbose(paste("Reading Control matrix-scan results table"), 1)
matrix.scan.results.control <- read.csv(file = matrix.scan.file.control, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results.control) <- c("seq_id", "ft_name", "bspos", "Pval")

#######################################
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

#############################
## Get the sequences + motif name's
seq.id.query <- unique(as.vector(matrix.scan.results.query$seq_id))
seq.id.control <- unique(as.vector(matrix.scan.results.control$seq_id))

matrix.names.query <- unique(as.vector(matrix.scan.results.query$ft_name))
matrix.names.control <- unique(as.vector(matrix.scan.results.control$ft_name))
matrix.names <- intersect(matrix.names.query, matrix.names.control)
verbose(paste(length(matrix.names), "motifs found in the intersection of Query vs Control sequences"), 1)

###################################
## Calculate the sequence limits
seq.length <- as.numeric(seq.length)

###################################
## Divide the sequences in bins 
# bin <- 50
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
matrix.scan.results.query$bspos <- matrix.scan.results.query$bspos + limits
matrix.scan.results.control$bspos <- matrix.scan.results.control$bspos + limits

ID.names.tab <- ID.to.names.correspondence.tab
ID.names <- read.table(ID.names.tab, sep = "\t")

feature.attributes <- vector("list", length(matrix.names)) 
profiles <- NULL
windows.labels <- NULL
setwd(results.folder)
print.formats <- c("pdf", "jpeg")

#########################################################################
## Create count table from matrix-scan results (if not exist in input) ##
#########################################################################

## Get the x-axis labels
xlab <- data.frame(windows)$start
xlab <- ifelse(xlab >= 0, xlab + bin, xlab)
# colnames(counts.per.bin.table) <- as.character(xlab)

input.count.table <- 0
seq.count.per.motif <- list()
all.counts.per.bin.query <- data.frame()
all.counts.per.bin.control <- data.frame()
query.over.control <- data.frame()
control.over.query <- data.frame()
feature.attributes <- vector("list", length(matrix.names))
feature.log2.ratio <- vector("list", length(matrix.names))

if(input.count.table == 0){
  
  verbose(paste("Creating counts and frequencies tables"), 1)
  thrash <- sapply(1:length(matrix.names), function(m){
    
    print(m)
    
    ## Select the matches of the query motif
    matrix.query <- matrix.names[m]
    
    matrix.query.selection.query.seq <- matrix.scan.results.query[matrix.scan.results.query$ft_name == matrix.query,]
    matrix.query.selection.control.seq <- matrix.scan.results.control[matrix.scan.results.control$ft_name == matrix.query,]
    
    ## Count the number of hits and the name of secuences with at least one hit
    nb.hits.query <- length(as.vector(matrix.query.selection.query.seq[matrix.query.selection.query.seq$ft_name == matrix.query,]$ft_name))
    nb.seq.query <- length(unique(as.vector(matrix.query.selection.query.seq[matrix.query.selection.query.seq$ft_name == matrix.query,]$seq_id)))
    fraction.query <- round(nb.seq.query/total.scanned.sequences.query, digits = 2)
    
    nb.hits.control <- length(as.vector(matrix.query.selection.control.seq[matrix.query.selection.control.seq$ft_name == matrix.query,]$ft_name))
    nb.seq.control <- length(unique(as.vector(matrix.query.selection.control.seq[matrix.query.selection.control.seq$ft_name == matrix.query,]$seq_id)))
    fraction.control <- round(nb.seq.control/total.scanned.sequences.control, digits = 2)

    
    ## To check
    if(off.set != 0){
      
      if (off.set.type == "start"){
        matrix.query.selection.query.seq$bspos <- matrix.query.selection.query.seq$bspos + seq.lengt - off.set
        matrix.query.selection.control.seq$bspos <- matrix.query.selection.control.seq$bspos + seq.lengt - off.set
      } else if (off.set.type == "end"){  
        matrix.query.selection$bspos <- matrix.query.selection$bspos + off.set
      }
      
    } else {
      matrix.query.selection.query.seq$bspos <- matrix.query.selection.query.seq$bspos# + limit.dw
      matrix.query.selection.control.seq$bspos <- matrix.query.selection.control.seq$bspos# + limit.dw
    }
    
    ## Convert the BSs in Ranges
    selection.IR.query <- IRanges(start = matrix.query.selection.query.seq$bspos, end = matrix.query.selection.query.seq$bspos)
    selection.IR.control <- IRanges(start = matrix.query.selection.control.seq$bspos, end = matrix.query.selection.control.seq$bspos)
    
    ## Count the overlap of BS in the bins
    counts.per.bin.query <- countOverlaps(windows, selection.IR.query)
    counts.per.bin.control <- countOverlaps(windows, selection.IR.control)
    
    ## 
    x1 <- round(counts.per.bin.query/sum(counts.per.bin.query), digits = 3)
    x2 <- round(counts.per.bin.control/sum(counts.per.bin.control), digits = 3)
    max.val <- max(c(x1,x2))
    
    query <- counts.per.bin.query
    control <- counts.per.bin.control
    names(query) <- xlab
    names(control) <- xlab
    
    #####################################
    ## Calculate the X2 in two senses:
    ## Query vs Control
    ## Control vs Query
    query.vs.control <- chisq.test(query, p = control/sum(control), correct = TRUE)
    control.vs.query <- chisq.test(control, p = query/sum(query), correct = TRUE)
    
    ## Get p-value and transform it in a pretty number
    query.vs.control.pval <- round(as.vector(query.vs.control[[3]]), digits = 100000000000000000)
    query.vs.control.pval.pretty <- prettyNum(query.vs.control.pval, scientific=TRUE, digits = 3)
    control.vs.query.pval <- round(as.vector(control.vs.query[[3]]), digits = 100000000000000000)
    control.vs.query.pval.pretty <- prettyNum(control.vs.query.pval, scientific=TRUE, digits = 3)
    
    ## Get e-value and transform it in a pretty number
    query.vs.control.eval <- round(as.vector(query.vs.control[[3]]), digits = 100000000000000000) * length(matrix.names)
    query.vs.control.eval.pretty <- prettyNum(query.vs.control.eval, scientific=TRUE, digits = 3)
    control.vs.query.eval <- round(as.vector(control.vs.query[[3]]), digits = 100000000000000000) * length(matrix.names)
    control.vs.query.eval.pretty <- prettyNum(control.vs.query.eval, scientific=TRUE, digits = 3)
    
    ## Get significances
    query.vs.control.sig <- round(-log10(query.vs.control.eval), digits = 2)
    control.vs.query.sig <- round(-log10(control.vs.query.eval), digits = 2)
    
    feature.attributes[[m]][["feature_id"]] <<- matrix.query
    feature.attributes[[m]][["Query_vs_Control_pval"]] <<- query.vs.control.pval.pretty 
    feature.attributes[[m]][["Query_vs_Control_eval"]] <<- query.vs.control.eval.pretty 
    feature.attributes[[m]][["Query_vs_Control_significance"]] <<- query.vs.control.sig
    feature.attributes[[m]][["Query_vs_Control_X2"]] <<- round(as.vector(query.vs.control[[1]]), digits = 2)
    feature.attributes[[m]][["Query_nb_hits"]] <<- nb.hits.query
    feature.attributes[[m]][["Query_nb_seq"]] <<- nb.seq.query
    feature.attributes[[m]][["Query_Fraction_of_sequences"]] <<- fraction.query
    feature.attributes[[m]][["Control_vs_Query_pval"]] <<- control.vs.query.pval.pretty
    feature.attributes[[m]][["Control_vs_Query_eval"]] <<- control.vs.query.eval.pretty 
    feature.attributes[[m]][["Control_vs_Query_significance"]] <<- control.vs.query.sig
    feature.attributes[[m]][["Control_vs_Query_X2"]] <<- round(as.vector(control.vs.query[[1]]), digits = 2)
    feature.attributes[[m]][["Control_nb_hits"]] <<- nb.hits.control
    feature.attributes[[m]][["Control_nb_seq"]] <<- nb.seq.control
    feature.attributes[[m]][["Control_Fraction_of_sequences"]] <<- fraction.control
    feature.attributes[[m]][["DF"]] <<- as.vector(query.vs.control[[2]])

    ########################################
    # ## Print the Profiles in JPEG and PDF
    for(pf in print.formats){

      if(pf == "pdf"){
        pdf.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", matrix.query, "_positional_profile.pdf", sep = "")
        pdf(pdf.file.name)
      } else {
        jpeg.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", matrix.query, "_positional_profile.jpeg", sep = "")
        jpeg(jpeg.file.name)
      }

      ## Draw the profile of the query sequences
      plot(y = query/sum(query),
           x = names(query),
           type = "l",
           col = "#00BFC4",
           lty = 1,
           lwd = 3,
           xlab = "Position (nt)",
           ylab = "Fraction of Binding Sites",
           main = paste(matrix.query, "Distribution of TFBSs\nQuery vs Control sequences"),
           panel.first=grid(col = "grey", lty = "solid"),
           ylim = c(0, max.val+0.05)
      )

      ## Draw the profile of the control sequences
      lines(y = control/sum(control),
            x = names(query),
            type = "l",
            col = "#F8766D",
            lty = 1,
            lwd = 3)

      ## Draw the legend
      legend("topleft", legend = c("Query (Q)", "Control (C)"), fill = c("#00BFC4", "#F8766D"), bty="o", bg="white")
      legend("bottomleft", legend = paste(c("E-value (Q vs C):", "E-value (C vs Q):"), c(query.vs.control.eval.pretty, control.vs.query.eval.pretty)), bty="o", bg="white")

      ## Add the logo in the plot
      logo.file <- paste(logo.folder, "/", matrix.query, "_logo.jpeg", sep = "")
      logo <- readJPEG(logo.file)
      rasterImage(logo,
                xleft = max(xlab) - 15 - round(seq.length/4.5),
                xright = max(xlab) - 15,
                ybottom = max.val,
                ytop = max.val+0.045)
      t <- dev.off()
    }
    
    ## Save the counts/fraction of the query and control sequences
    all.counts.per.bin.query <<- rbind(all.counts.per.bin.query, counts.per.bin.query)
    all.counts.per.bin.control <<- rbind(all.counts.per.bin.control, counts.per.bin.control)
    
    
    ## Save the log2 ratio between the query and control
    ## This will be used to plot the heatmap
    query.over.control <<- rbind(query.over.control, round(-log2(query/control), digits = 3))
    control.over.query <<- rbind(control.over.query, round(-log2(control/query), digits = 3))
    
  })
}
# counts.per.bin.table <- t(counts.per.bin)
#counts.per.bin.table <- round(counts.per.bin.table, digits = 3)
rm(thrash)
rownames(all.counts.per.bin.query) <- matrix.names
rownames(all.counts.per.bin.control) <- matrix.names
colnames(all.counts.per.bin.query) <- xlab
colnames(all.counts.per.bin.control) <- xlab

## Heatmap data rownames and colnames
rownames(query.over.control) <- matrix.names
rownames(control.over.query) <- matrix.names
colnames(query.over.control) <- xlab
colnames(control.over.query) <- xlab

all.fraction.query <- apply(all.counts.per.bin.query, 1, function(s){
  s/sum(s)
})
all.fraction.control <- apply(all.counts.per.bin.control, 1, function(s){
  s/sum(s)
})

##########################################
## Export Counts and Frequencies tables
query.counts.tab.file <- paste(basename, "_counts_per_bin_profiles_query.tab", sep = "") 
write.table(all.counts.per.bin.query, file = query.counts.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

control.counts.tab.file <- paste(basename, "_counts_per_bin_profiles_control.tab", sep = "") 
write.table(all.counts.per.bin.control, file = control.counts.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

fraction.tab.file.query <- paste(basename, "_tfbs_fraction_per_bin_profiles_query.tab", sep = "") 
write.table(all.fraction.query, file = fraction.tab.file.query, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

fraction.tab.file.control <- paste(basename, "_tfbs_fraction_per_bin_profiles_control.tab", sep = "") 
write.table(all.fraction.control, file = fraction.tab.file.control, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")

###############################
## Draw -log2 ratio heatmaps
## 1) Query vs Control
## 2) Control vs Query

## Create color palette
heatmap.color.classes <- 11
heatmap.color.palette <- "RdBu"
rgb.palette <- rev(colorRampPalette(brewer.pal(heatmap.color.classes, heatmap.color.palette), space="Lab")(heatmap.color.classes))

for(counter in 1:2){
  
  if(counter == 1){
    data.t <- query.over.control
    hm.main <- paste("Profile Heatmap\nQuery over Control")
  } else {
    hm.main <- paste("Profile Heatmap\nControl over Query")
    data.t <- control.over.query
  }
  
  for(pf in print.formats){

    if(pf == "pdf"){
      if(counter == 1){
        pdf.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", "Query_over_control_positional_profile.pdf", sep = "")
      } else {
        pdf.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", "Control_over_Query_positional_profile.pdf", sep = "")
      }
      pdf(pdf.file.name)
    } else {
      
      if(counter == 1){
        jpeg.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", "Query_over_control_positional_profile.jpeg", sep = "")
      } else {
        jpeg.file.name <- paste(basename(prefix), "_TFBSs_positional_profiles/", "Control_over_Query_positional_profile.jpeg", sep = "")
      }
      jpeg(jpeg.file.name)
    }
    
    data.t <- as.matrix(data.t)
    
    #############################
    ## Define profile clusters
    tree.profiles <- hclust(Dist(data.t, method = "correlation"), method = "complete")
    clusters.tree.profiles <- cutreeDynamic(tree.profiles, minClusterSize = 2, method = "tree")
    names(clusters.tree.profiles) <- tree.profiles[[4]]
    
    ## Generate a color palette
    nb.profile.clusters <- length(unique(clusters.tree.profiles))
    cluster.tree.profiles.palette <- colorRampPalette(brewer.pal(9, "Set1"), space="Lab")(nb.profile.clusters)
    
    ## Assign a different color to each cluster
    color.clusters.tree.profiles <- as.vector(sapply(clusters.tree.profiles+1, function(color){
      cluster.tree.profiles.palette[color]
    }))
  
    heatmap.2(data.t,
              
              main = hm.main,
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
              cexRow = 0.66,
              offsetCol = 0.25
    )
    t <- dev.off()
  }
}

##################################################
## Fraction of bound sequences Query vs Control
# x <- as.vector(feature.attributes.df$Query_fraction_of_sequences)
# y <- as.vector(feature.attributes.df$Control_fraction_of_sequences)
# plot(x,
#      y,
#      xlim = c(0,1),
#      ylim = c(0,1),
#      xlab = "Fraction of bound sequences (Query)",
#      ylab = "Fraction of bound sequences (Control)",
#      main = "Fraction of bound sequences\nQuery vs Control",
#      panel.first=grid(col = "grey", lty = "solid"),
#      col = "#00BFC4",
#      lty = 1, 
#      lwd = 3,
# )
# lines(x = c(0,1), y = c(0,1))

######################################
## Convert the list in a data frame
feature.attributes.df <- data.frame(t(
  matrix(as.vector(unlist(feature.attributes)), 
         ncol = length(feature.attributes))))
colnames(feature.attributes.df) <- c("Feature", "Pval_Q_vs_C", "Eval_Q_vs_C", "Significance_Q_vs_C", "Chi_Q_vs_C", "Query_nb_hits", "Query_nb_seq", "Query_fraction_of_sequences", "Pval_C_vs_Q", "Eval_C_vs_Q", "Significance_C_vs_Q", "Chi_C_vs_Q", "Control_nb_hits", "Control_nb_seq", "Control_fraction_of_sequences", "DF")

## Calculate q-values
## This step is executed once all the p-values were calculated
## The variable with class 'qvalue' is stored to its further exportation
p1 <- as.numeric(as.vector(feature.attributes.df$Pval_Q_vs_C))
features.qvalues.query.vs.control <- p.adjust(p1, method = "BH")
feature.attributes.df$Qval_Q_vs_C <- prettyNum(features.qvalues.query.vs.control, scientific=TRUE, digits = 2)

p2 <- as.numeric(as.vector(feature.attributes.df$Pval_C_vs_Q))
features.qvalues.control.vs.query <- p.adjust(p2, method = "BH")
feature.attributes.df$Qval_C_vs_Q <- prettyNum(features.qvalues.control.vs.query, scientific=TRUE, digits = 2)

#####################################
## Write the logo and profile path
TF.IDs <- as.vector(feature.attributes.df$Feature)
logos.F <- sapply(TF.IDs, function(i){
  paste(logo.folder, "/", i, "_logo.jpeg", sep = "")
})

logos.R <- sapply(TF.IDs, function(i){
  paste(logo.folder, "/", i, "_logo_rc.jpeg", sep = "")
})

## Write the Profile and TFBSs plots path
profiles.plots <- sapply(TF.IDs, function(i) {
  paste(basename(prefix), "_TFBSs_positional_profiles/", i, "_positional_profile.jpeg", sep = "")
})

###############################################
## Fill the Profile_cluster attribute
profile.clusters.names <- as.vector(sapply(matrix.names, function(m){
  profile.cluster <- as.vector(clusters.tree.profiles[m])
  paste("Profile_cluster_", profile.cluster, sep = "")
}))
profile.clusters.names.unique <- unique(profile.clusters.names)

############################
## Get the IDs of the TFs

## Set the motif names and IDs
TF.names <- as.vector(feature.attributes.df[,1])
TF.IDs <- as.vector(sapply(TF.names, function(x){
  as.vector(ID.names[which(ID.names[,2] == x),1])
}))

## As some IDs include a period ('.') in their text, we change it by
## an underscore ('_')
TF.IDs.cp <- gsub("\\.", "", TF.IDs)

feature.attributes.df$ID <- as.vector(TF.IDs.cp)
feature.attributes.df$Profile_cluster <- profile.clusters.names
feature.attributes.df$Profiles <- profiles.plots
feature.attributes.df$Logo <- logos.F
feature.attributes.df$Logo_RC <- logos.R

## Order the table according the Significance (-log10(E-value))
order.by.eval <- order(as.numeric(as.vector(feature.attributes.df$Eval_Q_vs_C)))
feature.attributes.df <- feature.attributes.df[order.by.eval,]
feature.attributes.df <- feature.attributes.df[,c(1,19,2,3,17,4,6,7,8,5,9,10,18,11,13,14,15,12,16,20,21,22,23)]

################################
## Create dynamic html report ##
################################

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

thrash <- sapply(order.by.eval, function(o){
  
  counter <<- counter + 1
  
  ## Get the counts and the fraction of TFBSs per bin in query and control
  counts.query <- all.counts.per.bin.query[o,]/sum(all.counts.per.bin.query[o,])
  counts.query <- round(counts.query, digits = 3)
  counts.control <- all.counts.per.bin.control[o,]/sum(all.counts.per.bin.control[o,])
  counts.control <- round(counts.control, digits = 3)
  
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
significance <- as.numeric(as.vector(feature.attributes.df$Significance_Q_vs_C))
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
      ID.names[which(ID.names[,2] == n),1]
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
# [1] "Feature"                       "Pval_Q_vs_C"                   "Eval_Q_vs_C"                   "Significance_Q_vs_C"          
# [5] "Chi_Q_vs_C"                    "Query_nb_hits"                 "Query_nb_seq"                  "Query_fraction_of_sequences"  
# [9] "Pval_C_vs_Q"                   "Eval_C_vs_Q"                   "Significance_C_vs_Q"           "Chi_C_vs_Q"                   
# [13] "Control_nb_hits"               "Control_nb_seq"                "Control_fraction_of_sequences" "DF"                           
# [17] "Qval_Q_vs_C"                   "Qval_C_vs_Q"                   "ID"                            "Profile_cluster"              
# [21] "Profiles"                      "Logo"                          "Logo_RC" 
profile.data.tab.html <- create.html.tab(feature.attributes.df, img = c(22,23), plot = c(21))
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
eval.rep <- repeat.n(as.vector(feature.attributes.df$Eval_Q_vs_C), times = 2)
evalues <- paste("evalues['", all.motifs, "'] = '", eval.rep, "';", sep = "")
evalues <- paste(evalues, collapse = "\n")
html.report <- gsub("--evalues--", evalues, html.report)

## Add the p-values (to display in the tooltip)
## They are inserted in the JS section
pval.rep <- repeat.n(as.vector(feature.attributes.df$Pval_Q_vs_C), times = 2)
pvalues <- paste("pvalues['", all.motifs, "'] = '", pval.rep, "';", sep = "")
pvalues <- paste(pvalues, collapse = "\n")
html.report <- gsub("--pvalues--", pvalues, html.report)

## Add the profile clusters (to display in the tooltip)
## They are inserted in the JS section
cl.rep <- repeat.n(as.vector(feature.attributes.df$Profile_cluster), times = 2)
profile.clusters.array <- paste(" profile_clusters['", all.motifs, "'] = '", cl.rep, "';", sep = "")
profile.clusters.array <- paste(profile.clusters.array, collapse = "\n")
html.report <- gsub("--profile_clusters_array--",  profile.clusters.array, html.report)

## Add the q-values (to display in the tooltip)
## They are inserted in the JS section
qval.rep <- repeat.n(as.vector(feature.attributes.df$Qval_Q_vs_C), times = 2)
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
  paste(logo.folder, i, "_logo.jpeg", sep = "")
})
logos.rep <- repeat.n(as.vector(feature.attributes.df$Logo), times = 2)
logos <- paste("pics['", all.motifs, "'] = '", logos.rep, "';", sep = "")
logos <- paste(logos, collapse = "\n")
html.report <- gsub("--pics--", logos, html.report)

## Logos in Reverse complement
logos.rc <- sapply(TF.IDs, function(i){
  paste(logo.folder, i, "_logo_rc.jpeg", sep = "")
})
logos.rc.rep <- repeat.n(as.vector(feature.attributes.df$Logo_RC), times = 2)
logos.rc <- paste("pics_rc['", all.motifs, "'] = '", logos.rc.rep, "';", sep = "")
logos.rc <- paste(logos.rc, collapse = "\n")
html.report <- gsub("--pics_rc--", logos.rc, html.report)

## Add the signficance (to display in the tooltip)
## They are inserted in the JS section
sig.rep <- repeat.n(as.vector(feature.attributes.df$Significance_Q_vs_C), times = 2)
sig <- paste("significances['", all.motifs, "'] = ", sig.rep, ";", sep = "")
sig <- paste(sig, collapse = "\n")
html.report <- gsub("--significances--", sig, html.report)

## Add the covertures (to display in the tooltip)
## They are inserted in the JS section
cc <- as.numeric(gsub("%", "", feature.attributes.df$Query_fraction_of_sequences))
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
max.fraction <- round(max(all.fraction.query, all.fraction.control), digits = 2)

max.y <- max.fraction + 0.02
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
html.report <- gsub("--heatmap_png--", paste(basename, "_profiles_heatmap.jpg", sep = ""), html.report)
html.report <- gsub("--heatmap_pdf--", paste(basename, "_profiles_heatmap.pdf", sep = ""), html.report)

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
min.sig <- min(as.numeric(as.vector(c(feature.attributes.df$Significance_Q_vs_C, feature.attributes.df$Significance_C_vs_Q))))
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