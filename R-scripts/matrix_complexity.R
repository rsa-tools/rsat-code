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

#################################################################################
## Export a XY_plot with the distribution of values of cor and Ncor
## resulting from the comparison between the input motif and X permuted motifs
## the color of each point corresponds to its offset
Ncor_vs_cor_distribution <- function(motif.data){
  
  motif <- as.vector(motif.data$ID)
  
  ## Load the logo (to be displayed in ggplot)
  logo.file <- as.vector(motif.data$Logo)
  logo <- readPNG(logo.file)
  logo.roster <- rasterGrob(logo, interpolate = TRUE)
  
  ## Open the comparison table
  motif.comparison.tab <- read.csv(as.vector(motif.data$Comparison), sep = "\t", header = TRUE, comment.char = ";")
  # View(motif.comparison.tab)
  
  ## Get the abs value of the offset
  motif.comparison.tab$offset <- as.factor(abs(motif.comparison.tab$offset))
  
  ## Get the range of offset values
  offset.values <- sort(unique(motif.comparison.tab$offset), decreasing = TRUE)
  offset.classes <- length(offset.values)
  
  ## Calculate a color palette for the offset values
  myColors <- colorRampPalette(brewer.pal(9, "YlGnBu"), space="Lab")(offset.classes)
  names(myColors) <- offset.values
  
  cor.ncor.plot.file <- paste(as.vector(motif.data$Results_Folder), "/", motif, "cor_Ncor_plot", sep = "")
  
  ## Plot the cor vs Ncor values
  ggplot(data=motif.comparison.tab, aes(y = cor, x=Ncor, colour = offset)) +
    geom_point(size = 0.75, stroke = 0.75) + 
    ylim(c(0,1)) +
    xlim(c(0,1)) +
    geom_rug(position='jitter') +
    geom_abline(intercept = 0, color='steelblue', size=0.25, alpha=0.4) + 
    scale_colour_manual(name = "offset", values = myColors, labels = rev(paste(as.character(offset.values)))) +
    labs(title=paste(motif,": Distribution of comparison scores (cor and Ncor)"), y = "Column-wise Pearson Correlation (cor)", x="Width-Normalized Column-wise Pearson Correlation (Ncor)") + 
    annotation_custom(logo.roster, xmax = 1-0.05, xmin = 0.5-0.05, ymin = 0+0.01, ymax = 0.25-0.05)
  
  suppressMessages(ggsave(paste(cor.ncor.plot.file, ".pdf", sep = "")))
  suppressMessages(ggsave(paste(cor.ncor.plot.file, ".jpeg", sep = "")))
  
  ## Plot the Ncor distribution
  sub.data <- motif.comparison.tab[,c("cor", "Ncor")]
  melt.sub.data <- melt(sub.data)
  colnames(melt.sub.data) <- c("Metric", "Value")
  melt.sub.data$Value <- round(melt.sub.data$Value, digits = 2)
  
  ggplot(data = melt.sub.data, aes(x=Value)) +
    geom_histogram(aes(fill = Metric), stat ="bin", alpha=0.45, binwidth = 0.01)
  
  # ggplot(data = melt.sub.data, aes(x= Value)) + 
  # stat_ecdf(aes(colour = Metric)) +
  #   ylim(c(0,1)) +
  #   xlim(c(0,1)) 

  
  ggplot(data = melt.sub.data, aes(x= Metric, y = Value, fill = Metric)) +
  geom_violin(aes(color = Metric), trim = FALSE)  + 
  geom_boxplot(width = 0.2)

  
  ggplot(data = motif.comparison.tab, aes(x=Ncor, fill = Ncor)) +
    geom_bar(aes(fill = Ncor), stat ="bin", alpha=0.45, binwidth = 0.01)
  
  ##
  x.test <- round(motif.comparison.tab$Ncor, digits = 2)
  xx <- hist(x.test, breaks = seq(0,1,0.01), plot = FALSE)
  xx.freq <- round(xx$counts/sum(xx$counts), digits = 3)
  
  Ncor.df <- data.frame(Ncor = seq(0.01,1,0.01), Freq = xx.freq, Freq.class = freq.class)

  max.freq <- max(Ncor.df$Freq)
  
  ggplot(data= Ncor.df, aes(y = Freq, x = Ncor, fill = Freq)) +
    geom_bar(stat="identity") +
    scale_fill_distiller(palette = "RdPu", direction = 1) +
    annotation_custom(logo.roster, xmax = 0.25, xmin = 0.05, ymin = 0.05, ymax = (max.freq/3)) +
    labs(title=paste(motif, "- Distribution of Ncor values"), y = "Frequency", x="Ncor")

  sum.66 <- sum(xx.freq[0:66])
  
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
if (!exists("motif.attributes.table")) {
  stop("Missing mandatory argument (Motif attributes table): motif.attributes.table ")
} else if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
}

###################################
## Step 3: Read attributes table ##
###################################

setwd("/home/jaimicore/Desktop/matrix_complexity_demo")
# prefix <- "/home/jaimicore/Desktop/matrix-complexity_demo/matrix_complexity"
motif.attributes.table <- "matrix-complexity_test/matrix_complexity_tables/pairwise_compa_matrix_descriptions.tab"

motif.description <- read.table(motif.attributes.table, sep = "\t", header = FALSE)
colnames(motif.description) <- c("Motif_nb", "ID", "Name", "Length", "Nb_permutations", "Consensus", "Consensus_RC", "IC", "File", "File_permutations", "Results_Folder", "Comparison", "Logo", "Logo_RC")

## Get the motif IDs
motif.IDs <- as.vector(motif.description$ID)

thrash <- sapply(motif.IDs[3], function(m){
  
  ## Get the information of the current motif
  motif.data <- subset(motif.description, ID == m)
  
  Ncor_vs_cor_distribution(motif.data)
  
  Ncor_distribution_hist(motif.data)
  
})
rm(thrash)
