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
create.html.tab <- function(tab, img = 0, plot = 0){
  
  full.tab <- NULL
  head.tab <- "<div id='Motif_tab' style='width:1600px;display:none' class='tab div_chart_sp'><p style='font-size:12px;padding:0px;border:0px'><b>Individual Motif View</b></p><table id='Motif_dyn_table' class='hover compact stripe' cellspacing='0' width='1190px' style='padding:15px;align:center;'><thead><tr><th class=\"tab_col\"> Motif_ID </th><th class=\"tab_col\"> Motif_name </th> <th class=\"tab_col\"> Consensus </th> <th class=\"tab_col\"> Consensus (RC) </th> <th class=\"tab_col\"> Motif length </th> <th class=\"tab_col\"> Nb permutations </th> <th class=\"tab_col\"> IC </th> <th class=\"tab_col\"> Logo </th> <th class=\"tab_col\"> Logo (RC) </th> <th class=\"tab_col\"> cor::Ncor scatterplot </th> <th class=\"tab_col\"> cor::Ncor histogram </th> <th class=\"tab_col\"> cor::Ncor violin </th> <th class=\"tab_col\"> Ncor distrib </th> <th class=\"tab_col\"> cor distrib </th> </tr></thead><tbody>"
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
    paste(row.head, rows.text, rows.pic.text, rows.plot.text, row.tail, sep = "")    
    
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
  
  ###############################################
  ## Load the logo (to be displayed in ggplot) ##
  ###############################################
  logo.file <- as.vector(motif.data$Logo)
  logo <- readPNG(logo.file)
  logo.roster <- rasterGrob(logo, interpolate = TRUE)
  
  ## Open the comparison table
  motif.comparison.tab <- read.csv(as.vector(motif.data$Comparison), sep = "\t", header = TRUE, comment.char = ";")
  # View(motif.comparison.tab)
  
  ## Get the abs value of the offset
  motif.comparison.tab$offset <- as.factor(abs(motif.comparison.tab$offset))
  
  ########################################
  ## Scatterplot for cor vs Ncor values ##
  ########################################
  ## Get the range of offset values
  offset.values <- sort(unique(motif.comparison.tab$offset), decreasing = TRUE)
  offset.classes <- length(offset.values)
  
  ## Calculate a color palette for the offset values
  myColors <- colorRampPalette(brewer.pal(9, "YlGnBu"), space="Lab")(offset.classes)
  names(myColors) <- offset.values
  
  cor.ncor.plot.file <- paste(as.vector(motif.data$Results_Folder), "/", motif, "_cor_Ncor_plot", sep = "")
  
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
  
  #####################################
  ## Plot the Ncor and cor histogram ##
  #####################################
  cor.ncor.plot.hist <- paste(as.vector(motif.data$Results_Folder), "/", motif, "_cor_Ncor_hist", sep = "")
  
  sub.data <- motif.comparison.tab[,c("cor", "Ncor")]
  melt.sub.data <- melt(sub.data)
  colnames(melt.sub.data) <- c("Metric", "Value")
  melt.sub.data$Value <- round(melt.sub.data$Value, digits = 2)
  
  ## Calculate the counts per bin
  max.cor.counts <- max(hist(sub.data$cor, breaks = seq(0,1,0.01), plot = FALSE)$counts)
  max.Ncor.counts <- max(hist(sub.data$Ncor, breaks = seq(0,1,0.01), plot = FALSE)$counts)
  max.count <- max(c(max.cor.counts, max.Ncor.counts))

  ggplot(data = melt.sub.data, aes(x=Value)) +
    geom_area(aes(fill = Metric), stat ="bin", alpha=0.5, binwidth = 0.01, position="identity") +
    annotation_custom(logo.roster, xmax = 0.95, xmin = 0.75, ymin = max.count - (max.count/4), ymax = max.count) +
    labs(title=paste(motif, "- cor and Ncor histogram"), y = "Frequency", x="Value (cor and Ncor)")

  suppressMessages(ggsave(paste(cor.ncor.plot.hist, ".pdf", sep = "")))
  suppressMessages(ggsave(paste(cor.ncor.plot.hist, ".jpeg", sep = "")))

  ############################################
  ## Plot the Ncor and cor Violin + boxplot ##
  ############################################
  cor.ncor.plot.vioplot <- paste(as.vector(motif.data$Results_Folder), "/", motif, "_cor_Ncor_violin_plot", sep = "")

  ggplot(data = melt.sub.data, aes(x= Metric, y = Value)) +
    geom_violin(aes(fill = Metric), trim = FALSE)  + 
    geom_boxplot(width = 0.2) +
    labs(title=paste(motif, "- cor and Ncor distributions"), y = "Values")
  
  suppressMessages(ggsave(paste(cor.ncor.plot.vioplot, ".pdf", sep = "")))
  suppressMessages(ggsave(paste(cor.ncor.plot.vioplot, ".jpeg", sep = "")))

  ################################################
  ## Plot the Ncor and cor individual histogram ##
  ################################################
  for(metric in c("cor", "Ncor")){
    
    metric.indiv.hist <- paste(as.vector(motif.data$Results_Folder), "/", motif, "_", metric, "_distribution", sep = "")
    
    values <- round(motif.comparison.tab[,metric], digits = 2)
    values.hist <- hist(values, breaks = seq(0,1,0.01), plot = FALSE)
    freq <- round(values.hist$counts/sum(values.hist$counts), digits = 3)
    
    metric.mean <- round(mean(values), digits = 2)
    metric.median <- round(median(values), digits = 2)
    
    metric.df <- data.frame(Metric = seq(0.01,1,0.01), Freq = freq)
    max.freq <- max(metric.df$Freq)
    
    ggplot(data= metric.df, aes(y = Freq, x = Metric, fill = Freq)) +
      geom_bar(stat="identity") +
      scale_fill_distiller(palette = "RdPu", direction = 1) +
      annotation_custom(logo.roster, xmax = 0.25, xmin = 0.05, ymin = 0.05, ymax = (max.freq/3)) +
      labs(title=paste(motif, "- Distribution of", metric, "values"), y = "Frequency", x=metric) +
      geom_vline(aes(xintercept = metric.mean, linetype = "Mean")) +
      geom_vline(aes(xintercept = metric.median, linetype = "Median")) +
      scale_linetype_manual(name = "", labels = c("Mean", "Median"), values = c("Mean" = 1, "Median" = 2))
    
    suppressMessages(ggsave(paste(metric.indiv.hist, ".pdf", sep = "")))
    suppressMessages(ggsave(paste(metric.indiv.hist, ".jpeg", sep = ""))) 
  }
  
  return(c(paste(cor.ncor.plot.file, ".jpeg", sep = ""), 
    paste(cor.ncor.plot.file, ".pdf", sep = ""),
    paste(cor.ncor.plot.hist, ".jpeg", sep = ""),
    paste(cor.ncor.plot.hist, ".pdf", sep = ""),
    paste(cor.ncor.plot.vioplot, ".jpeg", sep = ""),
    paste(cor.ncor.plot.vioplot, ".pdf", sep = ""),
    paste(as.vector(motif.data$Results_Folder), "/", motif, "_Ncor_distribution.jpeg", sep = ""),
    paste(as.vector(motif.data$Results_Folder), "/", motif, "_Ncor_distribution.pdf", sep = ""),
    paste(as.vector(motif.data$Results_Folder), "/", motif, "_cor_distribution.jpeg", sep = ""),
    paste(as.vector(motif.data$Results_Folder), "/", motif, "_cor_distribution.pdf", sep = "")
  ))

  # sum.66.cor <- sum(xx.freq[0:66])
}

#######################################
## Generate and fill the HTML report
generate.html.report <- function(){
  
  html.report <- readLines(html.template.file)
  verbose(paste("Fill html report"), 1)
  
  ###########################
  ## Fill parameters table ##
  ###########################
  html.report <- gsub("--motif_nb--", nb.motifs, html.report)
  
  ############################################
  ## Create and fill individual motif table ##
  ############################################
  # [1] "Motif_nb"          "ID"                "Name"             
  # [4] "Length"            "Nb_permutations"   "Consensus"        
  # [7] "Consensus_RC"      "IC"                "File"             
  # [10] "File_permutations" "Results_Folder"    "Comparison"       
  # [13] "Logo"              "Logo_RC"           "scatterplot_jpeg" 
  # [16] "scatterplot_pdf"   "hist_jpeg"         "hist_pdf"         
  # [19] "violin_jpeg"       "violin_jpeg"       "Ncor_dist_jpeg"   
  # [22] "Ncor_dist_pdf"     "cor_dist_jpeg"     "cor_dist_pdf"
  motif.tab.html <- create.html.tab(motif.description[,c(2,3,6,7,4,5,8,13,14,15,17,19,21,23)], img = c(8,9), plot = c(10,11,12,13,14))
  motif.tab.html <- paste(motif.tab.html, collapse = "\n")
  html.report <- gsub("--tab--", motif.tab.html, html.report)
  
  verbose(paste("Exporting html report"), 1)
  html.report.file <- paste(basename(prefix), "_scan_profile_report.html", sep = "")
  write(html.report, file = html.report.file)
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
} else if(!exists("html.template.file")){
  stop("Missing mandatory argument (html.template.file): html.template.file ")
}


###################################
## Step 3: Read attributes table ##
###################################

# setwd("/home/jaimicore/Desktop/matrix_complexity_demo")
# prefix <- "/home/jaimicore/Desktop/matrix_complexity_demo/matrix-complexity_demo/matrix_complexity"
# motif.attributes.table <- "/home/jaimicore/Desktop/matrix_complexity_demo/matrix-complexity_demo/matrix_complexity_tables/pairwise_compa_matrix_descriptions.tab"

motif.description <- read.table(motif.attributes.table, sep = "\t", header = FALSE)
colnames(motif.description) <- c("Motif_nb", "ID", "Name", "Length", "Nb_permutations", "Consensus", "Consensus_RC", "IC", "File", "File_permutations", "Results_Folder", "Comparison", "Logo", "Logo_RC")

setwd(dirname(prefix))

## Get the motif IDs
motif.IDs <- as.vector(motif.description$ID)
nb.motifs <- nrow(motif.description)

## Iterate over the motif IDs
exported.files <- sapply(motif.IDs, function(m){
  
  ## Get the information of the current motif
  motif.data <- subset(motif.description, ID == m)
  
  ## Generate the plots
  Ncor_vs_cor_distribution(motif.data)
})
exported.files.df <- as.data.frame(t(exported.files))

names(exported.files.df) <- c("scatterplot_jpeg", "scatterplot_pdf", "hist_jpeg", "hist_pdf", "violin_jpeg", "violin_jpeg", "Ncor_dist_jpeg", "Ncor_dist_pdf", "cor_dist_jpeg", "cor_dist_pdf")

## Complete the description table
motif.description <- cbind(motif.description, exported.files.df) 

## Export the HTML file
verbose(paste("Creating HTML dynamic report"), 1)
generate.html.report()

