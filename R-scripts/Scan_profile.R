## Required libraries
suppressPackageStartupMessages(library("IRanges", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("RColorBrewer", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("gplots", warn.conflicts=FALSE, character.only = TRUE))
suppressPackageStartupMessages(library("jpeg", warn.conflicts=FALSE, character.only = TRUE))

#################################################################################################
## Functions

create.html.tab <- function(tab, img = 0){
  
  full.tab <- NULL
  head.tab <- "<div id='individual_motif_tab' style='width:1200px;display:none' class='tab div_chart_sp'><p style='font-size:15px;padding:0px;border:0px'><b>Individual Motif View</b></p><table id='Motif_tab' class='hover compact stripe' cellspacing='0' width='1190px' style='padding:15px;align:center;'><thead><tr><th class=\"tab_col\"> Motif_name </th><th class=\"tab_col\"> Motif_ID </th> <th class=\"tab_col\"> P-value </th> <th class=\"tab_col\"> E-value </th> <th class=\"tab_col\"> Significance </th> <th class=\"tab_col\"> Seq with hits</th> <th class=\"tab_col\"> Nb of hits </th> <th class=\"tab_col\"> Logo </th></tr></thead><tbody>"
  
  
  content.tab <- apply(tab, 1, function(row){
    
    row.length <- length(row)
    rows.nb <- 1:row.length
    if(img != 0){
      
      ## Get the number of the columns with/without picture
      ## This is done becuase the tab require different arguments
      rows.no.pic <- rows.nb[!(rows.nb %in% img)]
      rows.pic <- rows.nb[rows.nb %in% img]
      
      row.head <- "<tr>"
      rows.no.pic.text <- paste("<td>", row[rows.no.pic], "</td>", collapse = "")
      rows.pic.text <- paste("<td><img class='logo_tab' src ='", row[rows.pic], "'/></td>", collapse = "")
      row.tail <- "</tr>"
      
      paste(row.head, rows.no.pic.text, rows.pic.text, row.tail, sep = "")
      
    } else{
      paste("<tr>", paste("<td>", row, "</td>", collapse = ""), "</tr>",sep = "")
    }
  })
  
  
  tail.tab <- "</tbody></table></div>"
  
  full.tab <- c(head.tab, content.tab , tail.tab)
  return(full.tab)
  
}

########################################################################################


###########################################
## Read arguments from the command line.
##
## Arguments passed on the command line
## will over-write the default arguments
## specified above.
message("Reading arguments from command-line")
args <- commandArgs(trailingOnly=TRUE);
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

message("Checking mandatory arguments")
if (!exists("matrix.scan.file")) {
  stop("Missing mandatory argument (matrix-scan results table): matrix.scan.file ")
} else if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
} else if (!exists("ID.to.names.correspondence.tab")) {
  stop("Missing mandatory argument (Correspondence of the motif IDs and names): ID.to.names.correspondence.tab ")
} else if (!exists("seq.length")) {
  
  ## For the moment we assume all the input sequences have the same size
  stop("Missing mandatory argument (Sequence length): seq.length ")
}

if (!exists("p.val")) {
  p.val <- 1e-3
} 
if (!exists("bin")) {
  bin <- 50
}
if (!exists("draw.area")) {
  draw.area <- 0
}
if (!exists("individual.plots")) {
  individual.plots <- 0
}

# setwd("/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/")
# matrix.scan.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/matrix_scan_pval_1e-3_GAF_OnTheFly_bg_mkv_2_PARSED.tab"
# matrix.scan.file <- "/home/jaimicore/Documents/PhD/Human_promoters_project/Drosophila_TFs_MArianne/Bin/matrix_scan_pval_1e-3_G_mkv_2_PARSED.tab"

#############################################
## Read matrix-scan table Active Promoters
message("Reading matrix-scan table")
matrix.scan.results <- read.csv(file = matrix.scan.file, sep = "\t", header = TRUE, comment.char = ";")
colnames(matrix.scan.results) <- c("seq_id", "ft_name", "bspos", "Pval")

#################
## Set p-value
p.val <- as.numeric(p.val)

#############################
## Get the matrices name's
matrix.names <- unique(as.vector(matrix.scan.results$ft_name))

#############################
## Get the sequences + motif name's
seq.id <- unique(as.vector(matrix.scan.results$seq_id))
matrix.names <- unique(as.vector(matrix.scan.results$ft_name))

###################################
## Calculate the sequence limits
seq.length <- as.numeric(seq.length)
limits <- seq.length/2

###################################
## Dive the sequences in bins 
# bin <- 50
message(paste("Setting bins of size", bin))
bin <- as.numeric(bin)
windows <- IRanges(start = seq(from = -limits, to = limits - bin + 1, by = bin), width = bin)


profiles <- NULL
all.pvalues <- NULL
all.evalues <- NULL
windows.labels <- NULL
all.nb.seq <- NULL
all.nb.hits <- NULL
all.max.pval.hits <- NULL
all.min.pval.hits <- NULL
all.chisq <- NULL
all.df <- NULL

ID.to.names.correspondence.tab <- "ID_names.txt"
ID.names.tab <- ID.to.names.correspondence.tab
ID.names <- read.table(ID.names.tab, sep = "\t")

# prefix <- paste("Jaspar_insects", "_pval_", p.val, sep = "")
message("Calculating P-value. Null Hypothesis: homogenous distribution of the TFBSs")
if(individual.plots == 1){
  message("Print the profiles in a PDF file")
}


pdf.file.name <- paste(prefix, "_positional_profiles.pdf", sep = "")
pdf(pdf.file.name)

##################################################
## Normalize the number of hits before plotting ##
##################################################
sapply(1:length(matrix.names), function(m){

#   print(m)
  
  ## Select the matches of the query motif
  matrix.query <- matrix.names[m]
  matrix.query.selection <- matrix.scan.results[matrix.scan.results$ft_name == matrix.query,]
  
  ## Get the number of sequences with at least one hit
  nb.seq.with.hits <- length(unique(as.vector(matrix.query.selection$seq_id)))
  all.nb.seq <<- append(all.nb.seq, nb.seq.with.hits)
  
  ## Get the total number of hits for each motif
  total.hits <- dim(matrix.query.selection)[1]
  all.nb.hits <<- append(all.nb.hits, total.hits)
  
  ## Get the max p-value among the p-values of the matches
  max.pval.hit <- max(matrix.query.selection$Pval)
  all.max.pval.hits <<- append(all.max.pval.hits, max.pval.hit)
  
  ## Get the min p-value among the p-values of the matches
  min.pval.hit <- min(matrix.query.selection$Pval)
  all.min.pval.hits <<- append(all.min.pval.hits, min.pval.hit)
  
  ## As the reference point in matrix-scan was the end of the sequence and as we are working with peaks
  ## we add 300 to the position to have -/+ position around the summit
  matrix.query.selection$bspos <- matrix.query.selection$bspos + 300
  
  ## Convert the BSs in Ranges
  selection.IR <- IRanges(start = matrix.query.selection$bspos, end = matrix.query.selection$bspos)
  
  ## Count the overlap of BS in the bins
  counts.per.bin <- countOverlaps(windows, selection.IR)
  
  ## Chi-square test to calculate a p-value
  ## H0 = the number of hits (TFBSs) follow an homogeneous distribution
  count.per.bin.HO <- round(sum(counts.per.bin)/length(counts.per.bin))
  HO.expected.counts <- rep(count.per.bin.HO, times = length(counts.per.bin))
  
  ## Goodnes-of-fit X2 test
  ## We don't require a probability vector because we assume the TFBSs (counts) 
  ## are distributed homogenously along the sequences
  chi <- chisq.test(counts.per.bin, correct = TRUE)
  
  ## Chi-squared
  cs.val <- chi[[1]]
  all.chisq <<- append(all.chisq, cs.val)
  
  ## Degrees of freedom
  df <- chi[[2]]
  all.df <<- append(all.df, df)
  
  ## Get p-value
  chi.pval <- round(as.numeric(chi[[3]]), digits = 100000)
  
  ## Caclculate e-value
  chi.eval <- chi.pval * length(matrix.names)
  
  chi.pval <- prettyNum(chi.pval, scientific=TRUE, digits = 2)
  chi.eval <- prettyNum(chi.eval, scientific=TRUE, digits = 2)

  windows.df <- data.frame(windows)

  names(counts.per.bin) <- windows.df$start
  
  profiles <<- cbind(profiles,(counts.per.bin/sum(counts.per.bin)))
  y.val <- counts.per.bin/sum(counts.per.bin)
  x.val <- as.numeric(names(counts.per.bin))

  all.pvalues <<- append(all.pvalues, chi.pval)
  all.evalues <<- append(all.evalues, chi.eval)
  
  windows.labels <<- names(counts.per.bin)
  
  if(individual.plots == 1){
    plot(x = c(-300,-300,50,50),
         y = c(-0.1,-0.4,-0.4,-0.1),
         type = "l",
  #        ylim = c(0, max(c(counts.per.range.inactive, counts.per.range.inactive))+50),
          ylim = c(0, 0.20),
         xlim = c(-300, 300),
         col = "#ffeda0",
         lwd = 1,
         ## Labels
         main = paste("Motif:", matrix.query),
         xlab = "Distance to peak summit",
         ylab = "Normalized Nb Hits",
         ## Hide x-axis
         xaxt='n', 
    )  
  
  # polygon(x = c(-250,-250, 50, 50), y = c(0, 1, 1, 0), col="#ffeda0", border = NA, lty = 0, )
      
    ## Draw the grid
    abline(v=(x.val), col="lightgray", lty="dotted")
    abline(h=(seq(from = 0, to = 1, by = 0.01)), col="lightgray", lty="dotted")
  
    ## Draw the TSS (+1) position
    abline(v = 0, col="#045a8d", lwd = 2, lty = 2)
    
    ## Set x-axis values 
    axis(side = c(1,2,3,4), at = as.character(x.val), labels = as.character(x.val))
  
    ## Draw the lines for the active promoters 
    lines(x = x.val, y = y.val, type = "l", col = "#00BFC4", lty = 1, lwd = 3)
  
  
    ## Draw the legend
     legend("topleft", legend = c(paste(matrix.query, "profile"), "Peak summit"), fill = c("#00BFC4", "#045a8d"), bty="o", bg="white")
  #   legend("topright", legend = paste("E-value:", TFBSs.enrichment.eval), bty="o", bg="white")
  
    matrix.ID <- as.vector(ID.names[which(ID.names[,2] == matrix.query),1])
    logo.file <- paste("logos/", matrix.ID, "_logo.jpeg", sep = "")
    logo <- readJPEG(logo.file)
    rasterImage(logo, 
              xleft = 60,
              xright = 275, 
              ybottom = 0.14,
              ytop = 0.195)
  }
})
if(individual.plots == 1){
  dev.off()
}


# profiles <- data.frame(t(profiles))
rownames(profiles) <- as.character(data.frame(windows)$start)
colnames(profiles) <- matrix.names
profiles <- as.matrix(profiles)
profiles <- round(profiles, digits = 3)


#######################################
## Create the dynamic report in HTML ##
#######################################

# aa <- profiles
#  profiles <- aa

## Set the motif names and IDs
TF.names <- colnames(profiles)

## Ordert the TF.names (and other variables) according the Significance (-log10(E-value))
order.by.eval <- rev(order(-log(as.numeric(all.evalues))))
TF.names <- TF.names[order.by.eval]
profiles <- profiles[,order.by.eval]
evalues.vector <- sort(as.numeric(all.evalues))
pvalues.vector <- sort(as.numeric(all.pvalues))
all.nb.seq <- all.nb.seq[order.by.eval]
all.nb.hits <- all.nb.hits[order.by.eval]
all.pval.match <- rep(p.val, times = length(TF.names))
all.max.pval.hits <- all.max.pval.hits[order.by.eval]
all.min.pval.hits <- all.min.pval.hits[order.by.eval]
all.chisq <- all.chisq[order.by.eval]
all.df <- all.df[order.by.eval]

## Get the IDs of the TF
TF.IDs <- as.vector(sapply(TF.names, function(x){
  as.vector(ID.names[which(ID.names[,2] == x),1])
}))
## As each motif ID is unique, we used it also in the html file,
## However as some IDs include a period ('.') in their text, we change it by
## an underscore ('_')
TF.IDs.cp <- gsub("\\.", "", TF.IDs)

## Set colors
set.colors <- colorRampPalette(brewer.pal(10,"Paired"))(length(TF.IDs))
# set.colors <- colorRampPalette(brewer.pal(9,"YlOrRd"))(length(TF.IDs))

counter <- 0
x.correspondence <- NULL
x.y <- NULL
plot.names <- NULL
area <- NULL
all.motifs <- NULL
all.motif.names <- NULL

## Each column of the variable profiles correspond to the counts per bin of each motif
thrash <- apply(profiles, 2, function(values){
  
  counter <<- counter + 1
  
  ## Here we create a unique ID without CSS speacial characters
#   motif <- paste(TF.IDs[counter], "_", TF.names[counter], sep ="")
  motif <- paste(TF.IDs.cp[counter], counter, sep = "")  
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
})


## Set the line width according the significance -log10(E-value)
## Higher significance means a wider line
significance <- rev(sort(-log(as.numeric(all.evalues))))
significance <- round(significance, digits = 2)
q <- quantile(significance)

line.w <- sapply(significance, function(s){
  if(s <= q[2]){
    w <- 1
  } else if (s >= q[2] & s <= q[3]){
    w <- 2
  } else if (s >= q[3] & s <= q[4]){
    w <- 3
  } else if (s >= q[4]){
    w <- 5
  }
})


## Write the logo's path
logos.F <- sapply(TF.IDs, function(i){
  paste("logos/",i, "_logo.jpeg", sep = "")
})

## Temporary not available
logos.R <- sapply(TF.IDs, function(i){
  paste("logos/",i, "_logo_rc.jpeg", sep = "")
})

## Create a Dataframe containing the information of all motifs
## This table will be exported and displayed as a dynamic table in the report
profile.data.tab <- data.frame(all.motif.names, TF.IDs, pvalues.vector, evalues.vector, significance, logos.F, all.nb.seq, all.nb.hits, all.pval.match, all.max.pval.hits, all.min.pval.hits, all.chisq, all.df)
colnames(profile.data.tab) <- c("Motif_name", "Motif_ID", "P-val", "E-val", "Significance", "Logo", "Nb_seq_with_hits", "Total_nb_hits", "P_val_threshold", "Max_pval", "Min_pval", "Chi-squared", "Degrees_Freedom")
rownames(profile.data.tab) <- NULL
profile.data.tab.file <- paste(prefix, "_profiles.tab", sep = "")
write.table(profile.data.tab, file = profile.data.tab.file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


##########################
## Draw Profiles heatmap
message("Drawing Heatmap profiles")
profiles.tab.file <- paste(prefix, "_counts_per_bin_profiles.tab", sep = "") 
write.table(t(profiles), file = profiles.tab.file, quote = FALSE, col.names = TRUE, row.names = TRUE)

## Color palette
rgb.palette <- colorRampPalette(brewer.pal(11, "RdBu"), space="Lab")

## Heatmap
out.format <- c("pdf", "jpg")
for (format in out.format){

  profiles.heatmap.file <- paste(prefix, "_profiles_heatmap.", format, sep = "") 
  
  if(format == "pdf"){
    pdf(profiles.heatmap.file)
  } else if (format == "jpg"){
    jpeg(profiles.heatmap.file)
  }
  
  heatmap.2(t(profiles),
            
            ## Dendrogram control
            dendrogram = c("none"),
            Rowv = TRUE,
            Colv = FALSE,
            
            main = "Profile Heatmap",
            xlab = "Position (bp)",
            ylab = "Motifs",
            
            #            hclustfun = function(d){hclust(d, method="ward")},
            
            ## Color
            col = rgb.palette,
            
            ## Trace
            trace = "none",
            
            ## Key control
            key = TRUE,
            keysize = 1,
            density.info = "none",
            key.xlab = "Density",
            key.ylab = "",
            key.title = "",
            offsetCol = 0.25,
            cexRow = 0.25,
  )
  dev.off()
  
}



## Fill the HTML template
## Substitute the words marked in the tamplate by the data
html.template.file <- "Template/index.html"
html.report <- readLines(html.template.file)

profile.data.tab.html <- create.html.tab(profile.data.tab[,c(1:8)], img = 6)
profile.data.tab.html <- paste(profile.data.tab.html, collapse = "\n")
html.report <- gsub("--tab--", profile.data.tab.html, html.report)

x.y <<- rbind(x.y, paste("['x',", paste(rownames(profiles), collapse = ","), "],", sep = ""))

## CSS section to set the line width
## Note: the width is proportional to the significance
line.w <- paste("#chart .c3-line-", all.motifs, "{ stroke-width: ", line.w, "px; }", sep = "")
line.w <- paste(line.w, collapse = "\n")
html.report <- gsub("--lines_w--", line.w, html.report)

## Add the e-values data
## They are inserted in the JS section
evalues <- paste("evalues['", all.motifs, "'] = '", evalues.vector, "';", sep = "")
evalues <- paste(evalues, collapse = "\n")
html.report <- gsub("--evalues--", evalues, html.report)

## Add the p-values (to display in the tooltip)
## They are inserted in the JS section
pvalues <- paste("pvalues['", all.motifs, "'] = '", pvalues.vector, "';", sep = "")
pvalues <- paste(pvalues, collapse = "\n")
html.report <- gsub("--pvalues--", pvalues, html.report)

## Add the real motif IDs (to display in the tooltip)
## They are inserted in the JS section
## I called them 'real' because are those found on the original motif file
IDs <- paste("IDs['", all.motifs, "'] = '", TF.IDs, "';", sep = "")
IDs <- paste(IDs, collapse = "\n")
html.report <- gsub("--IDs--", IDs, html.report)

## Add the real motif logo path (to display in the tooltip)
## They are inserted in the JS section
logos <- sapply(TF.IDs, function(i){
  paste("logos/",i, "_logo.jpeg", sep = "")
})
logos <- paste("pics['", all.motifs, "'] = '", logos, "';", sep = "")
logos <- paste(logos, collapse = "\n")
html.report <- gsub("--pics--", logos, html.report)

## Add the signficance (to display in the tooltip)
## They are inserted in the JS section
sig <- paste("significances['", all.motifs, "'] = '", significance, "';", sep = "")
sig <- paste(sig, collapse = "\n")
html.report <- gsub("--significances--", sig, html.report)

## Add x values (one row per motif)
## They are inserted in the C3 section
xx <- paste(x.y, collapse = "\n")
html.report <- gsub("--x_y--", xx, html.report)

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


## Insert the motif names (to hide/show all)
## They are inserted in the JQuery section
all.motifs <- paste(paste("'", all.motifs, "'", sep = ""), collapse = ",")
html.report <- gsub("--all--", all.motifs, html.report)

## Insert the Y axis limits
## They are inserted in the C3section
max.y <- max(profiles) + 0.02
html.report <- gsub("--y_axis--", max.y, html.report)

## Fill the parameters table
html.report <- gsub("--bin_l--", bin, html.report)
html.report <- gsub("--bin_nb--", length(windows.labels), html.report)
html.report <- gsub("--seq_nb--", length(seq.id), html.report)
html.report <- gsub("--motif_nb--", length(matrix.names), html.report)
html.report <- gsub("--p--", prettyNum(p.val), html.report)

## Fill the heatmap section
html.report <- gsub("--heatmap_png--", paste(prefix, "_profiles_heatmap.jpg", sep = ""), html.report)
html.report <- gsub("--heatmap_pdf--", paste(prefix, "_profiles_heatmap.pdf", sep = ""), html.report)

## Export the report
html.report.file <- paste(prefix, "_scan_profile_report.html", sep = "")
write(html.report, file = html.report.file)


# grep -v '^;' matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments.tab | awk -F '\t'  ' $2!="limit" && ($11 >= 4) {print $1"\t"$3"\t"($6+$5)/2"\t"$9} ' > matrix_scan_pval_1e-3_GAF_Jaspar_Insects_bg_mkv_2_random_fragments_PARSED.tab
# cat /home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/bin/enrichment_by_scan/Plot_matches_extended_promoters.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " matrix.scan.active = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; matrix.scan.inactive = '/home/jcastro/Documents/JaimeCastro/PhD/Human_promoters_project/test_metrics_with_yeast_data/CapStarrseq_Active_Prom_K562_merge_IP_extended_matrix_scan_pval_1e-3_HOCOMOCO_bg_mkv_2.tab'; p.val = '1e-4'; bin = '50'; pdf.file = './test.pdf'"