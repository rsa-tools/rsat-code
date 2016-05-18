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

# prefix
# TFs 
# Sequences
# html.template.file

if (!exists("prefix")) {
  stop("Missing mandatory argument (prefix): prefix ")
} else if (!exists("html.template.file")){
  stop("Missing mandatory argument (HTML template to draw the profiles): html.template.file ")
} else if (!exists("Sequences")){
  stop("Missing mandatory argument (Sequences names): Sequences ")
} else if (!exists("TFs")){
  stop("Missing mandatory argument (TF names): TFs ")
}

## Required libraries
library("RColorBrewer")

##########################
## Initialize variables
all.IDs <- NULL
# all.names <- all.profiles ## Remember!
hash.profile.ID <- list()
plot.names <- NULL
ID.counter <- 0
x.y <- NULL 
line.w <- 5

# prefix <- "/home/jaimicore/Documents/PhD/motif_enrichment/20160218/capstarr-seq"
# TFs <- c("DDIT", "FCP2_GRHL1", "JUN_FOS", "REST", "SP", "USF2", "YY")
# Sequences <- c("HELA", "K562", "inactive")

# html.template.file <- "/home/jaimicore/Documents/PhD/motif_enrichment/dynamic_OCC_profiles.html"
base.name <- basename(prefix)
dir.name <- dirname(prefix)

## Get sequences names
Sequences <- unlist(strsplit(Sequences, "---", perl = TRUE))

## Get TF names
TFs <- unlist(strsplit(TFs, "---", perl = TRUE))

# cat /home/jaimicore/Documents/PhD/motif_enrichment/20160218/plot_binomial_occ.R | /usr/bin/R --slave --no-save --no-restore --no-environ --args " prefix = '/home/jaimicore/Documents/PhD/motif_enrichment/20160218/capstarr-seq'; html.template.file = '/home/jaimicore/Documents/PhD/motif_enrichment/dynamic_OCC_profiles.html'; Sequences = 'HELA---K562---inactive---'; TFs = 'DDIT---FCP2_GRHL1---JUN_FOS---REST---SP---USF2---YY---';"



prof <- as.vector(outer(Sequences, TFs, paste))
all.profiles <<- as.vector(sapply(prof, function(p){gsub(" ", "_", p)}))

pdf.file <- paste(prefix, "_Binomial_occ_plots.pdf", sep = "")
pdf(pdf.file )
for(TF in TFs){
  
  ## Read Binomial OCC table
  binomial.occ.file <- paste(dir.name, "/", TF, "/", base.name, "_", TF, "_occ_proba_", TF,"_compare-scores.tab", sep = "")
  binomial.occ.tab <- read.csv(binomial.occ.file, header = TRUE, sep = "\t")
  
  col.nb <- dim(binomial.occ.tab)[2]
  colors <- colorRampPalette(brewer.pal(5, "Set1"), space="Lab")(col.nb-1)
  
  old.names <- names(binomial.occ.tab)
  names(binomial.occ.tab) <- c("values", as.character(2:col.nb))
  
  ## Convert the <NULL> values (coming from previous steps) to NA
  ## in order to have a numeric vector
  binomial.occ.tab.cp <- apply(binomial.occ.tab, 2, function(x){
    
    x <- as.vector(x)
    x <- as.numeric(gsub("<NULL>", NA, x))
    return(x)
  })
  
  calculate.max.y <- max(binomial.occ.tab.cp[,2:col.nb], na.rm = TRUE)
  calculate.min.y <- min(binomial.occ.tab.cp[,2:col.nb], na.rm = TRUE)
  
  x <- binomial.occ.tab.cp[,1]
  y <- binomial.occ.tab.cp[,2]
  
  ## X-axis values scaled at -log10
  x <- -log10(x)
  
  ## Plot the Binomial OCC separately of each TF
  plot(x,y, 
       type = "l", 
       col = colors[1],
       main = paste("Binomial OCC", TF, sep = " "),
       ylim = c(calculate.min.y - 1, calculate.max.y + 1),
       lwd = 2,
       xlab = "P-value",
       ylab = "Binomial OCC significance"
  )
  
  ## Add the lines only if there are 2 or more sequences
  if(col.nb > 2){
    for(c in 3:col.nb){
      lines(x = x, y = binomial.occ.tab.cp[,c],
            type = "l", 
            col = colors[c-1],
            lwd = 2)
    }
  }
  legend("topright", Sequences, fill = colors, cex = 0.75)
  
  ################################
  ## Prepare C3 html fields     
  ## (Before create the report) 
  
  ## Define all profiles of the current TF
  set.profiles <- paste(TF, "_", Sequences, sep = "")
  
  ## Set colors (one per profile)
  set.colors <- colorRampPalette(brewer.pal(10,"Paired"))(length(all.profiles))
  
  binomial.occ.tab <- binomial.occ.tab.cp
  binomial.occ.tab[,1] <- -log10(binomial.occ.tab.cp[,1])
  colnames(binomial.occ.tab) <- c("values", set.profiles)
  binomial.occ.tab <- t(binomial.occ.tab)
  
  ## Rename the columns
  colnames(binomial.occ.tab) <- binomial.occ.tab[1,]
  binomial.occ.tab <- binomial.occ.tab[2:(length(set.profiles)+1),]
  sub.tab.nb.col <- dim(binomial.occ.tab)[1]
  
  ###################################################
  ## Get the X and Y data
  ## As there are some NA values, they are removed 
  ## from the data to draw them properly in C3
  thrash <- apply(binomial.occ.tab, 1, function(values){
    
    ## Count the IDs
    ID.counter <<- ID.counter + 1
    
    ## Here we create a unique ID without CSS special characters
    ## Only to manipulate the objects in the HTML form
    ID <- all.profiles[ID.counter]
    ID <- gsub("_", "", ID)
    ID <- gsub("_", "", ID)
    ID <- gsub("-", "", ID)
    ID <- gsub("\\.", "", ID)
    ID <- gsub(":", "", ID)
    ID <- gsub("\\s+", "", ID, perl = TRUE)
    ID <- paste(ID.counter, ID, ID.counter, sep = "") 
    
    hash.profile.ID[[all.profiles[ID.counter]]] <<- ID
    all.IDs <<- append(all.IDs, ID)
    
    ## Create the name's data for the HTML file
    ## The name correspond to the motif name. 
    ## Note that two motifs for the same TF will have the same name and this name will appears twice 
    ## in the report. However their IDs are unique.
    plot.names <<- append(plot.names, paste("'", ID, "' : '",  all.profiles[ID.counter],"',", sep = ""))
    
    ## PArse the values vector
    ## Only extract those positions without NA
    values.wo.NA <- values[which(values != "NA")]
    
    a <- sum(values.wo.NA)
    b <- cumsum(values.wo.NA)
    print("- - - - - - - - - - - - - - - -- - - -  - - - -")
    print(paste("Profile: ", all.profiles[ID.counter]," - Sum: ",a, " - CumSum: ", b[length(b)], sep = ""))
    
    # Create the Y value
    y <- paste("['",
               ID,
               "',",
               paste(values.wo.NA,
                     collapse = ","),
               "],",
               sep = "")
    x.y <<- rbind(x.y, y)
    
    ## Add the X-value
    
    aaa <- as.numeric(names(values.wo.NA))
    aaa <- prettyNum(aaa, scientific=TRUE, digits = 6)
    
    
    x.y <<- rbind(x.y, paste("['x", ID.counter,"',", paste(aaa,collapse = ","), "],", sep = ""))
  })
}
dev.off()

#############################
## Create the HTML c3 plot ##
#############################
html.report <- readLines(html.template.file)

## Print Motif names array
html.names <- sapply(TFs, function(x){
  rep(x, times = length(Sequences))
})
html.names <- paste("Names['", all.IDs, "'] = '", html.names, "';", sep = "")
html.names <- paste(html.names, collapse = "\n")
html.report <- gsub("--names_vector--", html.names, html.report)

## Print Sequence names array
html.seq <- rep(Sequences, times = length(Sequences))
html.seq<- as.vector(html.seq)
html.seq <- paste("Seqs['", all.IDs, "'] = '", html.seq, "';", sep = "")
html.seq <- paste(html.seq, collapse = "\n")
html.report <- gsub("--seqs--", html.seq, html.report)

## Write the logo's path
## for both orientations
## Each path is repeated N-times where N is the 
## number of sequences
logos.F <- sapply(TFs, function(i){
  paste(i, "/", base.name , "_", i, "_logo_m1.png", sep = "")
})
logos.F <- as.vector(logos.F)
logos.F <- sapply(logos.F, function(x){
  rep(x, times = length(Sequences))
})
logos.F <- paste("pics['", all.IDs, "'] = '", logos.F, "';", sep = "")
logos.F <- paste(logos.F, collapse = "\n")
html.report <- gsub("--pics--", logos.F, html.report)


logos.R <- sapply(TFs, function(i){
  paste(i, "/", base.name , "_", i, "_logo_m1_rc.png", sep = "")
})
logos.R <- as.vector(logos.R)
logos.R <- sapply(logos.R, function(x){
  rep(x, times = length(Sequences))
})
logos.R <- paste("pics_rc['", all.IDs, "'] = '", logos.R, "';", sep = "")
logos.R <- paste(logos.R, collapse = "\n")
html.report <- gsub("--pics_rc--", logos.R, html.report)


## Add the color code (one color per motif)
## They are inserted in the C3 section
colors <- colorRampPalette(brewer.pal(5, "Set1"), space="Lab")(length(all.IDs))
set.colors <- paste(paste("'", colors, "'", sep = ""), collapse = ",")
html.report <- gsub("--color_pattern--", set.colors, html.report)

## CSS section to set the line width
## Note: the width is proportional cummulative sum of the profile (TO DO)
## For the moment, only a fixed value is working
line.w <- paste("#chart .c3-line-", all.IDs, "{ stroke-width: ", line.w, "px; }", sep = "")
line.w <- paste(line.w, collapse = "\n")
html.report <- gsub("--lines_w--", line.w, html.report)

## Add the TF_names data
TF.names <- paste("TF_names['", all.profiles, "'] = '", all.IDs, "';", sep = "")
TF.names <- paste(TF.names, collapse = "\n")
html.report <- gsub("--TF_names--", TF.names, html.report)

## Add the TF_names data
tfs <- paste(paste("'", all.profiles, "'", sep = ""), collapse = ",")
html.report <- gsub("--tfs--", tfs, html.report)

## Add the real motif IDs (to display in the tooltip)
## They are inserted in the JS section
## I called them 'real' because are those found on the original motif file
IDs <- paste("IDs['", all.IDs, "'] = '", all.profiles, "';", sep = "")
IDs <- paste(IDs, collapse = "\n")
html.report <- gsub("--IDs--", IDs, html.report)

## The plot heigth depends in the number of motifs
motif.total <- length(all.profiles)
chart.heigth <- 500
if(motif.total >= 300){
  chart.heigth <- 700
} else if(motif.total >= 400){
  chart.heigth <- 900
}
html.report <- gsub("--chart_h--", chart.heigth, html.report)

## Add xs values (one row per motif)
## They are inserted in the C3 section
xs <- paste("'", all.IDs, "' : 'x", 1:ID.counter, "',", sep = "")
xs <- paste(xs, collapse = "\n")
html.report <- gsub("--xs--", xs, html.report)

## Add x values (one row per motif)
## They are inserted in the C3 section
xx <- paste(x.y, collapse = "\n")
xx <- gsub(",NA,", ",NaN,", xx)
html.report <- gsub("--x_y--", xx, html.report)

## Insert the motif names
## They are inserted in the C3 section
plot.names <- paste(plot.names, collapse = "\n")
html.report <- gsub("--names--", plot.names, html.report)

## Insert the motif names (to hide/show all)
## They are inserted in the JQuery section
all.IDs <- paste(paste("'", all.IDs, "'", sep = ""), collapse = ",")
html.report <- gsub("--all--", all.IDs, html.report)

# ## Insert the JavaScript libraries path
# html.report <- gsub("--c3_css--", c3.css.base, html.report)
# html.report <- gsub("--c3--", c3.base, html.report)

## Export the report
html.report.file <- paste(prefix, "_binomial_occ_sig_profiles.html", sep = "")
write(html.report, file = html.report.file)