################################################################
## Post-processing of the diff-enrichemnt result file which contains
## the distribution of TFBM occurences in a query sequence set, and 
## the significance with respect to a reference set).
##
## Author: Jacques van Helden & Jaime Castro-Mondragon

## Default parameters (can be overwritten with args)
pval.threshold <- 1e-3 ## Upper threshold on the TFBS p-value, which will be converted to a threshold on the weight in a matrix-specific manner.
occ.sig.threshold <- 5 ## Lower threshold on the significance of the number of sites. Only report matrices reaching this significance threshold. 
plot.width <- 6 ## plot width in inches
plot.heigh <- 8 ## plot height in inches
dpi <- 72 ## Dots per inch for bitmap formats

## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")


## Load some custom libraries
source(file.path(dir.rsat, 'R-scripts/config.R'))

###########################################
## Read arguments from the command line.
##
## Arguments passed on the command line
## will over-write the default arguments
## specified above.
args <- commandArgs(trailingOnly=TRUE);
## TEMP
args <- c(
  "diff.enrich.file='/no_backup/polycomb_clusters/results/motif_enrichment/between_classes/PRE_cluster_6_vs_active_1-4/PRE_cluster_6_vs_active_1-4.tab'",
  "out.prefix='/no_backup/polycomb_clusters/results/motif_enrichment/between_classes/PRE_cluster_6_vs_active_1-4/PRE_cluster_6_vs_active_1-4'"
)

## Parse the arguments and create a variable for each one
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  verbose(args, 3)
}


## Check that input file has been specified
if (!exists("diff.table")) {
  stop("Missing mandatory argument: diff.table=[diff_enrichment_result.tab] ")
}
verbose(paste("Differential enrichment file", diff.table), 3)

## Check that output prefix has been specified
if (!exists("out.prefix")) {
  stop("Missing mandatory argument: out.prefix=[output_prefix] ")
}
verbose(paste("Output prefix", out.prefix), 3)


## Read the diff-enrichment file
diff.table <- read.delim(file = diff.enrich.file, comment.char = ";", header = TRUE)
## Suppress the X. in first column header, which comes from the parsing of RSAT # header character. 
names(diff.table)[1] <- sub(pattern = "^X.", replacement = "", x =  names(diff.table)[1])
# View(diff.table)
# names(diff.table)
# dim(diff.table)


matrices <- as.vector(unique(diff.table$matrix))

## Iterate over matrices
current.matrix <- matrices[1]
plot.dir <- file.path(out.prefix, "plots")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)
message("Plot dir: ", plot.dir)

for (current.matrix in matrices) {
  diff.one.matrix <- subset(x = diff.table, matrix == current.matrix)
  # dim(diff.one.matrix)
  plot.format <- ("pdf")  
  
  for (plot.format in c("pdf", "png")) {
    plot.file <- file.path(plot.dir, paste(sep="", current.matrix, "_enrich_plots.", plot.format))  
    
    
    if (plot.format == "pdf") {
      pdf(file = plot.file, width = plot.width, height = plot.heigh)
    } else if (plot.format == "png") {
      png(file = plot.file, width = plot.width*dpi, height = plot.heigh*dpi)
    }
    
    par(mfrow=c(2,2))
    par.ori <- par(no.readonly = TRUE) ## Take a copy of original plotting parameters
    
    
    ## Plot observed and expected distributions of TFBS occurrences
    par(mar=c(1.1, 4.1, 3.1, 1.1))
    ## Jaime: to fix: the ylab should be on the left, in order to leave space for big numbers (1 million)
    plot(diff.one.matrix$score, 
         diff.one.matrix$inv_cum,
         main=current.matrix, 
         ylab="Number of sites (log scale)",
         cex.axis=0.8,
         las=1, type="l", 
         col="#008800", 
         panel.first=grid(), 
         log="y")
    
    lines(diff.one.matrix$score, diff.one.matrix$exp_occ, 
          col="blue")
    abline(v=0)
    
    ## Significance plot
    par(mar=c(4.1, 4.1, 1.1, 1.1))
    plot(diff.one.matrix$score, diff.one.matrix$occ_sig, 
         xlab = "Weight score", 
         ylab="Significance",
         panel.first=grid(), las=1, 
         type="l", ylim=c(0,50), col="#BB00BB")
    abline(h=occ.sig.threshold, col="red")
    abline(v=0)
    
    par(mfrow=c(1,1))
    par <- par.ori ## Restore original plotting parameters
    
    ## Significance as a function of the site p-value
    #plot(diff.one.matrix$)
    
    ## PARA JAIME
    
    dev.off()
  }
}  

