################################################################
## Post-processing of the diff-enrichemnt result file which contains
## the distribution of TFBM occurences in a query sequence set, and 
## the significance with respect to a reference set).
##
## Author: Jacques van Helden & Jaime Castro-Mondragon



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
args <- "diff.enrich.file='/no_backup/polycomb_clusters/results/motif_enrichment/between_classes/PRE_cluster_6_vs_active_1-4/PRE_cluster_6_vs_active_1-4.tab'"

## Parse the arguments and create a variable for each one
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  verbose(args, 3)
}

## Read the diff-enrichment file
diff.table <- read.delim(file = diff.enrich.file, comment.char = ";", header = TRUE)
## Suppress the X. in first column header, which comes from the parsing of RSAT # header character. 
names(diff.table)[1] <- sub(pattern = "^X.", replacement = "", x =  names(diff.table)[1])
# View(diff.table)
# names(diff.table)


