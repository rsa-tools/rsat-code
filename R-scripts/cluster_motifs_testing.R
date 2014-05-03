################################################################
## Load paameters for the demo 1: motifs discovered with peak-motif demo.
demo.nb <- 1

## How to chance this path to make it general? 
rsat.dir <- Sys.getenv("RSAT")
source(file.path(rsat.dir, "R-scripts", "cluster_motifs_lib.R"))
dir.demo <- file.path(rsat.dir, "public_html", "demo_files")

demo.prefix <- "matrix-clustering_demo_peak-motifs"
file.prefix.peakmo <- file.path(dir.demo, demo.prefix)

if (demo.nb == 1) {
  file.prefix <- file.prefix.peakmo
}


## Define parameters
score <- "Ncor";
hclust.method <- "average"

## RDB matrices
#dir.results <- "/home/jaimecm/Documents/TAGC/Clustering_test/Test_different_hclust_methods/results/RDB_TF_Fam"
#file.prefix <- file.path(dir.results, "E_coli_LacI_LysR_")

infile <- paste(sep="", file.prefix, "_pairwise_compa.tab")
description.file <-  paste(sep="", file.prefix, "_matrix_descriptions.tab")

## Create a temporary result dir if required
dir.results <- file.path(Sys.getenv("HOME"), "rsat_demo", demo.prefix)
dir.create(dir.results, showWarnings = FALSE, recursive = TRUE)
setwd(dir.results)
out.prefix <- file.path(dir.results, "clustered_motifs")
