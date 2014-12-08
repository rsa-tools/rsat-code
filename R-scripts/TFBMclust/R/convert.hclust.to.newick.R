#####################################################
## Convert the hclust object to a character object
## with the lines ready to print a Newick file
convert.hclust.to.newick <- function(tree, decimals = 3){

  ## Require ctc if it is required
  if(!require("ctc")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("ctc")
  }

  ## Load library
  suppressPackageStartupMessages(library("ctc", warn.conflicts=FALSE))

  ## Round the decimals of the branches distances
  tree[[2]] <- round(tree[[2]], digits = decimals)

  ## Convert the hclust tree to Newick format
  ## If ‘flat=TRUE’ the result is a string can be written in a file                                     file).
  return(hc2Newick(tree, flat = TRUE))
}