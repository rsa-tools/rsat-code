################################################################
## Install R packages required for different RSAT applications.
##
## This script requires the write access to the folder RSAT/R-scripts,
## and should thus be executed only by the administator of the RSAT suite.

## Define the preferred CRAN repository.
## If not defined in an environment variable, take rstudio by default.
if (Sys.getenv("CRAN_REPOS") == "") {
  rcran.repos <- "http://cran.rstudio.com/"
} else {
  rcran.repos <- Sys.getenv("CRAN_REPOS")
}

## List of packages to install
required.packages = c("RJSONIO",
                      "dendextend",
                      "Rcpp",
                      "Rclusterpp",
                      "gplots",
                      "devtools")

## List of required packages from Bioconductor
required.packages.bioconductor <- c("ctc")

## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop("The environment variable RSAT is not defined.")
}
dir.rsat.rlib <- file.path(dir.rsat, "R-scripts", "Rpackages")

## Install R packages from the CRAN
for (pkg in required.packages) {
  if(!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    install.packages(pkg, repos=rcran.repos, dependencies=TRUE, lib=dir.rsat.rlib)
  }
}


## Check requirement for bioconductor packages
for (pkg in required.packages.bioconductor) {
  if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(lib=dir.rsat.rlib);
    biocLite(pkg, lib=dir.rsat.rlib)
    #  library("ctc", quietly=TRUE, character.only = TRUE)
  }
}

## Load required libraries
for (pkg in c(required.packages, required.packages.bioconductor)) {
  suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE))
}

## Update the package TFBMclust
print("Installing RSAT TFBM package")
install.packages(pkgs=file.path(dir.rsat.rscripts, "TFBMclust"), repos=NULL,  lib=dir.rsat.rlib, type="source", INSTALL_opts=c("no-multiarch"))
#system(paste("R CMD INSTALL --no-multiarch --with-keep.source  \"", file.path(dir.rsat, "R-scripts/TFBMclust"), "\"", sep =""))
# reload(file.path(dir.rsat, "R-scripts/TFBMclust"))
suppressPackageStartupMessages(library(TFBMclust, warn.conflicts=FALSE))




