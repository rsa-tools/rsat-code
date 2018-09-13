################################################################
## Update R CRAN packages required for different RSAT applications.
##
## This script requires write access to the folder RSAT/R-scripts,
## and should thus be executed only by the administator of the RSAT suite.

## Define the preferred CRAN repository.
## If not defined in an environment variable, take rstudio by default.
if (Sys.getenv("CRAN_REPOS") == "") {
  rcran.repos <- "http://cran.rstudio.com/"
} else {
  rcran.repos <- Sys.getenv("CRAN_REPOS")
}

## List of required packages from Bioconductor
required.packages.bioconductor <- c("ctc", "amap", "qvalue", "GenomicRanges")

## List of RSAT-specific packages to be compiled on the server
#required.packages.rsat <- c("TFBMclust")

## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop("The environment variable RSAT is not defined.")
}
dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")

## List of CRAN packages to install
required.packages = installed.packages(lib.loc=dir.rsat.rlib)

################################################################
## By default, install all packages in the first element of
## .libPaths(). However if you are not admin user, they can be
## installed in the install rsat libraries.
install.dir <- dir.rsat.rlib

dir.create(install.dir, showWarnings = FALSE, recursive = FALSE)

## Install R packages from the CRAN
message("Updating R packages from CRAN repository: ", rcran.repos)
message("Installing R packages in local directory: ", install.dir)

update.packages(instPkgs=required.packages, repos=rcran.repos, lib=install.dir, ask = FALSE, checkBuilt = TRUE)

message("Updating BioConductor packages")

for (pkg in required.packages.bioconductor) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg, dependencies=TRUE, lib=install.dir,  lib.loc=install.dir, suppressUpdates="")
}




