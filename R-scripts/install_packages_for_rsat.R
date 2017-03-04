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
required.packages = c("devtools",
                      "RJSONIO",
                      "gplots",
                      "jpeg",
                      "png",
                      "dynamicTreeCut",
                      "ggplot2",
                      "zoo",
                      "reshape2",
                      "dendextend",
                      "RColorBrewer",
                      "gridExtra",
                      "cowplot",
                      "flux")

#                       "Rcpp",
#                       "RcppEigen",
#                       "Rclusterpp",

## List of required packages from Bioconductor
required.packages.bioconductor <- c("ctc", "amap", "qvalue")

## List of RSAT-specific packages to be compiled on the server
required.packages.rsat <- c("TFBMclust")

## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop("The environment variable RSAT is not defined.")
}
dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")
dir.create(dir.rsat.rlib, showWarnings = FALSE, recursive = FALSE)

## Install R packages from the CRAN
message("Installing R packages from CRAN repository: ", rcran.repos)
#message(required.packages)
for (pkg in required.packages) {
  if(suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE, lib=c(.libPaths(),dir.rsat.rlib)))) {
    message(pkg, " CRAN package already installed. Skipping. ")
  } else {
    message("Installing CRAN package ", pkg, " in dir ", dir.rsat.rlib)
    install.packages(pkg, repos=rcran.repos, dependencies=TRUE, lib=dir.rsat.rlib)
    message(pkg, " CRAN package installed in dir ", dir.rsat.rlib)
  }
}

################################################################
## Check requirement for bioconductor packages
message("Installing BioConductor packages")
message(cat("Required BioConductor packages: ", required.packages.bioconductor))
for (pkg in required.packages.bioconductor) {
  if (suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE, lib=c(.libPaths(),dir.rsat.rlib)))) {
    message(pkg, " BioConductor package already installed. Skipping. ")
  } else {
    message("Installing Bioconductor package ", pkg, " in dir ", dir.rsat.rlib)
    .libPaths(c(dir.rsat.rlib, .libPaths())) ## this line fixes the problem at ENS (Morgane)
    source("http://bioconductor.org/biocLite.R")
#    biocLite(ask=FALSE, lib=dir.rsat.rlib,  lib.loc=dir.rsat.rlib)
    biocLite(lib=dir.rsat.rlib, lib.loc=c(.libPaths(),dir.rsat.rlib))
    biocLite(pkg, dependencies=TRUE, lib=dir.rsat.rlib,  lib.loc=dir.rsat.rlib)
    message(pkg, " BioConductor package installed in dir ", dir.rsat.rlib)
  }
}

# ## Load required libraries
# for (pkg in c(required.packages, required.packages.bioconductor)) {
#   suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE))
# }

################################################################
## Install RSAT-specific packages.  We force re-installation of these
## packages since they may have changed since the last installation.
message("Installing RSAT-specific packages")
message(required.packages.rsat)
for (package in required.packages.rsat) {
  message("Installing RSAT package ", package, " in folder ", dir.rsat.rlib)
  install.packages(pkgs=file.path(dir.rsat.rscripts, "TFBMclust"), repos=NULL,  lib=dir.rsat.rlib, type="source")
}




