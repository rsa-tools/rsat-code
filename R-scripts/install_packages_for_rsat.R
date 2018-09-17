################################################################
## Install R packages required for different RSAT applications.
##
## This script requires the write access to the folder RSAT/R-scripts,
## and should thus be executed only by the administator of the RSAT suite.

## Define the preferred CRAN repository.
## If not defined in an environment variable, take rstudio by default.
if (Sys.getenv("CRAN_REPOS") == "") {
  rcran.repos <- "https://cran.rstudio.com/"
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
                      "reshape2",
                      "dendextend",
                      "gridExtra",
                      "grid",
                      "egg",
                      "flux",
                      "zoo",
                      "RColorBrewer",
                      "changepoint",
                      "tidyverse" # includes ggplot2, dplyr
#                       "Rcpp",
#                       "RcppEigen",
#                       "Rclusterpp",
)

## List of required packages from Bioconductor
required.packages.bioconductor <- c("S4Vectors",
				    "ctc",
                                    "qvalue",
				    "IRanges", 
                                    "BiocGenerics",
                                    "amap", 
                                    "GenomicRanges"
)

## List of RSAT-specific packages to be compiled on the server
required.packages.rsat <- c("TFBMclust")

## Define the local directory for R librairies
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop("The environment variable RSAT is not defined.")
}
dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")

################################################################
## By default, install all packages in the first element of
## .libPaths(). However if you are not admin user, they can be
## installed in the install rsat libraries.
install.dir <- dir.rsat.rlib
#install.dir <- .libPaths()[1]

dir.create(install.dir, showWarnings = FALSE, recursive = FALSE)

## Install R packages from the CRAN (type="source" is default)
message("Installing R packages from CRAN repository: ", rcran.repos)
message("Installing R packages in local directory: ", install.dir)

##message(required.packages)
# pkg <- required.packages[1]
for (pkg in required.packages) {
    if (suppressPackageStartupMessages(
        require(pkg, quietly=TRUE,
                character.only = TRUE,
                lib=c(.libPaths(),install.dir)))) {
        message(pkg, "         CRAN package already installed. Skipping. ")
    } else {
        message("Installing CRAN package ", pkg, " in dir ", install.dir)
        install.packages(pkg, repos=rcran.repos, dependencies=TRUE, lib=install.dir)
        message(pkg, " CRAN package installed in dir ", install.dir)
    }
}



################################################################
## Check requirement for bioconductor packages
message("Installing BioConductor packages")
message(cat("Required BioConductor packages: ", required.packages.bioconductor))
for (pkg in required.packages.bioconductor) {
  if (suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE, lib=c(.libPaths(),install.dir)))) {
    message(pkg, " BioConductor package already installed. Skipping. ")
  } else {
    message("Installing Bioconductor package ", pkg, " in dir ", install.dir)
#    .libPaths(c(install.dir, .libPaths())) ## this line fixes the problem at ENS (Morgane) ## MAYBE BUT IT CREATES PROBLEMS IF THE USER HAS NOT DONE IT IN THE USE S
    source("http://bioconductor.org/biocLite.R")
#    biocLite(ask=FALSE, lib=install.dir,  lib.loc=install.dir)
    biocLite(lib=install.dir, lib.loc=c(.libPaths(),install.dir), suppressUpdates="")
    biocLite(pkg, dependencies=TRUE, lib=install.dir,  lib.loc=install.dir, suppressUpdates="")
    message(pkg, " BioConductor package installed in dir ", install.dir)
  }
}

# ## Load required libraries
# for (pkg in c(required.packages, required.packages.bioconductor)) {
#   suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE))
# }

################################################################
## Install RSAT-specific packages.  We force re-installation of these
## packages since they may have changed since the last installation.
message("Installing RSAT-specific packages in ", dir.rsat.rlib)

dir.create(dir.rsat.rlib, recursive=TRUE, showWarnings=FALSE)
message(required.packages.rsat)

for (package in required.packages.rsat) {
  message("Installing RSAT package ", package, " in folder ", dir.rsat.rlib)
  install.packages(pkgs=file.path(dir.rsat.rscripts, package), repos=NULL,  lib=dir.rsat.rlib, type="source")
}




