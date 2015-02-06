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
                      "dendextend",
                      "Rcpp",
                      "RcppEigen",
                      "Rclusterpp",
                      "gplots")

## List of required packages from Bioconductor
required.packages.bioconductor <- c("ctc")

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
print(paste("Installing R packages from CRAN repository", rcran.repos))
#print(required.packages)
for (pkg in required.packages) {
  if(!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
    install.packages(pkg, repos=rcran.repos, dependencies=TRUE, lib=dir.rsat.rlib)
    print(paste(pkg, "CRAN package installed in dir", dir.rsat.rlib))
  }
}


################################################################
####### THIS DOES NOT WORK AT ENS, where we are not sudoers !!!
## Check requirement for bioconductor packages
print("Installing BioConductor packages")
print(required.packages.bioconductor)
for (pkg in required.packages.bioconductor) {
  if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE, lib=c(.libPaths(),dir.rsat.rlib)))
      ) {
    source("http://bioconductor.org/biocLite.R")
#    biocLite(ask=FALSE, lib=dir.rsat.rlib,  lib.loc=dir.rsat.rlib)
    biocLite(lib=dir.rsat.rlib, lib.loc=c(.libPaths(),dir.rsat.rlib))
    biocLite(pkg, dependencies=TRUE, lib=dir.rsat.rlib,  lib.loc=dir.rsat.rlib)
    print(paste(pkg, "BioConductor package installed in dir", dir.rsat.rlib))
  }
}

# ## Load required libraries
# for (pkg in c(required.packages, required.packages.bioconductor)) {
#   suppressPackageStartupMessages(library(pkg, warn.conflicts=FALSE, character.only = TRUE))
# }


## Install RSAT-specific packages
print("Installing RSAT-specific packages")
print(required.packages.rsat)
for (package in required.packages.rsat) {
  print(paste("Installing RSAT package", package, "in folder", dir.rsat.rlib))
  install.packages(pkgs=file.path(dir.rsat.rscripts, "TFBMclust"), repos=NULL,  lib=dir.rsat.rlib, type="source")
}




