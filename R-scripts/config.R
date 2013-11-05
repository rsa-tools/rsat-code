################################################################
#
# This files defines the default configuration for the course
# Statistics Applied to Bioinformatics It should be loaded before
# starting the practicals. All the rest will be done via internet.  #
#
# If you install the course on your computer, you only need to change
# the variable dir.main

#### Load files from the web server
dir.rsat <- Sys.getenv("RSAT");
if (dir.rsat == "") {
  stop("The environment variable RSAT is not defined.")
}

## Directory with the R scripts for the course
dir.R.files <- file.path(dir.rsat, "R-scripts")
print(paste("R scripts", dir.R.files))
dir.util <- file.path(dir.R.files, 'util')

################################################################
## Load utilities
source(file.path(dir.util, 'util.R'))

################################################################
# Global parameters
verbosity <- 1
export.formats.plots <- c("eps","pdf", "png")
export.formats.obj <- c("table")

## ################################################################
## ## Graphic device Special fix for a bug in Mac OSX 10.5: X11() has
## ## different fonts from the export devices (pdf, png, eps), which
## ## creates problem with the legends (legends are cut bu the box
## ## limit). The problem is fixed by opening a Mac-specific graphical
## ## window using quartz() instead of X11().
## sysinfo <- Sys.info()
## if (sysinfo["sysname"] == "Darwin") {
##   X11 <- function (...) {
##     quartz(...)
##   }
## }
