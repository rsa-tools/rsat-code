\name{TFBMclust-package}
\alias{TFBMclust-package}
\alias{TFBMclust}
\docType{package}
\title{Transcription factor binding motif (TFBM) clustering}
\description{This package takes as input a comparison table between transcription-factor binding motifs (TFBMs), runs hierarchical clustering to define the relationships between motifs, and cuts the tree to idenfity clusters of similar motifs.}
\details{
\tabular{ll}{
Package: \tab TFBMclust\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-21\cr
License: \tab This package is part of the RSAT software suite and inherits its license.\cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Jaime Castro-Mondragon (jcastro@lcg.unam.mx)\cr
Jacques van Helden (Jacques.van-Helden@univ-amu.fr)\cr
Maintainer: Jaime Castro-Mondragon (jcastro@lcg.unam.mx)
}
\references{
%%~~ Literature or other references for background information ~~
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
\keyword{ package
bioinformatics
cis-regulation
}
\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(TFBMclust) ## Load the library

## load a demo dataset with 21 motifs discoverd in chip-seq for Oct4
data(oct4chipseq)
}
