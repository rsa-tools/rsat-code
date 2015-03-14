\name{oct4.tfbm.compa}
\alias{oct4.tfbm.compa}
\docType{data}
\title{Comparison matrix between 21 TFBM discovered by RSAT peak-motifs in Oct4 ChIP-seq peaks}
\description{This comparison matrix was produced by RSAT *compare-matrices*. It contains various comparison scores between each pair of motifs from a collection of position-specific scoring matrices discovered by *peak-motifs* in a set of 9,198 peaks for the mouse transcription factor Oct4.}
\usage{data(oct4chipseq)}
\format{
  A data frame with 231 rows, describing the similarity between 21 PSSM.
  Each row provides various comparison scores between to matrice.
  \describe{
    \item{\code{id1}}{Identifier of the first motif}
    \item{\code{id2}}{Identifier of the second motif.}
    \item{\code{name1}}{Name of the first motif.}
    \item{\code{name2}}{Name of the second motif.}
    \item{\code{label1}}{Label of the first motif (used for sequence logos).}
    \item{\code{label2}}{Label of the secondmotif (used for sequence logos).}
    \item{\code{cor}}{Correlation (Pearson's correlation coefficient).}
    \item{\code{Ncor}}{Normalized correlation.}
    \item{\code{logoDP}}{dot-porduct computed from sequence logo matrix.}
    \item{\code{NIcor}}{Normalized correlation computed on information content matrix.}
    \item{\code{NSW}}{Normalized Sandelin-Wasserman metrics.}
    \item{\code{w1}}{Width of motif 1.}
    \item{\code{w2}}{Width of motif 2.}
    \item{\code{w}}{Number of aligned columns between motifs 1 and 2.}
    \item{\code{W}}{Total width of the aligned matrices (including non-aligned columns) (W = w1+w2-w).}
    \item{\code{Wr}}{Relative width (mutual coverage) of the alignment. Wr = w / W}
    \item{\code{wr1}}{Relative width (coverage) of the alignment, respective to the first motif. wr1 = w1/W.}
    \item{\code{wr2}}{Relative width (coverage) of the alignment, respective to the second motif. wr2 = w2/W.}
    \item{\code{strand}}{Orientation of the second motif in the alignment \code{D} for direct, \code{R} for reverse complementary.}
    \item{\code{offset}}{Offset of the alignment.}
    \item{\code{consensus1}}{Aligned consensus for the first motif (including gaps for the offset).}
    \item{\code{consensus2}}{Aligned consensus for the second motif (including gaps for the offset).}
    \item{\code{rcor}}{correlation rank}
    \item{\code{rNcor}}{Ncor rank}
    \item{\code{rlogoDP}}{logoDP rank}
    \item{\code{rNIcor}}{NIcor rank}
    \item{\code{rNSW}}{NSW rank}
    \item{\code{rank.mean}}{mean rank (mean between the preceding rank scores)}
    \item{\code{match.rank}}{???}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Thomas-Chollier, M., Herrmann, C., Defrance, M., Sand, O., Thieffry, D. and van Helden, J. (2012) RSAT peak-motifs: motif analysis in full-size ChIP-seq datasets. Nucleic Acids Res, 40, e31–e31.
}
\examples{
library(TFBMclust)
data(oct4chipseq)

## The dataset contains 2 data frames

## (1) Matrix-to-matrix comparison table, with
## various scores to compare each pair of matrices
View(oct4.tfbm.compa)

## (2) Description table for the input PSSM
## (ID,name,consensus in both orientations, ...)
View(oct4.tfbm.desc)

}
\keyword{datasets}
