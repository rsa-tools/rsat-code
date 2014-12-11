\name{inverted.alignment}
\alias{inverted.alignment}
\title{ Invert the aligned consensus and orientation information for a group of aligned motifs. }
\description{ This function takes as input a set of motifs IDs and search these motifs in the motifs.list, takes the consensus and takes from the desc.table the consensus on its opposite orientation, then add the gaps ("-") and finally update the new information. The modified IDs are exported in a \code{list} with the same attributes than the \code{list} object \code{motifs.info}. }
\usage{ inverted.alignment(ids, motifs.list, desc.table) }
\arguments{
  \item{ids}{ A \code{vector} with the ID of the motif(s). }
  \item{motifs.list}{ An object of class \code{list} containing one \code{list} for each motif analyzed. It must contain the next attributes: \code{name}, \code{consensus}, \code{strand}, \code{number}, \code{spacer.up} and \code{spacer.dw}. The \code{names()} of each list must be the motif \code{ID}. This object is inherited from the function \code{align.motifs()} which is the function where this function is called. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name},     \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID.} }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a \code{list} object with the same attributes than \code{motifs.list} with the inverted consensus and orientation. }
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Load the data
data(oct4chipseq)

## Collect the motif distance matrix
oct4.dist <- build.distance.matrix(oct4.tfbm.compa, score = "Ncor")

## Save the matrix
oct4.dist.matrix <- oct4.dist$matrix

## Run the hierarchical clustering analysis
tree <- hclust.motifs(oct4.dist.matrix, hclust.method = "average")

## Align the motifs
oct4.alignment <- align.motifs(tree, oct4.tfbm.desc, oct4.tfbm.compa, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), score = "Ncor", method = "average", metric = "Ncor", nodes.attributes = TRUE, intermediate.alignments = FALSE)

## Save the alignment of all the motifs
oct4.alignment.motifs <- oct4.alignment$motifs.alignment

## Takes two motifs and print their consensuses
motifs.1.and.2 <- oct4.alignment.motifs[1:2]

sapply(motifs.1.and.2, function(x){print(x[["consensus"]])})

## Get the ids od the motifs 1 and 2
## Then invert the orientation and print the consensuses
ids <- names(oct4.alignment.motifs[1:2])

motifs.1.and.2.inverted <- inverted.alignment(ids, motifs.1.and.2, oct4.tfbm.desc)

sapply(motifs.1.and.2.inverted, function(x){print(x[["consensus"]])})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Inverted Alignment }
