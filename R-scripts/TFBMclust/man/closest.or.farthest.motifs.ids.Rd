\name{closest.or.farthest.motifs.ids}
\alias{closest.or.farthest.motifs.ids}
\title{ Identify the "central" or farthest leaves (motifs) between two clusters. }
\description{
Identify the "central" or farthest leaves (motifs) between two subtrees (clusters). i.e. the closest/farthest pair of motif among all the other motifs of these clusters. This functions take in consideration if the score (metric) selected is a distance or a correlation.
}
\usage{
central.motifs.ids(id1, id2, compa.table, metric = "Ncor", closest = TRUE)
}
\arguments{
  \item{id1}{ Cluster 1 IDs. }
  \item{id2}{ Cluster 2 IDs. }
  \item{compa.table}{ A comparison matrix containing at least 3 columns: \code{ID1}, \code{ID2} and one or more comparison metrics between each pair of motifs. }
  \item{metric}{ Column of the matrix to be used to build the hierarchical tree. Default: Ncor. }
  \item{closest}{ logical.  Should be searched the closest or farthest motifs?. Default: TRUE}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ The IDs of the pair of closest motifs between the two  input clusters. }
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
## Load the Oct4 chip-seq demo dataset
data(oct4chipseq)

## Collect the motif distance matrix
oct4.dist <- build.distance.matrix(oct4.tfbm.compa, metric = "Ncor")

## Save the matrix
oct4.dist.matrix <- oct4.dist$matrix

## Run the hierarchical clustering analysis
tree <- hclust.motifs(oct4.dist.matrix, hclust.method = "average")

## Save the list with the leaves
leaves.per.level <- leaves.per.node(tree)

## Get the IDs of the levels 10 and 12
IDs1 <- get.id(leaves.per.level[[10]], oct4.tfbm.desc)
IDs2 <- get.id(leaves.per.level[[12]], oct4.tfbm.desc)

## Get the closest pair of motifs
closest.or.farthest.motifs.ids(IDs1, IDs2, oct4.tfbm.compa, metric = "Ncor", closest = TRUE)

## Get the farthest pair of motifs
closest.or.farthest.motifs.ids(IDs1, IDs2, oct4.tfbm.compa, metric = "Ncor", closest = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Closest nodes }
\keyword{ Farthest nodes }
