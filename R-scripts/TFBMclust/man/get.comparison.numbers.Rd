\name{get.comparison.numbers}
\alias{get.comparison.numbers}
\title{ Given two sets of IDs, return the number of comparison from the given comparison table }
\description{
Return the comparison numbers for all the pairs of the leaves belonging to the two IDs (clusters) sets.
}
\usage{
get.comparison.numbers(id1, id2, comparison.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id1}{ Cluster 1 IDs }
  \item{id1}{ Cluster 2 IDs }
  \item{comparison.table}{ A comparison matrix containing at least 4 columns: \code{ID1}, \code{ID2}, \code{strand} and one or more comparison scores between each pair of motifs. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the comparison number between each pair of motifs among the two clusters.
}
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
oct4.dist <- build.distance.matrix(oct4.tfbm.compa, score = "Ncor")

## Save the matrix
oct4.dist.matrix <- oct4.dist$matrix

## Run the hierarchical clustering analysis
tree <- hclust.motifs(oct4.dist.matrix, hclust.method = "average")

## Save the list with the leaves
leaves.per.level <- leaves.per.node(tree)

## Get the IDs of the levels 10 and 12
IDs1 <- get.id(leaves.per.level[[10]], oct4.tfbm.desc)
IDs2 <- get.id(leaves.per.level[[12]], oct4.tfbm.desc)

## Get the comparison numbers
get.comparison.numbers(IDs1, IDs2, oct4.tfbm.compa)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
