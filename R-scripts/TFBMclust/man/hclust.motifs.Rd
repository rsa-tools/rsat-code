\name{hclust.motifs}
\alias{hclust.motifs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hierarchical clustering of the motifs.}
\description{This function takes as input a \code{dist} object with the distances between all the input motifs. The hierarchical clustering analysis is done with the function \code{Rclusterpp.hclust}.}
\usage{hclust.motifs(dist.matrix, hclust.method = "average")}
\arguments{
  \item{dist.table}{ This table must be of class \code{dist} and contain the distances between each pair of the input motifs. }
  \item{hclust.method}{ Hierachical clustering agglomeration rule. Default: average. Supported: average, complete, single.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{hclust.tree}{
  An object of class \code{hclust}, with the same fields as the \code{hclust} fuction. }
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hclust}
\keyword{Rcpp}
