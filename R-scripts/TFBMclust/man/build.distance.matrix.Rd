\name{build.distance.matrix}
\alias{build.distance.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build distance matrix from a motif comparison table}
\description{This function takes as input a motif comparison table, and collects the scores from a user-specified column in order to build a motif-to-motif distance matrix, which will be further used as input for hierarchical clustering. Importantly, the function also "interprets" the score as either a distance or a similarity, and, in case of similarity scores, converts it to a distance score. }
\usage{
build.distance.matrix(comparison.table, metric="Ncor")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comparison.table}{This table must contain at least three columns indicating the ID of the first motif (first column), the ID of the second matrix (second column), and the score (user-specified column).
  }
  \item{metric}{ Column of the comparison matrix to be used as score. Default: Ncor. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with two fields:
\item{$matrix}{A distance object providing the score for each pair of motif (lower triangle only, since the score is assumed to be symmetrical).}
\item{$table}{A matrix providing the score for each pair of motif, in both directions (the upper and lower triangles are filled up).}
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
oct4.dist <- build.distance.matrix(oct4.tfbm.compa, metric="Ncor")
}
\keyword{datasets}
