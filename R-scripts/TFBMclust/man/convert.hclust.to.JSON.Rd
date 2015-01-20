\name{convert.hclust.to.JSON}
\alias{convert.hclust.to.JSON}
\title{Convert \code{hclust} tree to Newick format}
\description{
This function takes one \code{hclust} object and convert it to an object of \code{character} class saving the tree in JSON format, it can be further exported in a text file.
}
\usage{
convert.hclust.to.JSON(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ One object of class \code{hclust} with the analyzis of the input motifs. The labels of the branches in the \code{hclust} tree will be the same in the JSON tree. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns an object of class \code{character} containing the Newick tree. }
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

## Save the string with the JSON tree
JSON.tree <- convert.hclust.to.JSON(tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ JSON }

