\name{convert.hclust.to.newick}
\alias{convert.hclust.to.newick}
\title{ Convert \code{hclust} tree to Newick format. }
\description{
This function takes one \code{hclust} object and convert it to an object of \code{character} class saving the tree in Newick format, it can be further exported in a text file. The number of decimals for the branch values can be specified by the user. }
\usage{ convert.to.newick(tree, decimals = 3) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ tree }{ One object of class \code{hclust} with the analyzis of the input motifs. The labels of the branches in the \code{hclust} tree will be the same in the Newick tree. }
  \item{ decimals }{ The number of decimals to round the distances. Default: 3. }
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

## Save the string with the Newick tree
newick.tree <- convert.hclust.to.newick(tree, decimals = 3)
}
\keyword{ Newick }
