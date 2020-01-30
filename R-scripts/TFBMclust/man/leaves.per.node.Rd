\name{leaves.per.node}
\alias{leaves.per.node}
\title{ Collect the nodes of each node of the tree. }
\description{ Collect the list of leaves (motifs) associated to each internal node (motif cluster) of the tree, using the attribute \code{merge} from the \code{hclust} object. The nodes are represented by numbers. By deafult \code{hclust} order the names of the object alphebetically and associate a number to this list.

This method works for any hclust result tree, irrespective of the
fact that it refers to motifs or anything else.
}
\usage{
leaves.per.node(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return an object of class \code{list} where each element contain the nodes of that merge level in the tree.
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

## Acces the leaves of level 2
leaves.per.level[[2]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hclust }
