\name{find.clusters}
\alias{find.clusters}
\title{ Traverse the tree in a bottom-up approach to separate the motifs clusters. }
\description{ This function takes information from the status of the alignment (aligned or non-aligned) on each merge step of one object \code{tree$merge}, traverse the tree in a bottom-up approach to group the motfis in clusters. Return the IDs of the motifs on the defined cluster. }
\usage{
find.clusters(attributes.list, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{attributes.list}{ An object of class \code{list} with the information of each agglomeration step of the hierarchcial tree. This object is the output of the function align.motifs() . }
  \item{tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns an object of class \code{list} where each element correspond to the number of the cluster defined using the bottom-up approach. }
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

## Save the alignment attributes at each level of the tree
oct4.alignment.attributes <- oct4.alignment$node.attributes

## Find the clusters
clusters.oct4 <- find.clusters(oct4.alignment.attributes, tree)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Cluster }
\keyword{ Bottom-up }
