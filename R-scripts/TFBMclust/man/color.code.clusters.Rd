\name{color.code.clusters}
\alias{color.code.clusters}
\title{ Assign one color to each cluster. }
\description{ Following the order of one object of class \code{hclust}, associate one number (1-4) corresponding to one color (black,red,blue,green). All the members of each cluster will be colored with the same color in the tree.This function produces a \code{vector} with number of 1 to 4. }
\usage{ color.code.clusters(clusters.ids, tree, desc.table) }
\arguments{
  \item{clusters.ids}{ One object of class \code{list} returned from the function define.clusters.bottom.up() where each elements is a \code{vector} with the Ids of each motif. }
  \item{tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name}, \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ An object of class \code{vector} with a serie of numbers with values 1-4  (one number for each motifs).
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
clusters.oct4 <- define.clusters.bottom.up(oct4.alignment.attributes, tree, oct4.tfbm.desc)

## Create the vector with the color code
color.code.clusters(clusters.oct4, tree, oct4.tfbm.desc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ color tree }
