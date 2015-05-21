\name{draw.heatmap.motifs}
\alias{draw.heatmap.motifs}
\title{ Draw a heatmap of the input motifs, using the same distance tab and the same aggolmeration method.}
\description{
This function takes a (motifs) distance table and draws a heatmap, to be consistent, the agglomeration method must be the same used to obtain the hierarchical tree using the function \code{hclust.motifs}.
}
\usage{
draw.heatmap.motifs(dist.table, method = "average", clusters.list, alignment.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.table}{ A symetrical table containing the distances between each pair of the input motifs. }
  \item{method}{ Hierachical clustering agglomeration rule. Default: average. Supported: average, complete, single.}
  \item{clusters.list}{ An object of class \code{list} where the name of the elements are the cluster's name and the values are the motifs Ids belonging to each cluster. }
  \item{alignment.list}{ An object of class \code{list} containing one \code{list} for each motif analyzed. It must contain the next attributes: \code{name}, \code{consensus}, \code{strand}, \code{number}, \code{spacer.up} and \code{spacer.dw}. The \code{names()} of each list must be the motif \code{ID}.}
  \item{metric}{ Column of the comparison matrix to be used to build the hierarchical tree. Default: Ncor. In this function the color of the heatmap is based on the metric definition [correlation or distance]. }
    \item{tree.pos}{ The position in the heatmap where the hclust tree will be displayed. Default: coulumn. Supported: both, column, none, row.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ This function returns the heatmap corresponding to the distance table, where the column names are to the motifs Ids and the row names are the aligned consensus. In addition each cluster is indicated with a bar of a different color in the heatmap.

The scale of colors for the heatmap is adapted relative to the metric definition, some metrics measure correlation and others do distances.
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
oct4.dist <- build.distance.matrix(oct4.tfbm.compa, metric = "Ncor")

## Save the table
oct4.dist.table <- oct4.dist$table

## Run the hierarchical clustering analysis
tree <- hclust.motifs(oct4.dist.matrix, hclust.method = "average")

## Align the motifs
oct4.alignment <- align.motifs(tree, oct4.tfbm.desc, oct4.tfbm.compa, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), metric = "Ncor", method = "average", metric = "Ncor", nodes.attributes = TRUE, intermediate.alignments = TRUE)

## Save the alignment of all the motifs
oct4.alignment.motifs <- oct4.alignment$motifs.alignment

## Save the alignment attributes at each level of the tree
oct4.alignment.attributes <- oct4.alignment$node.attributes

## Find the clusters
clusters.oct4 <- define.clusters.bottom.up(oct4.alignment.attributes, tree, oct4.tfbm.desc)

## Draw the heatmap
draw.heatmap.motifs(oct4.dist.table, method = "average", clusters.oct4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Heatmap }