\name{align.motifs}
\alias{align.motifs}
\title{ Align all the motifs, using as base the merge order of one hierarchical tree. }
\description{ This function does a progressive alignment of the motifs following the order of a hierarchical clustering tree. Identifies if the aligned motis are: two leaves, one leaf and one cluster or two clusters. Export a \code{list} with the alignment of each motifs, if it is indicated returns a \code{list} with the alignment status on each level of the tree. Also, if it is indicated, return the information of each intermediate alignment at each level of the tree. }
\usage{
align.motifs(hclust.tree, desc.table, compa.table, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", nodes.attributes = TRUE, intermediate.alignments = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hclust.tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name},     \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID. }
  \item{compa.table}{ A comparison matrix containing at least 4 columns: \code{ID1}, \code{ID2}, \code{strand} and one or more comparison metrics between each pair of motifs. }
  \item{thresholds}{ An object of class \code{list} containing as fields (\code{names()}) the metrics names and storing the values of the respective metrics. Default: Ncor <= 0.4; cor <= 0.5; w <= 5 . }
  \item{method}{ One agglomeration method used in the \code{hclust} function. Supported: average, complete, single. Default: average.

  \code{Average}: all the pairs of motifs in the two clusters are considered, for each metric in the threshold, the mean of all its values is calculated, the threshold value is evaluated on each metric, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned.

  \code{Complete}: using as reference the same metric used to build the hierarchical tree, select the farthest pair of motifs, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned.

  \code{Single}: using as reference the same metric used to build the hierarchical tree, select the closest pair of motifs, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned. }
  \item{metric}{ Column of the comparison matrix to be used as base to build the hierarchical clustering tree. Default: Ncor. }
  \item{nodes.attributes}{ logical. Should be returned the alignment status of each level tree?. Default: FALSE. }
  \item{intermediate.alignments}{ logical. Should be returned the intermediate alignment corresponding to each level of the tree?. Default: FALSE. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{list} with three elements:
  \item{motifs.alignment}{ An object of class \code{list} where each element correspond to the motifs Id and each value correspond to information of the motif: name, strand, number, spacer ustream and spacer downstream. }
  \item{node.attributes}{ An object of class \code{list} where each element correspond to the level of the tree and each value correspond to information of the alignment in that level: ids aligned, alignmet status. }
  \item{intermediate.alignments}{ An object of class \code{list} where each element correspond to the level of the tree, with one \code{list} for each motif aligned in that level, each value of the motif list correspond to information of the motif: name, strand, number, spacer ustream and spacer downstream.}
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

## Save the matrix
oct4.dist.matrix <- oct4.dist$matrix

## Run the hierarchical clustering analysis
tree <- hclust.motifs(oct4.dist.matrix, hclust.method = "average")

## Align the motifs
oct4.alignment <- align.motifs(tree, oct4.tfbm.desc, oct4.tfbm.compa, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", nodes.attributes = TRUE, intermediate.alignments = TRUE)

## Save the alignment of all the motifs
oct4.alignment.motifs <- oct4.alignment$motifs.alignment

## Save the alignment attributes at each level of the tree
oct4.alignment.attributes <- oct4.alignment$node.attributes

## Save the intermediate alignments at each level of the tree
oct4.alignment.intermediates <- oct4.alignment$intermediate.alignments

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Alignment }
