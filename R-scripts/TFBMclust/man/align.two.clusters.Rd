\name{align.two.clusters}
\alias{align.two.clusters}
\title{ Align two already aligned clusters of motifs. }
\description{ This function take the information of the level of the tree and aligns the leaves corresponding to that level if both nodes are clusters of previously aligned motifs. The alignment is stored in the variable \code{motifs.info} which is inherited from the function \code{align.motifs()} where this function is called. }
\usage{
align.two.clusters(child1, child2, desc.table, compa.table, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), method = "average", metric = "Ncor", hclust.tree, nodes.attributes = TRUE, motif.at.tree.level = leaves.per.node(hclust.tree))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{child1}{ This number indicates the node 1 that is aligned, see the attribute \code{merge} for the \code{hclust} object, where a negative number means a single leaf and a positive number means a group of leaves. }
  \item{child2}{ This number indicates the node 2 that is aligned, see the attribute \code{merge} for the \code{hclust} object, where a negative number means a single leaf and a positive number means a group of leaves. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name},     \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID. }
  \item{compa.table}{ A comparison matrix containing at least 4 columns: \code{ID1}, \code{ID2}, \code{strand} and one or more comparison metrics between each pair of motifs. }
  \item{metric}{ Column of the comparison matrix to be used to build the hierarchical tree. Default: Ncor. }
  \item{thresholds}{ An object of class \code{list} containing as fields (\code{names()}) the metrics names and storing the values of the respective metrics. Default: Ncor <= 0.4; cor <= 0.5; w <= 5 . }
  \item{method}{ One agglomeration method used in the \code{hclust} function. Supported: average, complete, single. Default: average.

  \code{Average}: all the pairs of motifs in the two clusters are considered, for each metric in the threshold, the mean of all its values is calculated, the threshold value is evaluated on each metric, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned.

  \code{Complete}: using as reference the same metric used to build the hierarchical tree, select the farthest pair of motifs, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned.

  \code{Single}: using as reference the same metric used to build the hierarchical tree, select the closest pair of motifs, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned. }
  \item{metric}{ Column of the comparison matrix to be used as base to build the hierarchical clustering tree. Default: Ncor. }
  \item{hclust.tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
  \item{nodes.attributes}{ logical. Should be returned the alignment status of each level tree?. Default: FALSE. }
  \item{motif.at.tree.level}{ A \code{list} with the number of motifs at each level of the tree. This object correspond to the output of the function \code{levels.per.node} . }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the \code{list}s called \code{motifs.list} and \code{internal.nodes.attributes} updated, with the alignment of those levels containing a pair of aligned clusters, and its corresponding alignment status, respectively. }
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Alignment }