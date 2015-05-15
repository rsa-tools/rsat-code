\name{check.alignment}
\alias{check.alignment}
\title{ Using a given agglomeration method, checks if the clusters will be aligned. }
\description{ Given a list with the metrics names and its threshold vakues,  evaluates if the clusters will be aligned according an agglomeration method (average, complete, single) selected by the user. The function considers if the threshold metrics are correlation or distances. }
\usage{
check.alignment(id1, id2, compa.table, thresholds, hclust.method, hclust.metric)
}

\arguments{
  \item{id1}{ Cluster 1 IDs. }
  \item{id2}{ Cluster 2 IDs. }
  \item{compa.table}{ A comparison matrix containing at least 3 columns: \code{ID1}, \code{ID2} and one or more comparison scores between each pair of motifs. }
  \item{thresholds}{ An object of class \code{list} containing as fields (names) the metrics names and storing the values of the respective metrics. Default: Ncor<=0.4; cor<=0.5; w=5 . }
  \item{hclust.method}{ One agglomeration method used in the \code{hclust} function. Supported: average, complete, single. Default: average.

  \code{Average}: all the pairs of motifs in the two clusters are considered, for each metric in the threshold, the mean of all its values is calculated, the threshold value is evaluated on each metric, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned.

  \code{Complete}: using as reference the same metric used to build the hierarchical tree, select the farthest pair of motifs, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned.

  \code{Single}: using as reference the same metric used to build the hierarchical tree, select the closest pair of motifs, if the values of all metrics selected in the threshold are satisfied, the cluster can be aligned. }

  \item{hclust.metric}{ One metric to measure (dis)similarity between motifs. For consistence, use the same metric selected to build the hierarchical clustering with the function \code{hclust.motifs} }
  }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns \code{1} if the cluster should be aligned under the selected method.
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

## Get the IDs of the levels 10 and 12
IDs1 <- get.id(leaves.per.level[[10]], oct4.tfbm.desc)
IDs2 <- get.id(leaves.per.level[[12]], oct4.tfbm.desc)

## Check if these clusters should be aligned with the average method
check.alignment(IDs1, IDs2, compa.table, thresholds, hclust.method = "average", hclust.metric = "Ncor")

## Check if these clusters should be aligned with the complete method
check.alignment(IDs1, IDs2, compa.table, thresholds, hclust.method = "complete", hclust.metric = "Ncor")

## Check if these clusters should be aligned with the single method
check.alignment(IDs1, IDs2, compa.table, thresholds, hclust.method = "single", hclust.metric = "Ncor")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Average method }