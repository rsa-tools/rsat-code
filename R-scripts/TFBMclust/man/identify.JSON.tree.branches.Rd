\name{identify.JSON.tree.branches}
\alias{identify.JSON.tree.branches}
\title{ Identify the order of the branches in the JSON tree, relative to the order of one hclust tree. }
\description{ As the JSON tree order is not the same as the hclust tree, the JSON tree must be explored in order to rename the branches and thus can assign a name which will be used to add the branch consensuses. }
\usage{
identify.JSON.tree.branches(tree, desc.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name}, \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return an object of class \code{dataframe} with 3 columns:
  \item{level}{ The level in the JSON tree. }
  \item{cluster}{ The motifs of the IDs corresponding to that level. }
  \item{merge_level}{ The level in the hclust tree. }
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

## Load the data
data(oct4chipseq)

## Collect the motif distance matrix
oct4.dist <- build.distance.matrix(oct4.tfbm.compa, score = "Ncor")

## Save the matrix
oct4.dist.matrix <- oct4.dist$matrix

## Run the hierarchical clustering analysis
tree <- hclust.motifs(oct4.dist.matrix, hclust.method = "average")

##
JSON.clusters.table <- identify.JSON.tree.branches(tree, oct4.tfbm.desc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ JSON tree }
