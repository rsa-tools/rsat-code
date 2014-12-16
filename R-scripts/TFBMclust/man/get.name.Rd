\name{get.name}
\alias{get.name}
\title{ Extract the name of a motif from a reference table. }
\description{
Given the ID of a (motif) leaf, search its corresponding motif name in the provided descrition table.

Note: The ids in the description table must be orderer alphabetically.
}
\usage{
get.name(id, desc.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name},     \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID. }

  \item{id}{ The ID of the motif. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{character} object with the motif name of the input motif ID.
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

## Get the ids of the nodes in level 2
ids.lev.2 <- get.id(leaves.per.level[[2]], oct4.tfbm.desc)

## Get the names of the nodes in level 2
## Note: In this example the names have the same values than the IDs
names.lev.2 <- get.name(ids.lev.2, oct4.tfbm.desc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
