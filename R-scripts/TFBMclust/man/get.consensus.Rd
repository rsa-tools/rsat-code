\name{get.consensus}
\alias{get.consensus}
\title{ Extract the motif IUPAC consensus (or its reverse complement) from a reference table. }
\description{
Given a motif ID, search in the reference table its corresponding IUPAC consensus or its reverse complement.
}
\usage{
get.consensus(id, desc.table, RC = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{ The ID of the motif. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name},     \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID. }
  \item{RC}{ logical.  Should be searched the consensus or its Reverse Complement?. Default: FALSE. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a \code{character} object with the motif consensus (or its RC) of the input motif ID. }
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

## Get the consensus of the nodes in level 2
get.consensus(ids.lev.2, oct4.tfbm.desc, RC = FALSE)

## Get the RC consensus of the nodes in level 2
get.consensus(ids.lev.2, oct4.tfbm.desc, RC = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IUPAC consensus }
