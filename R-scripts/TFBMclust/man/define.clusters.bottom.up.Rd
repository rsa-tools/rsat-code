\name{define.clusters.bottom.up}
\alias{define.clusters.bottom.up}
\title{ Traverse the tree in a bottom-up approach to separate the motifs clusters. }
\description{ This function takes information from the status of the alignment (aligned or non-aligned) on each merge step of one object \code{tree$merge}, traverse the tree in a bottom-up approach to group the motfis in clusters. Return the IDs of the motifs on the defined cluster. }
\usage{
define.clusters.bottom.up(attributes.list, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{attributes.list}{ An object of class \code{list} with the information of each agglomeration step of the hierarchcial tree. This object is the output of the function align.motifs() . }
  \item{tree}{ One object of class \code{hclust} with the analyzis of the input motifs. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns an object of class \code{list} where each element correspond to the IDs of the cluster defined using the bottom-up approach. }
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
