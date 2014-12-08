\name{get.spacer.nb}
\alias{get.spacer.nb}
\title{ Counts the number of gaps ("-") in both ends of a consensus string. }
\description{ Given an object of class \code{character} (a motif consensus), this function counts the number of gaps ("-") in both ends of the sequence (upstream and downstream) and return an objecto of class \code{list} with these values.}
\usage{ get.spacer.nb(consensus) }
\arguments{
  \item{consensus}{ An object of class \code{character}. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with two elements
\item{up.spacer}{ An object of class \code{numeric} indicating the number of gaps in the upstream end of the consensus}.
\item{dw.spacer}{ An object of class \code{numeric} indicating the number of gaps in the downstream end of the consensus}.
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
## Declare a motif consensus string
consensus <- "-----bwTwTKYWwATgh--"

## Count the number of gaps in both ends
get.spacer.nb(consensus)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Consensus }
