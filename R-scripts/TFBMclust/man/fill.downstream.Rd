\name{fill.downstream}
\alias{fill.downstream}
\title{ Add gaps to the downstream end of the consensus }
\description{ This functions takes as input a set of IDs and the global list \code{motifs.info}, takes all the consensuses corresponding to the input IDs, search the largest consensus and add gaps ('-') in the downstream side of the consensuses. }
\usage{ fill.downstream(motifs.ids, motifs.list) }
\arguments{
  \item{motifs.ids}{ The IDs of the motifs.}
  \item{motifs.list}{ An object of class \code{list} containing one \code{list} for each motif analyzed. It must contain the next attributes: \code{name}, \code{consensus}, \code{strand}, \code{number}, \code{spacer.up} and \code{spacer.dw}. The \code{names()} of each list must be the motif \code{ID}. This object is inherited from the function \code{align.motifs()} which is the function where this function is called.  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the \code{motifs.list} with the attribute \code{consensus} filled with gaps in the downstream end. (If they are required) }
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

