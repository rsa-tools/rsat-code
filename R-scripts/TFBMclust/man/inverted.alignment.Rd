\name{inverted.alignment}
\alias{inverted.alignment}
\title{ Invert the aligned consensus and orientation information for a group of aligned motifs. }
\description{ This function takes as input a set of motifs IDs and search these motifs in the motifs.list, takes the consensus and takes from the desc.table the consensus on its opposite orientation, then add the gaps ("-") and finally update the new information. The modified IDs are exported in a \code{list} with the same attributes than the \code{list} object \code{motifs.info}. }
\usage{ inverted.alignment(ids, motifs.list, desc.table) }
\arguments{
  \item{ids}{ A \code{vector} with the ID of the motif(s). }
  \item{motifs.list}{ An object of class \code{list} containing one \code{list} for each motif analyzed. It must contain the next attributes: \code{name}, \code{consensus}, \code{strand}, \code{number}, \code{spacer.up} and \code{spacer.dw}. The \code{names()} of each list must be the motif \code{ID}. This object is inherited from the function \code{align.motifs()} which is the function where this function is called. }
  \item{desc.table}{ A table with the motifs attributes with at least 4 columns: \code{id}, \code{name},     \code{consensus} and \code{rc.consensus}.

  The ID must be unique, the name can be the same as the ID.} }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a \code{list} object with the same attributes than \code{motifs.list} with the inverted consensus and orientation. }
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Inverted Alignment }
