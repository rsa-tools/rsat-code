#' @title Draw a simple heatmap from a matrix, with values written in the cells.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description  Draw a simple heatmap from a matrix, without any clustering or 
#' fancy features, but with the values written in the cells. This is a simple way to circumvent 
#' several difficulties with the heatmap() and heatmap.2() functions, which 
#' generate beautiful maps with many intersting features, which I generally spend
#' a lot of time to disactivate. 
#'
#' @details
#' First version: 2015-10
#' Last modification: 2015-10
#'
#' @param x A matrix or data frame
#' @param main=NA Main title, passed to image()
#' @param xlab=NA Label for X axis, passed to image()
#' @param ylab=NA Label for Y axis, passed to image()
#' @param display.values=TRUE Display the values on the cell.
#' @param text.color="black" Color to write the values.
#' @param col=gray.colors(256, start=1, end=0) Color map. 
#' @param zlim=range(m) passed to image()
#' @param round.digits=NA If specified, values are rounded with the function round().
#' @param signif.digits=NA If specified, values are rounded with the function signif().
#' @param ...  Additional parameters are passed to the function image()
#' 
#' @return No return
#' @examples
#' ## Draw the heatmap of a Markov model of order 1 trained with the
#' ##  whole Human genome.
#' ## This example highlights the strong avoidance of the CpG dinucleottide.
#' org <- "Homo_sapiens_GRCh38"
#' m <- 1 ## Order of the Markov model
#' k <- m + 1 ## Length of oligomers used to estimate the parameters
#' transition.file <- system.file("extdata", 
#'   file.path("markov_models", org, 
#'       paste(sep="", k, "nt_genomic_", org, "-1str-ovlp_transitions.tab.gz")), 
#'       package = "stats4bioinfo")
#' markov.transitions <- read.delim(transition.file, comment=";", row.names=1)
#' # markov.transitions <- round(digits=2, as.matrix(markov.transitions))
#' heatmap.simple(markov.transitions[, c("a","c","g","t")], zlim=c(0,1),
#'    main=paste(org, "genomic Markov, order =", m), signif.digits=10,
#'    xlab="Suffix", ylab="Prefix")
#' @export
heatmap.simple <- function(x,
                           main=NA,
                           xlab=NA,
                           ylab=NA,
                           display.values=TRUE,
                           text.color="auto",
                           col=gray.colors(256, start=1, end=0, gamma=1),
                           zlim=range(x, na.rm=TRUE),
                           round.digits=NA,
                           signif.digits=NA,
                           ...) {
  
  ## Invert the input matrix to display heatmap from topleft corner
  x <- as.matrix(x[nrow(x):1,])

  ## Rounding by significant digits
  if (!is.na(signif.digits)) {
    x <- signif(digits=signif.digits, x)
  }
  
  ## Rounding by fixed number of digits
  if (!is.na(round.digits)) {
    x <- round(digits=round.digits, x)
  }
  
  ## Compute margins
  max.label.len <- max(nchar(rownames(x)))
  par(mar=c(max.label.len*0.5 + 1.1, max.label.len*0.5+1.1,4.1,1.1))
  par(las=2)
  
  ## Generate the heat map
  image(1:ncol(x), 1:nrow(x), t(x), 
        main=main,
        xlab=xlab,
        ylab=ylab,
        col=col,
        zlim=zlim,
        axes = FALSE, ...)
  
  ## Row and columns labels
  axis(1, 1:ncol(x), colnames(x), tick=FALSE)
  axis(2, 1:nrow(x), rownames(x), las=2, tick=FALSE)
  
  ## Display values
  if (display.values) {
    for (i in 1:ncol(x))
      for (j in 1:nrow(x))
        if (sum(dim(as.data.frame(text.color)) == dim(x))==2) {
          ## Cell-specific text colors defined inby passing a text color table
          text(i, j, x[j,i], col=text.color[j,i])
        } else if (text.color == "auto") {
          ## Automatic selection of a color contrasting with the heatmap color
          if ((x[j,i] <= zlim[2]*0.5) || (x[j,i] > zlim[2])) {
            text(i, j, x[j,i], col=col[length(col)])
          } else {
            text(i, j, x[j,i], col=col[1])
          }
        } else {
          ## Same color for all cells, specified in the attributes
          text(i, j, x[j,i], col=text.color)
        }
  }
}
