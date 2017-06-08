#' @title Given two tables (formatted as data frames or matrices), apply a chi2 homogeneity test to each pair of rows.
#' @author Jacques van Helden
#' @description Given two data frames or matrices, apply a chi2 homogeneity test to each 
#' pair of rows.
#' @param x1 first table of counts. Must be absolute counts and not relative frequencies. 
#' @param x2 second vector of counts. Must be absolute counts and not relative frequencies. 
#' @examples 
#'
#' ##  TEMPORARY: we use some pre-loaded dataset from Marianne + Jaime
#' ## This will have to be replaced by some example of count tables in RSAT demo
#' x1 <- counts.per.bin.query
#' x2 <- counts.per.bin.control
#' 
#' ## Run the chi2 test of homogeneity on each row of a pair of count tables
#' stats.per.row <- Chi2TestPerRow(x1, x2)
#'        
#' ## Run the chi2 test of homogeneity on all the rows (motifs) of the count tables. 
#'    
#' @export
Chi2TestPerRow <- function(x1, x2) {

  ## Check parameter validity
  if (sum(dim(x1) != dim(x2)) > 0) {
    stop("Chi2TestPerRow() fatal error: ",
         "dimensions differ between x1 (", paste(dim(x1), collapse=","), ") ",
      "and x2 (", paste(dim(x2), collapse=","), ").")
  }
  
  ## Compute matrix dimensions
  rows <- nrow(x1)
  columns <- ncol(x1)

  
  ## Build a 3D array with observed counts
  obs.counts <- array(c(x1, x2), c(rows, columns, 2))
#  t(obs.counts[1,,]) %*%
  
    
  ## Compute column (bin) sums for each row and layer
  xzsum <- data.frame(apply(obs.counts, c(1,3), sum))
  rownames(xzsum) <- rownames(x1)
  colnames(xzsum) <- c("x1", "x2")
  # dim(xzsum)
  # View(xzsum)
  
  ## Compute layer sum for each row and column (bin)
  xysum <- data.frame(apply(obs.counts, c(1,2), sum))
  rownames(xysum) <- rownames(x1)
  colnames(xysum) <- colnames(x1)
  # dim(xysum)
  
  xsum <- apply(obs.counts, 1, sum)
  
  
  ## Compute expected occurrences
  count.sum <- sum(obs.counts)
  exp.counts <- array(NA, c(rows, columns, 2))
  for (k in c(1,2)) {
    exp.counts[,,k] <- as.matrix(xysum * xzsum[,k]) / xsum
  }
  View(exp.counts)
  
  ## Sum check
  # sum(obs.counts) == sum(xzsum)
  # sum(obs.counts) == sum(xysum)
  # sum(obs.counts) == round(sum(exp.counts))
  
  ## Compute a 3D array with the chi2 elements
  chi2.array <- (obs.counts - exp.counts)^2 / exp.counts
  
  result <- list()
  result$chi2.stat <- apply(chi2.array,1,sum)
  result$pval.log10 <- pchisq(result$chi2.stat, df=columns - 1, lower.tail = FALSE, log=TRUE)/log(10)
  result$pval <- 10 ^ result$pval.log10
  result$eval <- result$pval * rows
  result$sigbis <- -log10(result$eval)
  result$sig <- - (result$pval.log10 + log10(rows))
  # View(result)
  return(result)
}


