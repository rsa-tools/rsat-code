############################################################################
## Given two sets of IDs (clusters), return all the comparison numbers
## between all the pairs of the leaves corresponding to the two clusters
get.comparison.number <- function(id1, id2){

  compa.numbers <- sapply(id2, function(X){
    sapply(id1, function(Y){
      which( (global.compare.matrices.table[,"id2"] == X & global.compare.matrices.table[,"id1"] == Y) |
               (global.compare.matrices.table[,"id1"] == X & global.compare.matrices.table[,"id2"] == Y) )[1]
    })
  })
  compa.numbers <- as.vector(compa.numbers)
  return(compa.numbers)
}