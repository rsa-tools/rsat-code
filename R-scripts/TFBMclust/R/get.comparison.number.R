##############################################################
## Given two sets of Ids, return the comparison numbers for
## all the pairs of the leaves belonging the two clusters
get.comparison.number <- function(id1, id2, compa.table){

  compa.numbers <- sapply(id2, function(X){
    sapply(id1, function(Y){
      which( (compa.table[,"id2"] == X & compa.table[,"id1"] == Y) | (compa.table[,"id1"] == X & compa.table[,"id2"] == Y) )
    })
  })
  compa.numbers <- as.vector(compa.numbers)
  return(compa.numbers)
}