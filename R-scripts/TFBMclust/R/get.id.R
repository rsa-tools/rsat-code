##########################################################
## Given the number of a leaf, search the id of the
## corresponding motif in the provided descrition table
get.id <- function(num, desc.table){
  return(as.vector(sapply(num, function(X){
    desc.table[X,"id"]
  })))
}