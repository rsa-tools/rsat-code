##########################################################
## Given the number of a leaf, search the id of the
## corresponding motif in the provided descrition table
get.id <- function(num){
  return(as.vector(sapply(num, function(X){
    global.description.table[X,"id"]
  })))
}