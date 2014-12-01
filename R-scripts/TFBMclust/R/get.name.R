########################################################
## Given the id of a motif, return its name
get.name <- function(id, desc.table){
  return(as.vector(sapply(id, function(X){
    desc.table[which(desc.table$id == X),][,"name"]
  })))
}