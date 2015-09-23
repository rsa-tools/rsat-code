########################################################
## Given the id of a motif, return its name
get.name <- function(id){
  return(as.vector(sapply(id, function(X){
    global.description.table[which(global.description.table$id == X),][,"name"]
  })))
}