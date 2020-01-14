##########################################################
## Given the number of a leaf, search the id of the
## corresponding motif in the provided descrition table
get.consensus <- function(id, RC = FALSE){

  return(as.vector(sapply(id, function(X){
    if(RC == FALSE){
      global.description.table[which(global.description.table$id == X),][,"consensus"]
    } else if(RC == TRUE){
      global.description.table[which(global.description.table$id == X),][,"rc.consensus"]
    }
  })))
}