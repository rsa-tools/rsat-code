##########################################################
## Given the number of a leaf, search the id of the
## corresponding motif in the provided descrition table
get.consensus <- function(id, desc.table, RC = FALSE){

  return(as.vector(sapply(id, function(X){
    if(RC == FALSE){
      desc.table[which(desc.table$id == X),][,"consensus"]
    } else if(RC == TRUE){
      desc.table[which(desc.table$id == X),][,"rc.consensus"]
    }
  })))
}