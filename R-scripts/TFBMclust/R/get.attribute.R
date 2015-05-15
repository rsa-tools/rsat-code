##############################################################
## Given a motif number, returns a selected attribute
## Supported attributes: id, name, consensus_d, consensus_r
get.attribute <- function(n, desc.table, attribute = "id"){

  ## Get the motifs ID
  id <- as.vector(sapply(n, function(X){
    desc.table[X,"id"]
  }))

  ## Return the IDs
  if(attribute == "id"){
    return(id)

  ## Get the motif name
  } else if(attribute == "name"){
    return(as.vector(sapply(id, function(X){
      desc.table[which(desc.table$id == X),][,"name"]
    })))

  ## Get the motif direct consensus
  } else if(attribute == "consensus_d"){
    return(as.vector(sapply(id, function(X){
        desc.table[which(desc.table$id == X),][,"consensus"]
    })))

  ## Get the motif reverse consensus
  } else if(attribute == "consensus_r"){
    return(as.vector(sapply(id, function(X){
        desc.table[which(desc.table$id == X),][,"rc.consensus"]
    })))

  ## Other attributes are not allowed
  } else{
    stop("Non valid attribute. Supported: id, name, consensus_d, consensus_r .")
  }
}