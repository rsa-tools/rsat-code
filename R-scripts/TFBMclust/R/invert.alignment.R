########################################################
## Given a vector with IDs, return a list with the
## information (consensus, number, id, strand, spacer)
## of the inverted alignment
inverted.alignment <- function(ids, motifs.list, desc.table){

  temporal.list <- list()

  ## Get the higher length among the consensuses
  consensuses <- sapply(motifs.info[ids], function(X){
    nchar(X[["consensus"]])
  })
  max.cons.length <- max(consensuses)

  ## Create the new consensuses with the inverted orientation
  temporal.list <- sapply(ids, function(X){

    inverted.aligment.list <- list()
    inverted.consensus <- NULL
    spacer.length <- NULL
    spacer <- NULL

    ## Create the spacer
    spacer.length <- as.numeric(motifs.list[[X]][["spacer"]])
    spacer <- paste(rep("-", times = spacer.length), collapse = "")

    ## Get the consensus in the inverted orientation
    if(motifs.info[[X]][["strand"]] == "D"){
      inverted.consensus <- as.vector(description.table[[as.numeric(motifs.info[[X]][["number"]]), "rc_consensus"]])
      inverted.aligment.list[[X]][["strand"]] <- "R"
    }else{
      inverted.consensus <- as.vector(description.table[[as.numeric(motifs.info[[X]][["number"]]), "consensus"]])
      inverted.aligment.list[[X]][["strand"]] <- "D"
    }

    ## Fill the spaces with "-"
    new.consensus <- paste(inverted.consensus, spacer, sep = "")
    new.consensus <- paste(paste(rep("-", times = max.cons.length - nchar(new.consensus)), collapse = ""), new.consensus, sep = "")
    new.consensus <- paste(unlist(strsplit(new.consensus, "[A-Za-z]+", perl = TRUE))[1], inverted.consensus, sep = "")

    ## Fill the new list
    inverted.aligment.list[[X]][["consensus"]] <- new.consensus
    inverted.aligment.list[[X]][["number"]] <- as.numeric(motifs.info[[X]][["number"]])
    inverted.aligment.list[[X]][["name"]] <- motifs.info[[X]][["name"]]
    inverted.aligment.list[[X]][["spacer"]] <- length(unlist(strsplit(new.consensus, "-")))-1

    return(inverted.aligment.list)
  })

  #names(temporal.list) <- sapply(strsplit(names(temporal.list), "\\."), function(X){return(X[1])})
  names(temporal.list) <- ids
  return(temporal.list)
}