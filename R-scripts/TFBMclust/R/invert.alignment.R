########################################################
## Given a vector with IDs, return a list with the
## information (consensus, number, id, strand, spacer)
## of the inverted alignment
inverted.alignment <- function(ids.inv, motifs.list){

  motifs.list.temp <- motifs.list[ids.inv]
  ## Create the new consensuses with the inverted orientation
  sapply(ids.inv, function(X){

    inverted.aligment.list <- list()
    inverted.consensus <- NULL
    spacer.up.length <- NULL
    spacer.up <- NULL
    spacer.dw.length <- NULL
    spacer.dw <- NULL

    ## Create the spacers
    spacer.up.length <- get.spacer.nb(motifs.list.temp[[X]][["consensus_d"]])$up.spacer
    spacer.dw.length <- get.spacer.nb(motifs.list.temp[[X]][["consensus_d"]])$dw.spacer

    spacer.up <- paste(rep("-", times = spacer.up.length), collapse = "")
    spacer.dw <- paste(rep("-", times = spacer.dw.length), collapse = "")

    inverted.consensus.1 <- NULL
    inverted.consensus.2 <- NULL

    ## Depending on the current orientation, get the consensus in the inverted orientation
    if(motifs.list.temp[[X]][["strand"]] == "D"){
      inverted.consensus.1 <- get.consensus(X, RC = TRUE)
      inverted.consensus.2 <- get.consensus(X, RC = FALSE)
      inverted.aligment.list[[X]][["strand"]] <- "R"
    } else{
      inverted.consensus.1 <- get.consensus(X, RC = FALSE)
      inverted.consensus.2 <- get.consensus(X, RC = TRUE)
      inverted.aligment.list[[X]][["strand"]] <- "D"
    }

    ## Create the new consensus
    new.consensus.d <- paste(spacer.dw, inverted.consensus.1, spacer.up, sep = "")
    new.consensus.r <- paste(spacer.up, inverted.consensus.2, spacer.dw, sep = "")

    ## Fill the new list
    inverted.aligment.list[[X]][["consensus_d"]] <- new.consensus.d
    inverted.aligment.list[[X]][["consensus_rc"]] <- new.consensus.r
    inverted.aligment.list[[X]][["name"]] <- motifs.list.temp[[X]][["name"]]
    inverted.aligment.list[[X]][["number"]] <- as.numeric(motifs.list.temp[[X]][["number"]])
    inverted.aligment.list[[X]][["spacer.up"]] <- get.spacer.nb(new.consensus.d)$up.spacer
    inverted.aligment.list[[X]][["spacer.dw"]] <- get.spacer.nb(new.consensus.d)$dw.spacer

    motifs.list.temp[names(inverted.aligment.list)] <<- inverted.aligment.list[names(inverted.aligment.list)]
  })
  return(motifs.list.temp)
}