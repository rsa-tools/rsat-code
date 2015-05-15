######################################################
## Given a list of IDs, the function fill up with
## gaps all the consensuses of the aligned motifs
## realtive to the largest aligned consensus + gaps
fill.downstream <- function(ids.fill, motifs.list){

  motifs.temp <- motifs.list[ids.fill]

  ## Get the highest length among the consensuses
  consensuses.d <- sapply(ids.fill, function(X){
    nchar(motifs.temp[[X]][["consensus_d"]])
  })
  max.width.d <- max(unlist(consensuses.d))

  ## Get the highest length among the consensuses
  consensuses.r <- sapply(ids.fill, function(X){
    nchar(motifs.temp[[X]][["consensus_rc"]])
  })
  max.width.r <- max(unlist(consensuses.r))

  ## Add the downstream spacer
  sapply(ids.fill, function(x){

    temp.list <- list()

    ## Update the D consensus
    cons.new.d <- NULL
    cons.new.d <- paste(motifs.temp[[x]][["consensus_d"]],
      paste(rep("-",
        times = max.width.d - nchar(motifs.temp[[x]][["consensus_d"]])),
        collapse = ""
      ),
      sep = ""
    )

    ## Update the R consensus
    cons.new.r <- NULL
    sp <- paste(rep("-", times = max.width.r - nchar(motifs.temp[[x]][["consensus_rc"]])), collapse = "")
    cons.new.r <- paste(sp, motifs.temp[[x]][["consensus_rc"]], sep = "")

    temp.list[[x]][["strand"]] <- motifs.temp[[x]][["strand"]]
    temp.list[[x]][["consensus_d"]] <- cons.new.d
    temp.list[[x]][["consensus_rc"]] <- cons.new.r
    temp.list[[x]][["name"]] <- motifs.temp[[x]][["name"]]
    temp.list[[x]][["number"]] <- motifs.temp[[x]][["number"]]
    temp.list[[x]][["spacer.up"]] <- get.spacer.nb(temp.list[[x]][["consensus_d"]])$up.spacer
    temp.list[[x]][["spacer.dw"]] <- get.spacer.nb(temp.list[[x]][["consensus_d"]])$dw.spacer

    ## Update the motifs list
    motifs.temp[names(temp.list)] <<- temp.list[names(temp.list)]
  })
  return(motifs.temp)
}