fill.downstream <- function(ids.fill, motifs.list){

  motifs.temp <- motifs.list[ids.fill]

  ## Get the highest length among the consensuses
  consensuses <- sapply(ids.fill, function(X){
    nchar(motifs.temp[[X]][["consensus"]])
  })
  max.width <- max(unlist(consensuses))

  ## Add the downstream spacer
  sapply(ids.fill, function(x){

    temp.list <- list()
    ## Update the consensus
    cons.new <- NULL
    cons.new <- paste(motifs.temp[[x]][["consensus"]],
      paste(rep("-",
        times = max.width - nchar(motifs.temp[[x]][["consensus"]])),
        collapse = ""
      ),
      sep = ""
    )

    temp.list[[x]][["strand"]] <- motifs.temp[[x]][["strand"]]
    temp.list[[x]][["consensus"]] <- cons.new
    temp.list[[x]][["name"]] <- motifs.temp[[x]][["name"]]
    temp.list[[x]][["number"]] <- motifs.temp[[x]][["number"]]
    temp.list[[x]][["spacer.up"]] <- get.spacer.nb(temp.list[[x]][["consensus"]])$up.spacer
    temp.list[[x]][["spacer.dw"]] <- get.spacer.nb(temp.list[[x]][["consensus"]])$dw.spacer

    ## Update the motifs list
    motifs.temp[names(temp.list)] <<- temp.list[names(temp.list)]
  })
  #print(motifs.temp)
  return(motifs.temp)
}