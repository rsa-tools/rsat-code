fill.downstream <- function(motifs.ids, motifs.list){

  ## Get the highest length among the consensuses
  consensuses <- sapply(motifs.ids, function(X){
    nchar(motifs.list[[X]][["consensus"]])
  })
  max.width <- (max(unlist(consensuses)))

  ## Add the downstream spacer
  motifs.info.temp <- lapply(motifs.info[motifs.ids], function(x){

      ## Update the downstream spacer value
      x[["spacer.dw"]] <- max.width - nchar(x[["consensus"]])

      ## Update the consensus
      x[["consensus"]] <- paste(x[["consensus"]],
          paste(rep("-",
                    times = max.width - nchar(x[["consensus"]])),
                collapse = ""
          ),
          sep = "")
      return(x)
  })
  return(motifs.info.temp)
}