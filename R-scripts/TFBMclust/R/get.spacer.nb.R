##############################################################
## Given a consensus string, count the number of gaps ("-")
## in the upstream and downstream ends. Return a list with
## these values.
get.spacer.nb <- function(consensus){

  ## Update the downstream spacer value
  ## Check if the consensus contains gaps ('-')
  if(grepl("-", consensus)){

    ## Check if contains gaps in both ends, if so calculate the downstream offset
    if(grepl("^-", consensus, perl = TRUE) & grepl("-$", consensus, perl = TRUE)){
        spacer.up <- nchar(unlist(strsplit(consensus, "\\w+"))[1])
        spacer.dw <- nchar(unlist(strsplit(consensus, "\\w+"))[2])
    }

    ## If there are gaps in the upstream but not in the downstream end, the downstream offset is zero
    if( (grepl("^-", consensus, perl = TRUE)) & !(grepl("-$", consensus, perl = TRUE))){
      spacer.up <- nchar(unlist(strsplit(consensus, "\\w+"))[1])
      spacer.dw <-  0
    }

    ## If there are gaps in the upstream end, but not in the downstream end, the downstream offset is zero
    if( !(grepl("^-", consensus, perl = TRUE)) & (grepl("-$", consensus, perl = TRUE))){
      spacer.up <-  0
      spacer.dw <- nchar(unlist(strsplit(consensus, "\\w+"))[2])
    }
  } else{
    spacer.up <-  0
    spacer.dw <-  0
  }
  return(list(up.spacer = spacer.up, dw.spacer = spacer.dw))
}