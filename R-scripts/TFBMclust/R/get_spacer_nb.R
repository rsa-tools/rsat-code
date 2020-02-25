##############################################################
## Given a consensus string, count the number of gaps ("+")
## in the upstream and downstream ends. Return a list with
## these values.
##
## Originally, the gaps were represented with '-', however, this was not optimal when
## there are empty positions in the input matrices (e.g., assemblies), because the motif
## consensus in such empty positions was annotated also with '-'.
## This caused a problem to shift the alignment when one branch was aligned with a single motif
## and when two branches were aligned.
## To avoid such problem from now on, the gaps are represented with '+'
##
## Updated by Jaime Castro (25-02-2020)
get.spacer.nb <- function(consensus){

  ##
  consensus.mod <- gsub(consensus, pattern = "[a-zA-z-]", replacement = "A")
  s.up <- NULL
  s.dw <- NULL

  ## Update the downstream spacer value
  ## Check if the consensus contains gaps ('+')
  if( grepl("\\+", consensus.mod) ) {

    ## Check if contains gaps in both ends, if so calculate the downstream offset
    if( grepl("^\\+", consensus.mod, perl = TRUE) & grepl("\\+$", consensus.mod, perl = TRUE) ) {
      s.up <- nchar(unlist(strsplit(consensus.mod, "\\w+"))[1])
      s.dw <- nchar(unlist(strsplit(consensus.mod, "\\w+"))[2])
    }

    ## If there are gaps in the upstream but not in the downstream end, the downstream offset is zero
    if( (grepl("^\\+", consensus.mod, perl = TRUE)) & !(grepl("\\+$", consensus.mod, perl = TRUE))){
      s.up <- nchar(unlist(strsplit(consensus.mod, "\\w+"))[1])
      s.dw <-  0
    }

    ## If there are gaps in the upstream end, but not in the downstream end, the downstream offset is zero
    if( !(grepl("^\\+", consensus.mod, perl = TRUE)) & (grepl("\\+$", consensus.mod, perl = TRUE))){
      s.up <-  0
      s.dw <- nchar(unlist(strsplit(consensus.mod, "\\w+"))[2])
    }
  } else{
    s.up <-  0
    s.dw <-  0
  }

  return(list(up.spacer = s.up,
              dw.spacer = s.dw))
}


## OLD version
# get.spacer.nb <- function(consensus){
#
#   # consensus.mod <- gsub(consensus, pattern = "[a-zA-z-]", replacement = "A")
#
#   ## Update the downstream spacer value
#   ## Check if the consensus contains gaps ('+')
#   if(grepl("\\+", consensus)){
#
#     ## Check if contains gaps in both ends, if so calculate the downstream offset
#     if(grepl("^\\+", consensus, perl = TRUE) & grepl("\\+$", consensus, perl = TRUE)){
#       s.up <- nchar(unlist(strsplit(consensus, "\\w+"))[1])
#       s.dw <- nchar(unlist(strsplit(consensus, "\\w+"))[2])
#     }
#
#     ## If there are gaps in the upstream but not in the downstream end, the downstream offset is zero
#     if( (grepl("^\\+", consensus, perl = TRUE)) & !(grepl("\\+$", consensus, perl = TRUE))){
#       s.up <- nchar(unlist(strsplit(consensus, "\\w+"))[1])
#       s.dw <-  0
#     }
#
#     ## If there are gaps in the upstream end, but not in the downstream end, the downstream offset is zero
#     if( !(grepl("^\\+", consensus, perl = TRUE)) & (grepl("\\+$", consensus, perl = TRUE))){
#       s.up <-  0
#       s.dw <- nchar(unlist(strsplit(consensus, "\\w+"))[2])
#     }
#   } else{
#     s.up <-  0
#     s.dw <-  0
#   }
#   return(list(up.spacer = s.up, dw.spacer = s.dw))
# }