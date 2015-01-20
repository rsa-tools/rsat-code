#################################################################
## Identify the "central" leaf (motif) of a subtree (cluster),
## i.e. the pair of motif with the smallest average distance to
## among all the other motifs of these clusters.
closest.or.farthest.motifs.ids <- function(id1, id2, compa.table, score = "Ncor", closest = TRUE){

  ## Get the comparison number between the pairs of motifs of both clusters
  compa.numbers <- as.vector(get.comparison.number(id1, id2, compa.table))

  ## Get the ids of the less distant nodes

  #########################################
  ## For correlations (Higher is better)
  if ((score == "Ncor")
      || (score=="cor")
      || (score=="logocor")
      || (score=="Nlocogor")
      || (score=="Icor")
      || (score=="NIcor")
  ) {
    if(closest == TRUE){
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, score] == max(compa.table[compa.numbers, score])),c("id1","id2")][1,]
    }else{
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, score] == min(compa.table[compa.numbers, score])),c("id1","id2")][1,]
    }

  #################################
  ## For distances (Lower is better)
  }else if ((score == "dEucl")
            || (score == "NsEucl")
            || (score == "SSD")
            || (score == "SW")
            || (score == "NSW")
  ){
    if(closest == TRUE){
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, score] == min(compa.table[compa.numbers, score])),c("id1","id2")][1,]
    } else {
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, score] == max(compa.table[compa.numbers, score])),c("id1","id2")][1,]
    }
  }
  return(as.vector(unlist(compa.info)))
}