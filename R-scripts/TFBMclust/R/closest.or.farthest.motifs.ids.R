#########################################################################
## Identify the "central"/"distal" leaf (motif) of a subtree (cluster)
## i.e. the pair of motif with the smallest/largest distance to
## among all the other motifs of these clusters.
closest.or.farthest.motifs.ids <- function(id1, id2, compa.table, metric = "Ncor", closest = TRUE){

  ## Get the comparison number between the pairs of motifs of both clusters
  compa.numbers <- as.vector(get.comparison.number(id1, id2, compa.table))

  ## Get the ids of the less distant nodes

  #########################################
  ## For correlations (Higher is better)
  if ((metric == "Ncor")
      || (metric=="cor")
      || (metric=="logocor")
      || (metric=="Nlocogor")
      || (metric=="Icor")
      || (metric=="NIcor")
  ) {
    if(closest == TRUE){
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, metric] == max(compa.table[compa.numbers, metric])),c("id1","id2")][1,]
    }else{
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, metric] == min(compa.table[compa.numbers, metric])),c("id1","id2")][1,]
    }

  #################################
  ## For distances (Lower is better)
  }else if ((metric == "dEucl")
            || (metric == "NsEucl")
            || (metric == "SSD")
            || (metric == "SW")
            || (metric == "NSW")
  ){
    if(closest == TRUE){
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, metric] == min(compa.table[compa.numbers, metric])),c("id1","id2")][1,]
    } else {
      compa.info <- compa.table[compa.numbers,][which(compa.table[compa.numbers, metric] == max(compa.table[compa.numbers, metric])),c("id1","id2")][1,]
    }
  }
  return(as.vector(unlist(compa.info)))
}