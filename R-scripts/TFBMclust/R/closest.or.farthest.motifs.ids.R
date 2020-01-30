#########################################################################
## Identify the "central"/"distal" leaf (motif) of a subtree (cluster)
## i.e. the pair of motif with the smallest/largest distance to
## among all the other motifs of these clusters.
closest.or.farthest.motifs.ids <- function(id1, id2, metric = "Ncor", closest = TRUE){

  ## Get the comparison number between the pairs of motifs of both clusters
  compa.numbers <- as.vector(get.comparison.number(id1, id2))

  ## Get the ids of the less distant nodes

  #########################################
  ## For correlations (Higher is better)
  if ((metric == "Ncor")
      || (metric=="cor")
      || (metric=="NcorS")
      || (metric=="logocor")
      || (metric=="Nlogocor")
      || (metric=="Icor")
      || (metric=="NIcor")
      || (metric=="logoDP")
      || (metric == "mean_zscore")
  ) {
    if(closest == TRUE){
      compa.info <- global.compare.matrices.table[compa.numbers,][which(global.compare.matrices.table[compa.numbers, metric] == max(global.compare.matrices.table[compa.numbers, metric])),c("id1","id2")][1,]
    }else{
      compa.info <- global.compare.matrices.table[compa.numbers,][which(global.compare.matrices.table[compa.numbers, metric] == min(global.compare.matrices.table[compa.numbers, metric])),c("id1","id2")][1,]
    }

  #################################
  ## For distances (Lower is better)
  }else if ((metric == "dEucl")
            || (metric == "NdEucl")
            || (metric == "SSD")
            || (metric == "SW")
            || (metric == "NSW")
            || (metric == "rank_mean")
  ){
    if(closest == TRUE){
      compa.info <- global.compare.matrices.table[compa.numbers,][which(global.compare.matrices.table[compa.numbers, metric] == min(global.compare.matrices.table[compa.numbers, metric])),c("id1","id2")][1,]
    } else {
      compa.info <- global.compare.matrices.table[compa.numbers,][which(global.compare.matrices.table[compa.numbers, metric] == max(global.compare.matrices.table[compa.numbers, metric])),c("id1","id2")][1,]
    }
  }
  return(as.vector(unlist(compa.info)))
}