#####################################################################
## Evaluate if the mean of scores for all pairs of leaves between
## the two clusters is upper/lower than the threshold,
## if so, then the two clusters should be aligned
check.alignment <- function(id1,
                            id2,
                            compa.table,
                            thresholds = list(Ncor = 0.4, cor = 0.6, w = 5),
                            hclust.method = "average",
                            metric = "Ncor"){

  ## Hclust method = average
  ## Calculate the mean of the scores for all the pairs of motifs
  if(hclust.method == "average"){
    compa.numbers <- get.comparison.number(id1, id2, compa.table)

    scores <- compa.table[compa.numbers, names(thresholds)]

    ## Calculate the median of the data
    mean.scores <- apply(scores, 2, mean)

  ## Hclust method = complete
  ## Calculate the farthest couple of motifs between all the pairs of motifs
  } else if(hclust.method == "complete"){

    farthest.motifs <- closest.or.farthest.motifs.ids(id1,
                                                      id2,
                                                      compa.table,
                                                      metric = metric,
                                                      closest = FALSE)
    id1.far <- farthest.motifs[1]
    id2.far <- farthest.motifs[2]

    compa.numbers <- get.comparison.number(id1.far, id2.far, compa.table)

    ## Get the scores of the comparisons
    farthest.scores <- compa.table[compa.numbers, names(thresholds)]

  ## Hclust method = single
  ## Calculate the closest pair of motifs between all the pairs of motifs
  } else if(hclust.method == "single"){

    closest.motifs <- closest.or.farthest.motifs.ids(id1,
                                                     id2,
                                                     compa.table,
                                                     metric = metric,
                                                     closest = TRUE)
    id1.close <- closest.motifs[1]
    id2.close <- closest.motifs[2]

    compa.numbers <- get.comparison.number(id1.close, id2.close, compa.table)

    ## Get the scores of the comparisons
    closest.scores <- compa.table[compa.numbers, names(thresholds)]
  }

  th <- list()

  ## According to the kind of metric selected
  ## evaluates if the clusters will be aligned
  th <- sapply(names(thresholds), function(names.th){

    if ((names.th == "Ncor")
        || (names.th =="cor")
        || (names.th =="logocor")
        || (names.th =="Nlocogor")
        || (names.th =="Icor")
        || (names.th =="NIcor")
        || (names.th =="w")
    ) {
      if(hclust.method == "average"){
        mean.scores[names.th] >= thresholds[names.th]
      } else if(hclust.method == "complete"){
        farthest.scores[names.th] >= thresholds[names.th]
      } else if(hclust.method == "single"){
        closest.scores[names.th] >= thresholds[names.th]
      }
    } else if ((names.th == "dEucl")
              || (names.th == "NsEucl")
              || (names.th == "SSD")
              || (names.th == "SW")
              || (names.th == "NSW")
    ){
      if(hclust.method == "average"){
        mean.scores[names.th] <= thresholds[names.th]
      } else if(hclust.method == "complete"){
        farthest.scores[names.th] <= thresholds[names.th]
      } else if(hclust.method == "single"){
        closest.scores[names.th] <= thresholds[names.th]
      }
    }
  })

  ## The condition to align the motifs at the current level is that all
  ## ALL the user-set thresholds must be satisfied.
  if(sum(unlist(th)) == length(thresholds)){
    return(1)
  }else{
    return(0)
  }
}
