#####################################################################
## Evaluate if the mean of scores for all pairs of leaves between
## the two clusters is upper/lower than the threshold,
## if so, then the two clusters should be aligned

alignment.test <- function(id1, id2, thresholds = list(Ncor = 0.4, cor = 0.6, w = 5), hclust.method = "average", hclust.metric = "Ncor"){

  ## Calculate the mean of the scores for all the pairs of motifs
  if(hclust.method == "average" | hclust.method == "centroid"){
    compa.numbers <- get.comparison.number(id1, id2)

    scores <- global.compare.matrices.table[compa.numbers, names(thresholds)]

    ## Calculate the median of the data
    mean.scores <- apply(scores, 2, mean)

  ## Calculate the median between all the pairs of motifs
  } else if(hclust.method == "median"){

    scores <- global.compare.matrices.table[compa.numbers, names(thresholds)]

    median.scores <- apply(scores, 2, median)

    ## Calculate the closest motifs between all the pairs of motifs
  } else if(hclust.method == "complete"){

    farthest.motifs <- closest.or.farthest.motifs.ids(id1, id2, metric = hclust.metric, closest = FALSE)
    id1.far <- farthest.motifs[1]
    id2.far <- farthest.motifs[2]

    compa.numbers <- get.comparison.number(id1.far, id2.far)

    ## Get the scores of the comparisons
    farthest.scores <- global.compare.matrices.table[compa.numbers, names(thresholds)]

  ## Calculate the closest motifs between all the pairs of motifs
  } else if(hclust.method == "single"){

    closest.motifs <- closest.or.farthest.motifs.ids(id1, id2, metric = hclust.metric, closest = TRUE)
    id1.close <- closest.motifs[1]
    id2.close <- closest.motifs[2]

    compa.numbers <- get.comparison.number(id1.close, id2.close)

    ## Get the scores of the comparisons
    closest.scores <- global.compare.matrices.table[compa.numbers, names(thresholds)]
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
        || (names.th == "logoDP")
    ) {
      if(hclust.method == "average" | hclust.method == "centroid"){
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
      if(hclust.method == "average" | hclust.method == "centroid"){
        mean.scores[names.th] <= thresholds[names.th]
      } else if(hclust.method == "complete"){
        farthest.scores[names.th] <= thresholds[names.th]
      } else if(hclust.method == "median"){
        median.scores[names.th] <= thresholds[names.th]
      }
      else if(hclust.method == "single"){
        closest.scores[names.th] <= thresholds[names.th]
      }
    }
  })


  if(sum(unlist(th)) == length(thresholds)){
    return(1)
  }else{
    return(0)
  }
}
