###########################################################
## Build a distance matrix from the distance score list
build.distance.matrix <- function(comparison.table,
                                  score="Ncor"){


  dist.table <- NULL
  distances.objects <- list()

  ## Extract score values
  score.values <- comparison.table[,score]

  ## If required, convert similarities to distances
  ## Similarity sores bounded to 1
  if ((score == "Ncor")
      || (score=="cor")
      || (score=="logocor")
      || (score=="Nlocogor")
      || (score=="Icor")
      || (score=="NIcor")
  ) {
    ## cor       Pearson correlation (computed on residue occurrences in aligned columns)
    ## Ncor 			Relative width-normalized Pearson correlation
    ## logocor 			correlation computed on sequence logos
    ## Nlogocor 			Relative width-normalized logocor
    ## Icor 			Pearson correlation computed on Information content
    ## NIcor 			Relative width-normalized Icor
    score.dist <- 1 - score.values

  } else if ((score == "logoDP")
             || (score == "cov")) {
    ## logoDP 			dot product of sequence logos
    ## cov 			covariance between residues in aligned columns

    stop("logoDP and cov scores are not supported yet")

  } else if ((score == "dEucl")
             || (score == "NdEucl")
             || (score == "SSD")
             || (score == "SW")
             || (score == "NSW")
  ) {
    ## dEucl 			Euclidian distance between residue occurrences in aligned columns
    ## NdEucl 			Relative width-normalized dEucl
    ## NsEucl 			similarity derived from Relative width-normalized Euclidian distance
    ## SSD 			Sum of square deviations
    ## SW 			Sandelin-Wasserman
    ## NSW 			Relative width-normalized Sandelin-Wasserman

    score.dist <- score.values

  } else {
    stop(paste(score, "is not a valid score", sep="\t"))
  }


  ## Add a column with score column to the compare matrices table, will
  ## be used to generate a cross-table
  comparison.table$score.dist <- score.dist

  ## Build the distance table from the column score
  dist.table <- xtabs(score.dist ~ id1+id2, comparison.table)


  ## Ensure that symmetrical distances are defined
  dist.table.sym <- pmax(dist.table,t(dist.table))

  ## Cast the distance table into an object of class "dist"
  dist.matrix <- as.dist(dist.table.sym)

  return(list(table = dist.table.sym, matrix = dist.matrix))
}

