#################################
## Build a distance matrix from the distance score list
build.distance.matrix <- function(comparison.table,
                                  score="Ncor"){

  dist.table <- NULL
  distances.objects <- list()

  ## Extract score values
  score.values <- comparison.table[,score]

  ## If required, convert similarities to distances
  if (score == "Ncor") {
    score.dist <- 1 - score.values
  } else {
    stop('Jaime still needs to call convert.scores() from here')
  }

  ## Add a column with score column to the compare matrices table, will
  ## be used to generate a cross-table
  comparison.table$score.dist <- score.dist

  ## Build the distance table from the column score
  dist.table <- xtabs(score.dist ~ id1+id2, comparison.table)


  #  dist.table <- matrix(unlist(dist.table),ncol=ncol(dist.table))

  ## Ensure that symmetrical distances are defined
  dist.table.sym <- pmax(dist.table,t(dist.table))

  ## Cast the distance table into an object of class "dist"
  dist.matrix <- as.dist(dist.table)

  return(list(dist.table, dist.matrix))
}

