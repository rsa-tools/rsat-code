########################################################
## As the JSON tree order is not the same as the hclust tree
## the JSON tree must be explored in order to rename the branches
## and thus can assign a name which will be used to add the
## branch consensuses (in the perl code)
identify.JSON.tree.branches <- function(tree){

  ## Change the labels for the numbers in the descripton table
  tree$labels <<- NULL
  tree$labels <<- as.vector(global.description.table$n)

  ## Erase the (',')
  jsonTree <- convert.hclust.to.JSON(tree)
  jsonTree <- gsub("\\],", "\\]", jsonTree, perl = TRUE)
  jsonTree <- gsub("\n\"order\":\\s+\\d+", "", jsonTree, perl = TRUE)

  ## Save just the numbers on each branch of the tree
  copy.Json <- jsonTree
  copy.Json <- gsub("[^ \\d+ \\[ \\]  \\.]", "", copy.Json, perl = TRUE)
  copy.Json <- gsub("(\\d+)", "|\\1", copy.Json, perl = TRUE)

  ## Convert the string in a vector, save those positions with
  ## either a number, "[" or "]"
  copy.Json.vector <- unlist(strsplit(copy.Json, " "))
  copy.Json.vector <- copy.Json.vector[2:length(copy.Json.vector)]
  copy.Json.vector <- copy.Json.vector[which(copy.Json.vector != "")]

  ## Saves the positions with a "[", except the first one,
  ## they will be used to find the subgroups within the tree
  symbol.counter <- length(which(copy.Json.vector == "["))
  symbol <- which(copy.Json.vector == "[")[2:symbol.counter]


  ## A dataframe will be created using this variables (col1 and col2)
  ## They are created outside the loops and are treated as global variables
  col1 <- NULL
  col2 <- NULL
  up.pos <- 1
  sapply(1:(symbol.counter-1), function(j){

    ## Counter of condition (i == j)
    cond.counter <<- 0
    col1 <<- append(col1, j)
    up.count <- 0
    down.count <- 0
    sapply(symbol[j]:length(copy.Json.vector) ,function(i){

      ## Here the code search for the [] in order to identify the
      ## subgroups and then can re-name them similarly to the hclust object
      if(copy.Json.vector[i] == "["){
        up.count <<- up.count + 1

        ## Positions up
        if(up.count == 1){
          up.pos <<- i
        }
      } else if(copy.Json.vector[i] == "]"){
        down.count <<- down.count + 1
      }

      ## Condition
      if(up.count == down.count){

        cond.counter <<- cond.counter + 1

        if(cond.counter == 1){
          ## Position down
          down.pos <- i
          temp <- NULL
          temp <- gsub("[^\\d+ \\|]", "" ,copy.Json.vector[up.pos:down.pos], perl = TRUE)
          numb <- NULL
          numb <- as.integer(unlist((strsplit(paste(temp, collapse = ""), "\\|"))))
          numb <- numb[which(numb != "NA")]
          col2 <<- append(col2, paste(get.id(numb), collapse = ","))
        }
      }
    })
  })

  ## Create the dataframe and rename the columns
  JSON.clusters.table <- data.frame(col1)
  JSON.clusters.table$cluster <- col2
  colnames(JSON.clusters.table) <- c(";level", "cluster")

  cluster <- NULL
  sapply(1:nrow(JSON.clusters.table), function(x){

    ## Save the Ids corresponding to each branch of the JSON tree
    leaves.JSON <- unlist(strsplit(JSON.clusters.table[x,2], ","))

    sapply(1:nrow(tree$merge), function(y){

      cond.counter <<- 0
      ## Save the Ids corresponding to each branch of the hclust tree
      leaves.merge <- NULL
      leaves.merge <- sort(leaves.per.node(tree)[[y]])

      if(length(leaves.JSON) == length(leaves.merge)){

        if(sort(get.id(leaves.merge))[1] == sort(leaves.JSON)[1]){

          cond.counter <<- cond.counter + 1

          if(cond.counter == 1){
            cluster <<- append(cluster, paste("node_", y, sep = ""))
          }
        }
      }
    })
  })

  JSON.clusters.table$node <- cluster
  return(JSON.clusters.table)
}
