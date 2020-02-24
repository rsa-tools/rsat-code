#####################################################
## Convert the hclust object to a character object
## with the lines ready to print a JSON file
convert.hclust.to.JSON <- function(tree){

  ## Require ctc if it is required
  if(!require("RJSONIO")){
    install.packages("RJSONIO")
  }

  ## Load library
  suppressPackageStartupMessages(library("RJSONIO", warn.conflicts=FALSE))

  #######################################################
  ## Extract the tree from an hclust object
  createLeafNode <- function(hclust, i) {
    list(label = hclust$labels[[i]],
         order = hclust$order[[i]])
  }

  ########################################################
  ## Convert an hclust tree into a JSON format tree
  hclustToTree <- function(hclust) {
    if (length(hclust$merge) == 0)
      return(NULL)
    merges <- list()
    for (index in 1:nrow(hclust$merge)) {
      left <- hclust$merge[index, 1]
      right <- hclust$merge[index, 2]
      if (left < 0)
        left <- createLeafNode(hclust, -left)
      else
        left <- merges[[left]]
      if (right < 0)
        right <- createLeafNode(hclust, -right)
      else
        right <- merges[[right]]
      if (left$order > right$order) {
        tmp <- left
        left <- right
        right <- tmp
      }
      merges[[index]] <- list(
        children = list(
          left,
          right
        ),
        order = left$order
      )
    }
    return(merges[nrow(hclust$merge)])
  }

  ### Creates and parse the json string
  halfway.tree <- hclustToTree(tree)
  jsonTree <- toJSON(halfway.tree)


  ## Fix some little technical issues for JSON compatibility with the tree display javascript
  jsonTree <- gsub("\\],", "\\]", jsonTree, perl = TRUE)
  jsonTree <- paste("{\n\"name\": \"\",\n\"children\":", jsonTree, "}", sep = "")
  jsonTree <- gsub("\n\"order\":\\s+\\d+", "", jsonTree, perl = TRUE)

  return(jsonTree)
}