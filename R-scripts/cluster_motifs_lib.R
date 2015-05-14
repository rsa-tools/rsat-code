###################################################################
## Library with R functions required in R-scripts/cluster_motifs.R

################################################################
## Check parameter values
check.param <- function() {
  ## Check that input file has been specified
  if (!exists("infile")) {
    stop("Missing mandatory argument: infile=[matrix_comparison_table] ")
  }
  verbose(paste("Input file", infile), 3)

  ## Check that description file
  if (!exists("description.file")) {
    stop("Missing mandatory argument: description.file=[matrix_description_table] ")
  }
  verbose(paste("Description file", description.file), 3)

  ## Check that output file has been specified
  if (!exists("out.prefix")) {
    stop("Missing mandatory argument: out.prefix=[output_prefix] ")
  }
  verbose(paste("Output prefix", out.prefix), 3)

  ## Check that distance table file has been specified
  if (!exists("distance.table")) {
    distance.table <<- paste(sep="", out.prefix, "_dist_table.tab")
  }
  verbose(paste("Distance table", distance.table), 3)

  ## Default score is the normalized correlation
  if (!exists("metric")) {
    metric <<- "Ncor";
  }

  ## Default hclust method is the complete method
  if (!exists("hclust.method")) {
    hclust.method <<- "average";
  }

  ## Option indicating if the heatmap must be computed
  if (!exists("draw.heatmap")) {
    draw.heatmap <<- "1";
  }

  ## Option indicating if the tree
  ## in Newick format should be exported
  if (!exists("export.newick")) {
    export.newick <<- "0";
  }

  ## Option indicating if the tree with aligned consensuses must be computed
  if (!exists("draw.consensus")) {
    draw.consensus <<- "1";
  }

  ## Option indicating if the heatmap must be computed
  if (!exists("pos.hclust.in.heatmap")) {
    pos.hclust.in.heatmap <<- "column";
  }

  ## Check that when the option -heatmap_position_tree is indicated in the
  ## perl script, the option -heatmap should also indicated, conversely
  ## the program ends
  if (exists("pos.hclust.in.heatmap") & !exists("draw.heatmap")) {
    stop("You must indicate the option -heatmap in order to select the option -heatmap_position_tree.")
  }

  ## When this option is activated
  ## Only the hclust is done, skipping the
  ## trees and heatmap computation
  if (!exists("only.hclust")) {
    only.hclust <<- "0";
  }

  ## Define the kind of metric used: scores or distances
  supported.scores <- c("cor", "Ncor")
  supported.distances <- c("dEucl", "NdEucl")

  if(metric %in% supported.scores){

    ## Default lower and upper thresholds equals to zero
    if (!exists("lth")) {
      lth <<- list()
      lth[["Ncor"]] <<- 0;
      lth[["cor"]] <<- 0;
      lth[["w"]] <<- 0;

      lth.values <<- unlist(lth)
      lth.scores <<- names(lth.values)
    }


  } else if(score %in% supported.distances){

    if (!exists("uth")) {
      uth <<- list()
      uth[["dEucl"]] <<- 1;
      uth[["NdEucl"]] <<- 1;
      uth[["w"]] <<- 0;

      uth.values <<- unlist(uth)
      uth.scores <<- names(uth.values)
    }
  }


  ###############################
  ## Read the lower thresholds
  if(exists("lthsp")){
    lthsp <- unlist(strsplit(lthsp, "_"))
    lth.scores <<- lthsp[seq(1,length(lthsp), by = 2)]
    lth.values <<- as.numeric(lthsp[seq(2,length(lthsp), by = 2)])

    for(i in 1:length(lth.scores)){
      thresholds[[lth.scores[i]]] <<- lth.values[i]
    }

    supported <- c("cor", "Ncor", "w")
    if(length(setdiff(supported, lth.scores)) > 0){
      for(add in setdiff(supported, lth.scores)){
        lth.scores <<- append(lth.scores, add)
        lth.values <<- append(lth.values, 0)
      }
    }
  }


  ###############################
  ## Read the upper thresholds
  if(exists("uthsp")){
    uthsp <- unlist(strsplit(uthsp, "_"))
    uth.scores <<- uthsp[seq(1,length(uthsp), by = 2)]
    uth.values <<- as.numeric(uthsp[seq(2,length(uthsp), by = 2)])

    for(i in 1:length(uth.scores)){
      thresholds[[uth.scores[i]]] <<- uth.values[i]
    }

    supported <- c("dEucl", "NdEucl")
    if(length(setdiff(supported, uth.scores)) > 0){
      for(add in setdiff(supported, uth.scores)){
        uth.scores <<- append(uth.scores, add)
        uth.values <<- append(uth.values, 0)
      }
    }
  }
}


#####################################################
## Creates the folders where the branch-motifs
## for each merge level of the tree will be stored
create.dir.merge <- function(level){

  ## Create the folder with the merged consensuses
  merge.dir <- paste("merged_consensuses/", level, sep = "")
  dir.create(file.path(cluster.folder, merge.dir), showWarnings = FALSE, recursive = TRUE)
  new.dir <- file.path(cluster.folder, merge.dir)

  flag <- system(paste("ls ", new.dir, "/ | wc -l", sep = ""), intern = TRUE)
  if(flag >= 1){
    system(paste("rm -r ", new.dir, "/*", sep = ""))
  }
}
