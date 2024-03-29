################################################################
## Run hierarchical clustering on the positional profiles of word
## occurrences returned by the RSAT program position-analysis.
##
## Mandatory argument:
##    file.pos=[position_analysis_output_file]
##  The input file (file.pos) must be in the same format as the
##  tab-delimited file returned by position-analysis.
##
## Usage:
## cat ${RSAT}/R-scripts/cluster_position_profiles.R \
##	| R --slave --no-save --no-restore --no-environ \
##          --args "file.pos=position_result.tab"
##
##
## Optional arguments:
##    dir.clusters (default: position_clusters sub-directory in the directory of input file)
##    prefix (default: taken from input file)
##    clust.nb (default: 8)
##    clust.method (default: "complete")
##    sig.threshold (default: 1)
##    rank.threshold (default: 0, meaning "no threshold")
## Example:
##    [...] --args "clust.nb=4; rank.threshold=100; file.pos=position_result.tab"


## Redefine the main directory (this should be adapted to local configuration)
dir.main <- getwd()

dir.rsat <- Sys.getenv("RSAT");
if (dir.rsat == "") {
  stop("The environment variable RSAT is not defined.")
}

## Load some libraries
source(file.path(dir.rsat, 'R-scripts/config.R'))
source(file.path(dir.util, 'util.R'))
source(file.path(dir.util, 'util_plots.R'))
source(file.path(dir.util, 'util_chip_analysis.R'))

################################################################
## Default parameters.  These parameters can be over-written by
## calling the script with command-line arguments.

## Threshold for selecting patterns
sig.threshold <- 1 ## Min level of chi2 significance
rank.threshold <- 0 ## Max number of patterns for the clustering. 0 is interpreted as "no threshold"
clust.nb <- 8 ## Number of clusters
clust.method <- "complete"
clust.suffix <- 'clusters'
heatmap.font.family <- 'mono' ## Font family for the heatmap; Supported: mono (default) | serif | sans

## Drawing preferences
export.plots <- TRUE
display.plots <- FALSE
export.detailed.tables <- FALSE ## Export all the details (normalized frequency profiles, correlation matrices)

display.intervals <- FALSE ## Display intervals rather than class centers
round.positions <- TRUE
all.plots <- FALSE
plot.sep.cluster.profiles <- FALSE ## Set to TRUE to export one image file per cluster (generates many image files)
xlab.by <- 4
shuffle.profiles <- FALSE
#verbosity <- 2

colors <- c("pos"="#00BB00",
            "shuffled"="#BBBBBB")
export.formats.plots <- c("png", "pdf")
heatmap.palette <- blue.to.yellow()


## TEST FILE - For debugging only (JvH)
## file.pos <- "analysis/motifs/position_analysis/SWEMBL_ES_indiff_C3_BN_vs_input_R0.002_summits_sorted_6nt_ci25-1str-noov/SWEMBL_ES_indiff_C3_BN_vs_input_R0.002_summits_sorted_6nt_ci25-1str-noov.tab"
## file.pos <- 'analysis/motifs/position_analysis/SWEMBL_ES_indiff_C3_BN_vs_input_R0.002_summits_sorted_6nt_ci50-1str-noov/SWEMBL_ES_indiff_C3_BN_vs_input_R0.002_summits_sorted_6nt_ci50-1str-noov.tab'
## file.pos <- 'analysis/motifs/position_analysis/SWEMBL_ES_indiff_C3_BN_vs_input_R0.002_summits_sorted_4nt_ci50-1str-noov/SWEMBL_ES_indiff_C3_BN_vs_input_R0.002_summits_sorted_4nt_ci50-1str-noov.tab'
## file.pos <- '/Users/jvanheld/mechali/analysis/motifs/position_analysis/SWEMBL_ES_indiff_C3_BN_vs_F4_RNAse_R0.002_summits_sorted_4nt_ci50-1str-noov/SWEMBL_ES_indiff_C3_BN_vs_F4_RNAse_R0.002_summits_sorted_4nt_ci50-1str-noov.tab'
## file.pos <- '/Users/jvanheld/mechali/analysis/motifs/position_analysis/SWEMBL_ES_indiff_C3_BN_vs_F4_RNAse_R0.002_SICERmatch_summits_4nt_ci50-1str-noov_bg_mkv-2/SWEMBL_ES_indiff_C3_BN_vs_F4_RNAse_R0.002_SICERmatch_summits_4nt_ci50-1str-noov_bg_mkv-2.tab'
## pos.offset <- -25.5
## file.pos <- '/Users/jvanheld/test/positions/results/positions/SWEMBL_mmus_CEBPA_vs_mmus_Input_peaks_R0.05_nof_5nt_ci20-2str-ovlp_top1000.tab'
## file.pos <- '/Users/jvanheld/test/positions/results/positions/SWEMBL_mmus_CEBPA_vs_mmus_Input_peaks_R0.05_nof_5nt_ci20-2str-ovlp_top0.tab'
## file.pos <- '/Users/jvanheld/test/positions/results/positions/SWEMBL_mmus_CEBPA_vs_mmus_Input_peaks_R0.05_nof_skip0_last0_4nt_ci20-2str-ovlp_sig1_mkv-2.tab'
## file.pos <- ''

################################################################
## Read arguments from the command line.
##
## Arguments passed on the command line will over-write the default
## arguments specified above.
args = commandArgs(trailingOnly=TRUE);
if(length(args)==0){
  stop("No arguments supplied. Mandatory: file.pos=[position_analysis_output_file] ")
}else{
#  print("Parsing command-line arguments")
#  print(args)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## Define if the plots have to be generated
draw.plots <- export.plots || display.plots

## Check that input file has been specified
if (!exists("file.pos")) {
  stop("Missing mandatory argument: file.pos=[position_analysis_output_file] ")
}
verbose(paste("Input file", file.pos), 2)

## Define path of directories and files relative to the main directory
dir.pos <- dirname(file.pos)
verbose(paste("Input directory", dir.pos), 2)

## Define prefix for output files
if (!exists("prefix")) {
  prefix <- basename(file.pos)
  prefix <- sub(".tab$", "", basename(file.pos), ignore.case = TRUE, perl = TRUE)
  prefix <- sub(".tsv$", "", basename(file.pos), ignore.case = TRUE, perl = TRUE)
}
verbose(paste("Prefix for output files", prefix), 2)

## Output directory for cluster files (xwe export them to a separate directory because there may be many files)
if (!exists("dir.clusters")) {
  dir.clusters <- file.path(dir.pos, paste(sep='_', prefix, clust.suffix))
}
dir.create(dir.clusters, showWarnings=FALSE, recursive=TRUE)
verbose(paste("Output directory for clusters", dir.clusters), 2)

################################################################
## Read position-analysis result file
verbose(paste("Reading position profiles", file.pos), 2)
pos.data <- read.delim(file.pos, comment.char=";", sep="\t", row.names=NULL)
names(pos.data)[1] <- "seq"
pos.data$seq <- toupper(pos.data$seq)
pos.data$id <- toupper(pos.data$id)

## Extract k-mer sequences
kmer.seq <- as.vector(pos.data[,1])

## define k-mer length
kmer.lengths <- nchar(kmer.seq)
if (max(kmer.lengths) == min(kmer.lengths)) {
  kmer.len <- max(kmer.lengths)
} else {
  kmer.len <- "k"
}


## Extract k-mer description (may include reverse complements)
kmer.desc <- as.vector(pos.data[,2])
row.names(pos.data) <- kmer.desc

## Define the columns containing occurrences per position window
if (!exists("column.offset")) {
  ## Detect columns not preceding the distribution (those returning other types of statistics)
  column.offset <- length(
      intersect(toupper(colnames(pos.data)),
                toupper(c("X.seq","seq", "id","occ","over","chi2","df","Pval","Eval", "sig","rank"))))
}
profile.col <- (column.offset+1):ncol(pos.data) ## columns containing the position profiles (occurrences per window)

# message("colnames(pos.data)\t", paste(collapse=", ", colnames(pos.data)))
# message("column.offset\t", column.offset)

## Fix the headers of window occurrence columns
colnames(pos.data)[1] <- sub('X.', '', colnames(pos.data)[1])
colnames(pos.data)[profile.col] <- sub("X\\.", "-", colnames(pos.data)[profile.col], perl="TRUE")
colnames(pos.data)[profile.col] <- sub("X", "+", colnames(pos.data)[profile.col], perl="TRUE")

## Fix a bug with previous version of position-analysis (exported a column of 0s)
last.profile.col <- profile.col[length(profile.col)]
last.profile.col.header <- as.numeric(colnames(pos.data)[last.profile.col])
if (is.na(last.profile.col.header)) {
  profile.col <- profile.col[1:(length(profile.col)-1)]
}

nb.windows <- length(profile.col)
# message('nb.windows = ', nb.windows)


## Define the selected patterns
selected.patterns <- (pos.data$sig >= sig.threshold) 

## Apply rank threshold if specified
if (rank.threshold > 0) {
  selected.patterns <- selected.patterns & (pos.data$rank <= rank.threshold)
}
nb.patterns <- sum(selected.patterns)
verbose(paste(nb.patterns , "selected patterns"), 2)

if (nb.patterns < 2) {
  stop("Clustering is irrelevant with less than two selected patterns")
}

## ##############################################################
## Create a separate data frame for the position profiles.
## This is redundant (-> memory-inefficient) but convenient.
pos.profiles <- pos.data[selected.patterns, profile.col]
ref.positions <- as.numeric(colnames(pos.profiles))

# message("ref.positions\t", paste(collapse=", ", ref.positions))

if (exists("pos.offset")) {
  ref.positions <- ref.positions + pos.offset
}

## Round positions to avoid .5
if (round.positions) {
  ref.positions <- round(ref.positions)
}

if (display.intervals) {
  class.mid <- ref.positions
  class.interval <- class.mid[2] - class.mid[1]
  class.min <- class.mid - (class.interval-1)/2
  class.max <- class.mid + (class.interval-1)/2
  colnames(pos.profiles) <- paste(sep='', '[', class.min, '-', class.max, '[')
} else {
  colnames(pos.profiles) <- ref.positions
}

## ##############################################################
## Normalize position profiles by row (convert k-mer occurrences in
## relative frequencies)
sum.per.kmer <- apply(pos.profiles, 1, sum)
pos.profiles.freq <- pos.profiles / sum.per.kmer
pos.profiles.freq.norm <- pos.profiles.freq - apply(pos.profiles.freq, 1, median)
if (export.detailed.tables) {
  export.object(pos.profiles.freq.norm, file=file.path(dir.clusters, paste(sep='_', prefix, 'position_profiles_norm_freq')), export.format='table')
}

## ##############################################################
## Perform a random shuffling of the position profiles, as a negative
## control for correlation and clustering.
if (shuffle.profiles) {
  verbose("Shuffling position profiles for negative controls", 2)
  shuffled.values <- sample(unlist(pos.profiles), replace=FALSE)
  # message("dim(pos.profiles)", paste(collapse = ", ", dim(pos.profiles)))
  shuffled.profiles <- as.data.frame(matrix(nrow=nrow(pos.profiles), ncol=ncol(pos.profiles), shuffled.values))
  rownames(shuffled.profiles) <- rownames(pos.profiles)
  colnames(shuffled.profiles) <- colnames(pos.profiles)
  shuffled.profiles.freq <- shuffled.profiles / sum.per.kmer
  if (export.detailed.tables) {
    export.object(shuffled.profiles.freq, file=file.path(dir.clusters, paste(sep='_', prefix, 'shuffled_profiles_norm_freq')), export.format='table')
  }
}

## ##############################################################
## Correlation analysis

## Compute correlations between profiles
verbose("Computing correlation between profiles", 2)
pos.cor <- cor(t(pos.profiles))
if (shuffle.profiles) {
  shuffled.cor <- cor(t(shuffled.profiles))
}
if (export.detailed.tables) {
  export.object(round(pos.cor, digits=3), file=file.path(dir.clusters, paste(sep='_', prefix, 'profiles_correlations')), export.format='table')
  if (shuffle.profiles) {
    export.object(round(shuffled.cor, digits=3), file=file.path(dir.clusters, paste(sep='_', prefix, 'shuffled_correlations')), export.format='table')
  }
}

## ##############################################################
## Hierarchical clustering
verbose(paste(sep="", "Hierarchical clustering; method=", clust.method), 2)
pos.tree <- hclust(as.dist(1-pos.cor), method=clust.method)
if (shuffle.profiles) {
  shuffled.tree <- hclust(as.dist(1-shuffled.cor), method=clust.method)
}

## Cut the tree in k clusters
clust.nb <- min(clust.nb, nb.patterns) ## Check that number of clusters does not exceed the number of selected patterns
verbose(paste(sep="", "Cutting the tree; k=", clust.nb), 2)
clusters <- cutree(pos.tree,k=clust.nb)
clusters <- clusters[order(clusters)]
cluster.table <- data.frame('seq'=pos.data[names(clusters),'seq'],
                            'cluster'=clusters,
                            'chi2'=pos.data[names(clusters),'chi2'],
                            'sig'=pos.data[names(clusters),'sig'],
                            'Eval'=pos.data[names(clusters),'Eval'],
                            'identifier'=names(clusters),
                            row.names=NULL)

## Export the cluster file in the same directory as the position
## profiles. All other files are exported in the separate cluster
## directory, to avoid confusion between the numerous output file.
cluster.file <- file.path(dir.pos, paste(sep='', prefix, '_',clust.suffix,'.tab'))
write.table(cluster.table, file=cluster.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

## Report the number of elements per cluster
verbose(paste("Cluster file:", cluster.file), 2)
cluster.sizes <- table(clusters)
cluster.sizes <- as.data.frame(cluster.sizes)
names(cluster.sizes) <- c("cluster", "n")
write.table(cluster.sizes, file=file.path(dir.clusters, paste(sep='', prefix, '_cluster_sizes.tab')), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

################################################################
## Compute median profile for each cluster
median.profiles <- data.frame(matrix(nrow=clust.nb, ncol=ncol(pos.profiles.freq.norm)))
colnames(median.profiles) <- colnames(pos.profiles.freq.norm)
rownames(median.profiles) <- paste("cl", 1:clust.nb, sep="")
for (i in 1:clust.nb) {
  cluster.size <- sum(clusters==i)
  cluster.kmers <- names(clusters[clusters==i])
  median.profiles[i,] <- apply(pos.profiles.freq.norm[cluster.kmers,], 2, median)
}
export.object(median.profiles, file=file.path(dir.clusters, paste(sep='', prefix, '_median_profiles_per_cluster')), export.format='table')

################################################################
## Draw plots if required
if (display.plots) {
  plot.device.format <- "x11"
} else {
  plot.device.format <- "pdf"
}

if (draw.plots) {

  ## ##############################################################
  ## Principal plots  
  ## Draw a heatmap of the clustered position profiles  
  file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_profile_heatmap'))
                                        #  X11(width=12, height=12)
  open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=12,height=12)  
  zmax <- max(abs(range(pos.profiles.freq.norm)))*1.01
  par(family=heatmap.font.family)
  par(cex=0.1)
  par(cex.axis=0.1)
  par(cex.lab=0.1)
  par(cex.main=1)
  ymargin=min(8,max(kmer.lengths)*2)
  xmargin=5
  heatmap(as.matrix(pos.profiles.freq.norm),
          scale="none",
          main=paste("Clustered position profiles; ", 'd=1-cor; ', clust.method, "linkage"),
          zlim=c(-zmax,zmax),
          col=green.white.red(),
          xlab='position', 
          ylab=paste(kmer.len, '-mers', sep=''), 
          Colv=NA, Rowv=as.dendrogram(pos.tree),
          mar=c(xmargin,ymargin)
          )
  if (plot.device.format == "x11") {
    export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=12,height=12)
  } else {
    dev.off()
  }

  ## ##############################################################
  ## Plot frequency profiles for each cluster
  file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_freq_profiles_per_cluster'))
  if (sqrt(clust.nb) == round(sqrt(clust.nb))) {
    mfrow.rows <- sqrt(clust.nb)
    mfrow.cols <- sqrt(clust.nb)
  } else {
    mfrow.rows <- min(clust.nb, 4)
    mfrow.cols <- max(1, ceiling(clust.nb/mfrow.rows))
  }

  ## x11(width=6*mfrow.cols, height=4*mfrow.rows)
  open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=6*mfrow.cols, height=4*mfrow.rows)
  par(mfrow=c(mfrow.rows, mfrow.cols))
  max.freq <- max(as.vector(as.matrix(pos.profiles.freq.norm)))
  min.freq <- min(as.vector(as.matrix(pos.profiles.freq.norm)))
  for (i in 1:clust.nb) {
    cluster.size <- sum(clusters==i)
    profile.stats <- plot.profiles(as.data.frame(pos.profiles.freq.norm[names(clusters[clusters==i]),]),
                                        #                ylim=c(min(pos.profiles.freq.norm[clusters==i,]),max(pos.profiles.freq.norm[clusters==i,])),
                                   main=paste(sep='', kmer.len, '-mer; method=', clust.method,'; clust ', i, '/', clust.nb),
                                        #                                   main=paste(sep='', 'cluster ', i, '/', clust.nb),
                                   col.profiles="gray",
#                                   col.profiles=rainbow(n=cluster.size),
                                   plot.median.profile=T, plot.mean.profile=F, plot.sd.profile=F, xlab.by=xlab.by,las=2,
                                   ylim=c(min.freq, max.freq),
                                   ylab='Longitudinal frequencies'
                                   )
  }
  par(mfrow=c(1,1))
  if (plot.device.format == "x11") {
    export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=6*mfrow.cols, height=4*mfrow.rows)
  } else {
    dev.off()
  }

  ################################################################
  ## A single plot with the median profile for each cluster
  file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_median_profiles_per_cluster'))
  ##  x11(width=12,height=7)
  open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=12,height=7)
  median.profile.stats <- plot.profiles(median.profiles,
                                        main=paste(sep='', "Median profiles per ",kmer.len,"-mer cluster"),
                                        #                                        col.profiles="gray",
                                        col.profiles=rainbow(n=clust.nb+2),
                                        plot.median.profile=F, plot.mean.profile=F, plot.sd.profile=F,
                                        #                                        lty=c("solid", "dashed"),
                                        xlab='Position', lwd=2, xlab.by=xlab.by,las=2,
                                        ylab='Longitudinal frequencies'
                                        )
  if (plot.device.format == "x11") {
    export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=12, height=7)
  } else {
    dev.off()
  }
  

  ## Plot occurrence profiles for all clusters together
  file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_occ_profiles_per_cluster'))
  ##  x11(width=6*mfrow.cols, height=4*mfrow.rows)
  open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=6*mfrow.cols, height=4*mfrow.rows)
  par(mfrow=c(mfrow.rows, mfrow.cols))
  max.occ <- max(as.vector(as.matrix(pos.profiles)))
  min.occ <- min(as.vector(as.matrix(pos.profiles)))
  for (i in 1:clust.nb) {
    cluster.size <- sum(clusters==i)
    plot.profiles(as.data.frame(pos.profiles[names(clusters[clusters==i]),]),
                                        #                ylim=c(0,max(pos.profiles[clusters==i,])),
                  main=paste(sep='', kmer.len, '-mer cluster ', i, '/', clust.nb),
                                   col.profiles="gray",
#                  col.profiles=rainbow(n=cluster.size),
                  plot.median.profile=T, plot.mean.profile=F, plot.sd.profile=F, xlab.by=xlab.by,las=2,
                  ylim=c(min.occ, max.occ),
                  ylab='Occurrences'
                  )
  }
  if (plot.device.format == "x11") {
    export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=6*mfrow.cols, height=4*mfrow.rows)
  } else {
    dev.off()
  }


  ## ##############################################################
  ## Optional plots
  if (all.plots) {

    ## Draw distribution of chi2 significance scores
    file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_sig_distrib'))
    open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=7,height=5)
    hist(pos.data$sig, breaks=100, main=paste("chi2 sig distribution"), xlab="sig of chi2 test", ylab=paste(sep='', "Number of ",kmer.len,"-mers"), col="#DDDDDD")    
    if (display.plots) {
      export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=7,height=5)
    } else {
      dev.off()
    }
    
    ## Draw the distributions of correlation values (position profiles + shuffled data)
    file.prefix <- file.path(dir.clusters, paste(sep='', '_profile_correlations'))
    if (shuffle.profiles) {
      open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=7,height=8)
      par(mfrow=c(2,1))
    } else {
      open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=7,height=5)
      par(mfrow=c(1,1))
    }
    hist(pos.cor, breaks=100,
         xlim=c(-1,1), col=colors["pos"],
         main=paste(sep='', "correlations between ",kmer.len,"-mer profiles"),
         xlab="Correlation",
         ylab=paste(sep='', "Pairs of ",kmer.len,"-mers"))
    if (shuffle.profiles) {
        hist(shuffled.cor, breaks=100, xlim=c(-1,1), col=colors["shuffled"],
         main=paste("correlations between shuffled profiles"),
         xlab="Correlation",
         ylab=paste(sep="", "Pairs of ",kmer.len,"-mers"))
    }
    par(mfrow=c(1,1))
    if (display.plots) {
      if (shuffle.profiles) {
        export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=7, height=8)
      } else {
        export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=7, height=5)
      }
    } else {
      dev.off()
    }

    ## Draw a heatmap of the correlation matrix
    file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_profile_correlations_heatmap'))
    open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=12,height=12)
    heatmap(pos.cor, scale="none",
            main='Correlations between position profiles',
            zlim=c(-1,1),
            col=heatmap.palette)
    if (display.plots) {
      export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=12,height=12)
    } else {
      dev.off()
    }
    
    ## Draw a heatmap of the correlation matrix of shuffled profiles
    if (shuffle.profiles) {
      file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_shuffled_correlations_heatmap'))
      open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=12,height=12)
      heatmap(shuffled.cor, scale="none",
               main='Correlations between shuffled profiles',
              zlim=c(-1,1),
              col=heatmap.palette)
      export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=12,height=12)
       if (display.plots) {
        export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=12,height=12)
      } else {
        dev.off()
      }
    }    

    ## Plot the tree structure for position profiles and shuffled data
    file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_profile_tree')) 
    if (shuffle.profiles) {
      open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=16,height=8)
      par(mfrow=c(2,1))
    } else {
      open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=16,height=5)
      par(mfrow=c(1,1))
    }
    plot(pos.tree, main=paste(prefix, ": clustering of position profiles;", "dist= 1 - cor; ", clust.method, "linkage"), xlab=NA)
    if (shuffle.profiles) {
      plot(shuffled.tree, main=paste(prefix, ": clustering of shuffled profiles;", "dist= 1 - cor; ", clust.method, "linkage"), xlab=NA)
    }
    par(mfrow=c(1,1))
    if (display.plots) {
      if (shuffle.profiles) {
        export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=16,height=8)
      } else {
        export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=16,height=5)
      }
    } else {
      dev.off()
    }
  
    ## ##############################################################
    ## If requested, print one separate figure of word occurrence and
    ## frequency profiles for each cluster
    if (plot.sep.cluster.profiles) {

      ## Plot normalized frequency profiles for each cluster
      for (i in 1:clust.nb) {
        file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_freq_profile_clust', i, '_'))
        open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=10,height=5)
        par(cex=0.7)
        par(cex.lab=0.7)
        cluster.size <- sum(clusters==i)
        plot.profiles(as.data.frame(pos.profiles.freq.norm[names(clusters[clusters==i]),]),
                                        #                ylim=c(0,max(pos.profiles.freq.norm[clusters==i,])),
                      main=paste(sep='', prefix, '; cluster ', i, '/', clust.nb),
                                   col.profiles="gray",
#                      col.profiles=rainbow(n=cluster.size),
                      plot.median.profile=T, plot.mean.profile=F, plot.sd.profile=F, xlab.by=xlab.by,las=2,
                      ylab='Longitudinal frequencies'
                      )
        if (display.plots) {
          export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=10, height=5)
        } else {
          dev.off()
        }
      }
      
      ## Plot occurrence profiles for each cluster
      for (i in 1:length(table(clusters))) {
        file.prefix <- file.path(dir.clusters, paste(sep='', prefix, '_occ_profile_clust', i, '_'))
        open.plot.device(file.prefix=file.prefix, format=plot.device.format, width=10,height=5)
        par(cex=0.7)
        par(cex.lab=0.7)
        cluster.size <- sum(clusters==i)
        plot.profiles(as.data.frame(pos.profiles[names(clusters[clusters==i]),]),
                                        #                ylim=c(0,max(pos.profiles[clusters==i,])),
                      main=paste(sep='', prefix, '; cluster ', i, '/', clust.nb),
                                   col.profiles="gray",
#                      col.profiles=rainbow(n=cluster.size),
                      plot.median.profile=T, plot.mean.profile=F, plot.sd.profile=F, xlab.by=xlab.by,las=2,
                      ylab='Longitudinal frequencies'
                      )
        if (display.plots) {
          export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=10, height=5)
        } else {
          dev.off()
        }
      }
    }
  }
}


verbose("job done", 2)
verbose(dir.clusters, 2)
