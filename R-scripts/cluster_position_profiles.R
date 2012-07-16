################################################################
## Run hierarchical clustering on the positional profiles of word
## occurrences returned by the RSAT program position-analysis.
##
## Mandatory argument:
##    file.pos=[position_analysis_output_file]
##  The input file (file.pos) must be in the same format as the
##  tab-delimited file returned by position-analysis.
##
## Optional arguments:
##    dir.clusters (default: position_clusters sub-directory in the directory of input file)
##    suffix (default: taken from input file)
##    nb.clusters (default: 5)
##    sig.threshold (default: 1)
##    rank.threshold (default: 100)

## Load some libraries
source('http://www.bigre.ulb.ac.be/Users/jvanheld/statistics_bioinformatics/R-files/config.R')
source(file.path(dir.util, 'util_chip_analysis.R'))
source(file.path(dir.util, 'microarray_util.R'))

## Redefine the main directory (this should be adapted to local configuration)
dir.main <- getwd()

################################################################
## Default parameters.  These parameters can be over-written by
## calling the script with command-line arguments.

## Threshold for selecting patterns
sig.threshold <- 1 ## Min level of chi2 significance
rank.threshold <- 100 ## Max number of patterns for the clustering
nb.clusters <- 5 ## Number of clusters

## Drawing preferences
colors <- c("pos"="#00BB00",
            "shuffled"="#BBBBBB")
export.formats.plots <- c("png", "pdf")
heatmap.palette <- blue.to.yellow()

################################################################
## Read arguments from the command line.
##
## Arguments passed on the command line will over-write the default
## arguments specified above.
args = commandArgs(trailingOnly=TRUE);

if(length(args)==0){
  stop("No arguments supplied. Mandatory: file.pos=[position_analysis_output_file] ")
}else{
  print("Parsing command-line arguments")
  print(args)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## Check that input file has been specified
if (!exists("file.pos")) {
  stop("Missing mandatory argument: file.pos=[position_analysis_output_file] ")
}
verbose(paste("Input file", file.pos), 1)

## Define path of directories and files relative to the main directory
dir.pos <- dirname(file.pos)
verbose(paste("Input directory", dir.pos), 1)

## Output directory
if (!exists("dir.clusters")) {
  dir.clusters <- file.path(dir.pos, "position_clusters")
}
dir.create(dir.clusters, showWarnings=FALSE, recursive=TRUE)
verbose(paste("Output directory for clusters", dir.clusters), 1)


## Define suffix for output files
if (!exists("suffix")) {
  suffix <- basename(file.pos)
  suffix <- sub(".tab$", "", basename(file.pos), ignore.case = TRUE, perl = TRUE)
}
verbose(paste("Suffix for output files", suffix), 1)

################################################################
## Read position-analysis result file
pos.data <- read.delim(file.pos, comment.char=";", sep="\t", row.names=1)

## TEMPORARY: fix a bug with position-analysis output file (header
## contains one tab at the end -> considered by R as a column of 0
pos.data <- pos.data[, 1:(ncol(pos.data)-1)]

profile.col <- 10:ncol(pos.data) ## columns containing the position profiles (occurrences per window)
nb.windows <- length(profile.col)
colnames(pos.data)[1] <- sub('X.', '', colnames(pos.data)[1])
colnames(pos.data)[profile.col] <- sub("X\\.", "-", colnames(pos.data)[profile.col], perl="TRUE")
colnames(pos.data)[profile.col] <- sub("X", "+", colnames(pos.data)[profile.col], perl="TRUE")

## Draw distribution of chi2 significance scores
x11(width=7, height=5)
hist(pos.data$sig, breaks=100, main="Distribution of chi2 significance", xlab="sig of chi2 test", ylab="Number of k-mers", col="#DDDDDD")
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'sig_distrib', suffix)), export.formats=export.formats.plots, width=7,height=5)

## Define the selected patterns
selected.patterns <- (pos.data$sig >= sig.threshold) & (pos.data$rank <= rank.threshold)
nb.patterns <- sum(selected.patterns)
verbose(paste(nb.patterns , "selected patterns"))

## Create a separate data frame for the position profiles.
## Redundant -> memory-inefficient but convenient.
pos.profiles <- pos.data[selected.patterns, profile.col]

## Normalize position profiles by row (convert k-mer occurrences in
## relative frequencies)
sum.per.kmer <- apply(pos.profiles, 1, sum)
pos.profiles.freq <- pos.profiles / sum.per.kmer
pos.profiles.freq.norm <- pos.profiles.freq - apply(pos.profiles.freq, 1, median)
export.object(pos.profiles.freq.norm, file=file.path(dir.clusters, paste(sep='_', 'position_profiles_norm_freq', suffix)), export.format='table')

## perform a random shuffling of the position profiles, as a negative
## control for correlation and clustering
shuffled.profiles <- as.data.frame(matrix(nrow=nb.patterns, ncol=nb.windows, sample(as.vector(as.matrix(pos.profiles, replace=FALSE)))))
rownames(shuffled.profiles) <- rownames(pos.profiles)
colnames(shuffled.profiles) <- colnames(pos.profiles)
shuffled.profiles.freq <- shuffled.profiles / sum.per.kmer
export.object(shuffled.profiles.freq, file=file.path(dir.clusters, paste(sep='_', 'shuffled_profiles_norm_freq', suffix)), export.format='table')

################################################################
## Correlation analysis

## Compute correlations between profiles
pos.cor <- cor(t(pos.profiles))
export.object(pos.cor, file=file.path(dir.clusters, paste(sep='_', 'profiles_correlations', suffix)), export.format='table')
shuffled.cor <- cor(t(shuffled.profiles))
export.object(shuffled.cor, file=file.path(dir.clusters, paste(sep='_', 'shuffled_correlations', suffix)), export.format='table')

## Draw the distributions of correlation values (position profiles + shuffled data)
x11(width=7, height=8)
par(mfrow=c(2,1))
hist(pos.cor, breaks=100, xlim=c(-1,1), col=colors["pos"], main="Correlations between position profiles", xlab="Correlation", ylab="Pairs of k-mers")
hist(shuffled.cor, breaks=100, xlim=c(-1,1), col=colors["shuffled"], main="Correlations between shuffled profiles", xlab="Correlation", ylab="Pairs of k-mers")
par(mfrow=c(1,1))
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'profile_correlations', suffix)), export.formats=export.formats.plots, width=7,height=8)

## Draw a heatmap of the correlation matrix
X11(width=12, height=12)
heatmap(pos.cor, scale="none",
        main='Correlations between position profiles',
        zlim=c(-1,1),
        col=heatmap.palette)
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'profile_correlations_heatmap', suffix)), export.formats=export.formats.plots, width=12,height=12)

## Draw a heatmap of the correlation matrix of shuffled profiles
X11(width=12, height=12)
heatmap(shuffled.cor, scale="none",
        main='Correlations between shuffled profiles',
        zlim=c(-1,1),
        col=heatmap.palette)
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'shuffled_correlations_heatmap', suffix)), export.formats=export.formats.plots, width=12,height=12)


################################################################
## Hierarchical clustering
clust.method <- "complete"
pos.tree <- hclust(as.dist(1-pos.cor), method=clust.method)
shuffled.tree <- hclust(as.dist(1-shuffled.cor), method=clust.method)

## Plot the tree structure for position profiles and shuffled data
x11(width=16, height=8)
par(mfrow=c(2,1))
plot(pos.tree, main=paste("Clustering of position profiles;", "dist= 1 - cor; ", clust.method, "linkage"), xlab=NA)
plot(shuffled.tree, main=paste("Clustering of shuffled profiles;", "dist= 1 - cor; ", clust.method, "linkage"), xlab=NA)
par(mfrow=c(1,1))
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'profile_tree', suffix)), export.formats=export.formats.plots, width=16,height=8)


## Draw a heatmap of the clustered position profiles
X11(width=12, height=12)
zmax <- max(abs(range(pos.profiles.freq.norm)))
heatmap(as.matrix(pos.profiles.freq.norm),
        scale="none",
        main=paste("Clustered position profiles; ", 'd=1-cor; ', clust.method, "linkage"),
        zlim=c(-zmax,zmax),
        col=green.white.red(),
        xlab='position relative to summit', 
        ylab='k-mers', 
        Colv=NA, Rowv=as.dendrogram(pos.tree))
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'profile_heatmap', suffix)), export.formats=export.formats.plots, width=12,height=12)


## Cut the tree in k clusters
clusters <- cutree(pos.tree,k=nb.clusters)
write.table(as.data.frame(clusters), file=file.path(dir.clusters, paste(sep='', 'clusters_', suffix, '.tab')), row.names=TRUE, col.names=FALSE, quote=FALSE)
table(clusters)

## Plot profiles for all clusters
x11(width=8, height=16)
par(mfrow=c(nb.clusters,1))
for (i in 1:nb.clusters) {
  cluster.size <- sum(clusters==i)
  plot.profiles(as.data.frame(pos.profiles.freq.norm[clusters==i,]), ,
                main=paste(sep='', 'cluster ', i, '/', nb.clusters),
                col.profiles=rainbow(n=cluster.size),
                plot.median.profile=F, plot.mean.profile=T, plot.sd.profile=F
                )
}
par(mfrow=c(1,1))
export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'profiles_all_clusters', suffix)), export.formats=export.formats.plots, width=8, height=16)

## Plot profiles for each cluster
x11(width=14, height=5)
par(cex=0.7)
par(cex.lab=0.7)
for (i in 1:length(table(clusters))) {
  cluster.size <- sum(clusters==i)
  plot.profiles(as.data.frame(pos.profiles.freq.norm[clusters==i,]), ,
                main=paste(sep='', 'cluster ', i, '/', nb.clusters),
                col.profiles=rainbow(n=cluster.size),
                plot.median.profile=F, plot.mean.profile=T, plot.sd.profile=F
                )
  export.plot(file.prefix=file.path(dir.clusters, paste(sep='_', 'profiles_cluster', i, suffix)), export.formats=export.formats.plots, width=14, height=5)
}

verbose ("job done")
verbose(dir.clusters)
