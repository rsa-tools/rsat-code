## ##############################################################
## Analyse pattern count distribution obtained with
## calibrate-oligos.pl
source('~/rsa-tools/R-scripts/negbin.R')
source("~/research/R-files/config.R")
setwd(dir.util)
source('util_test_fitting_poisson.R')

## ##############################################################
## Initialize variables
# organism <- "Homo_sapiens"
#organism <- "Homo_sapiens"
organisms <- c("Homo_sapiens", "Saccharomyces_cerevisiae")

organism <- "Homo_sapiens"
repet <- 10000
seq.nb <- 10
seq.len <- 500

organism <- "Saccharomyces_cerevisiae"
repet <- 10000
seq.nb <- 11
seq.len <- 1000


oligo.len <- 6
str <- '-2str'
ov <- '-noov'

## Directories
dir.main <- "~/motif_discovery_competition_2003"

dir.calib <- paste(dir.main, 
		    "results", 
		    organism, 
		    "rand_gene_selections", 
		    paste(oligo.len, "nt", str, ov, '_N', seq.nb, '_L', seq.len, '_R', repet, sep=''),
		    sep="/")
setwd(dir.calib)

## files
file.distrib <- paste(organism, "_",oligo.len,"nt_", str, ov, '_n',seq.nb,'_l',seq.len,'_r',repet,'_distrib.tab', sep="")

file.stats <- paste(organism, "_",oligo.len,"nt_", str, ov, '_n',seq.nb,'_l',seq.len,'_r',repet,'_stats.tab', sep="")

file.fitting.poisson <- paste(organism, "_",oligo.len,"nt_", str, ov, '_n',seq.nb,'_l',seq.len,'_r',repet,'_fitting_poisson.tab', sep="")

file.fitting.negbin <- paste(organism, "_",oligo.len,"nt_", str, ov, '_n',seq.nb,'_l',seq.len,'_r',repet,'_fitting_negbin.tab', sep="")


## pattern count distribution
distrib <- read.table(file.distrib, sep="\t", comment.char=";",row.names = 1)
values <- 0:(ncol(distrib)-1)
names(distrib) <- values

theor <- "poisson"
file.fitting <- paste(organism, "_",oligo.len,"nt_", str, ov, '_n',seq.nb,'_l',seq.len,'_r',repet,'_', theor, '.tab', sep="")

header <- c('avg',
	    'std',
	    'var',
	    'repet',
	    'fitted.distrib',
	    'chi2',
	    'df',
	    'left.group',
	    'right.group',
	    'obs.grouped',
	    'exp.grouped');
if (theor == 'negbin') {
    header <- c(header, 'p', 'k')
}

## Read statistics calculated by calibrate-oligos.pl
#stats <- read.table(file.stats, sep="\t", comment.char=";",row.names = 1)
fitting <- read.table(file.fitting, sep="\t", comment.char=";",row.names = 1)
names(fitting) <- header

## Calculate P-value, E-value and significance
fitting$p.value <- pchisq(fitting$chi2, fitting$df,lower.tail=F)
fitting$e.value <- fitting$p.value * nrow(fitting)
fitting$sig <- -log(fitting$e.value, base=10)

## Sort the table according to the P-value
fitting.sorted <- fitting[order(fitting$p.value),]

## Report the most significant patterns
fitting.sorted[1:10,c("avg","chi2","df","p.value","e.value","sig")]

## Report the less significant patterns
fitting.sorted[(nrow(fitting.sorted)-9):nrow(fitting.sorted),c("avg","chi2","df","p.value","e.value","sig")]

## Count number of significant motifs
table(fitting$sig >= 0)
    
## count the number of patterns for which the variance is greater than the mean
table(fitting$var > fitting$avg)

## plot a graph with variance as a function of avg
plot(fitting$avg,fitting$var,
     main=file.distrib,
     xlab='average',
     ylab='variance',
     log="xy",
     panel.first=c(grid(col="blue"),
		   abline(a=0,b=1,col="red")
		   )
     )


## Select the most significant patterns
selected.patterns <- row.names(fitting.sorted)[1:16]

## select the less significant patterns
# selected.patterns <- rev(row.names(fitting.sorted))[1:16]

## Plot the distrribution of the selected patterns
par(mfrow=c(4,4))
for (pattern in selected.patterns) {
    test.fitting.poisson.histo(values,distrib[pattern,],subtitle=paste(toupper(pattern), fitting[pattern, "sig"]))
    lines(values,fitting[pattern,"repet"]*dnegbin(values, m=fitting[pattern,"avg"], v=fitting[pattern,"var"]), col="#00AA00",lwd=2)
}
par(mfrow=c(1,1))

#  ################################################################
#  ## Calculate statistics from the distribution


stats.2 <- data.frame(row.names=row.names(distrib), 
		      sum=rep(NA,nrow(distrib)),
		      ssq=rep(NA,nrow(distrib)),
		      avg=rep(NA,nrow(distrib)),
		      var=rep(NA,nrow(distrib)),
		      std=rep(NA,nrow(distrib))
		      )
stats.2$sum <-  t(t(as.matrix(values)) %*% t(as.matrix(distrib)))
stats.2$ssq <-  t(t(as.matrix(values^2)) %*% t(as.matrix(distrib)))
stats.2$avg <- stats.2$sum/repet
stats.2$var <- stats.2$ssq/repet - (stats.2$avg)^2		
stats.2$std <- sqrt(stats.2$var)


## plot a graph with std as a function of avg
plot(stats.2$avg,stats.2$var,
     main=file.distrib,
     xlab='average',
     ylab='standard deviation',
     log="xy",
     panel.first=c(grid(col="blue"),
		   abline(a=0,b=1,col="red")
		   )
     )




