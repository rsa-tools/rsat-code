source("~/research/R-files/config.R")
source(file.path(dir.util, 'util.R'))
source(file.path(dir.util, 'util_test_fitting_poisson.R'))

################################################################
## Calculate P-value, E-value and sig for a fitting file and draw
## distribution histograms for a selected of oligonucleotides
analyze.fitting <- function(file.prefix, ## file resulting from the perl script fit-distribution
			    theor, ## theoretcial distribution; supported: negbin, poisson
			    dir.input='.', ## Input directory (containing the input file)
			    dir.output=NA, ## Directory where the fitting files have to be stored)
			    dir.figures=NA, ## Directorb where the figrues have to be stored
                            export.plots=T, ## Export plots
                            export.small.plots=T, ## Export plots
                            close.windows=F
                            ) {

  
  ## generate an object with the result
  result <- list()
  result$file.prefix <- file.prefix
  result$theor <- theor
  result$dir.input <- dir.input
  if (is.na(dir.output)) {
    dir.output <- dir.input
  }
  result$dir.output <- dir.output
  
  ## Check output directories
  if (is.na(dir.figures)) {
    dir.figures <- file.path(dir.output, 'figures')
    dir.create(dir.figures, showWarnings=F)
  }
  verbose (dir.output)
  for (dir in c(dir.output, dir.figures)) {
    if (!file.exists(dir)) {
      dir.create(dir)
    }
  }

  file.fitting <- paste(file.prefix, '_', theor, '.tab', sep='')
  verbose (file.fitting)

  ## Column headers
  header <- c(
	      'avg',
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
  setwd(dir.input) 
  fitting <- read.table(file.fitting, sep="\t", comment.char=";",row.names = 1)
  names(fitting) <- header
  
  ################################################################
  ## Calculate P-value, E-value and significance
  fitting$p.value <- pchisq(fitting$chi2, fitting$df,lower.tail=F)
  fitting$e.value <- fitting$p.value * nrow(fitting)
  fitting$sig <- -log(fitting$e.value, base=10)

  result$nrow <- nrow(fitting)
  result$fitting <- fitting
  
  ## Sort the table according to the P-value
  fitting.sorted <- fitting[order(fitting$p.value),]

  ## Export the resulting table
  setwd(dir.output); export.object(fitting.sorted, file.prefix=paste(file.prefix, theor, 'sig', sep='_'), export.formats='table')
  
  ## Report the most significant patterns
  fitting.sorted[1:10,c("avg","chi2","df","p.value","e.value","sig")]
  
  ## Report the less significant patterns
  fitting.sorted[(nrow(fitting.sorted)-9):nrow(fitting.sorted),c("avg","chi2","df","p.value","e.value","sig")]
  
  ## Count number of significant motifs
  verbose(paste(theor, 'fitting; paterns with sig > 0'))
  print(table(fitting$sig >= 0))
  result$rows.fitting <- sum(fitting$sig < 0)
      
  ## count the number of patterns for which the variance is greater than the mean
  verbose(paste(theor, 'fitting; paterns with var > avg'))
  print(table(fitting$var > fitting$avg))
  result$rows.var.gt.avg <- sum(fitting$var > fitting$avg)
  

  ## ##############################################################
  ## Plot some distributions
  if (export.plots) {
  
    ## plot a graph with variance as a function of avg
    X11(width=8,height=8)
    plot(fitting$avg,fitting$var,
         main=file.prefix,
                                        #       sub='variance versus average',
         xlab='average',
         ylab='variance',
         log="xy",
         panel.first=c(grid(col="blue"),
           abline(a=0,b=1,col="red")
           )
         )
    setwd(dir.figures); export.plot(file.prefix=paste(file.prefix, 'avg', 'var', sep='_'), export.formats=export.formats.plots, width=8,height=8)
    if (close.windows) {dev.off()}
    
    ## Select the 16 most significant patterns
    selected.patterns <- row.names(fitting.sorted)[1:16]
    
    ## select the less significant patterns
    ## selected.patterns <- rev(row.names(fitting.sorted))[1:16]
    
    ## pattern count distribution
    file.distrib <- paste(file.prefix, 'distrib.tab', sep='_')
    setwd(dir.input); distrib <- read.table(file.distrib, sep="\t", comment.char=";",row.names = 1)
    distrib <- distrib[, c(2:ncol(distrib))]
    values <- 0:(ncol(distrib)-1)
    names(distrib) <- values
    
    ## Plot the distribution of the selected patterns - one big image with 16 frames
    X11(width=12,height=9)
    par(mfrow=c(4,4))
    for (pattern in selected.patterns) {
      test.fitting.poisson.histo(values,distrib[pattern,],subtitle=paste(toupper(pattern), "; chi2=",round(fitting[pattern, "sig"],digits=2)))
      grid(col='black')
      lines(values,fitting[pattern,"repet"]*dnegbin(values, m=fitting[pattern,"avg"], v=fitting[pattern,"var"]), col="#008800",lwd=2)
    }
    par(mfrow=c(1,1))
    setwd(dir.figures); export.plot(file.prefix=paste(file.prefix, 'fitting', 'top16', sep='_'), export.formats=export.formats.plots, width=12,height=9)
    if (close.windows) {dev.off()}

        
    if (export.small.plots) {      
      ## Plot the distribution of the selected patterns - one image per word
      X11(width=8,height=8)
      for (pattern in selected.patterns) {
        test.fitting.poisson.histo(values,distrib[pattern,],subtitle=paste(toupper(pattern), "; chi2=",round(fitting[pattern, "sig"],digits=2)))
        grid(col='black')
        lines(values,fitting[pattern,"repet"]*dnegbin(values, m=fitting[pattern,"avg"], v=fitting[pattern,"var"]), col="#008800",lwd=2)
        setwd(dir.figures); export.plot(file.prefix=paste(file.prefix, theor, 'fitting', toupper(pattern), sep='_'), export.formats=export.formats.plots, width=8,height=8)
      }
      dev.off()
    }

  }

  return(result)
}



################################################################
## Run calibration for one set of parameters
one.calibration <- function(subdir="results",
                            organism="Saccharomyces_cerevisiae",
                            oligo.len=6,
                            seq.nb=10,
                            repet=1000,
                            seq.len=600,
                            str='-2str',
                            ov='-noov',
                            theor='poisson',
                            ...
                            ) {

  dir.calib <- file.path(dir.main, 
                         subdir,
                         organism, 
                         "rand_gene_selections", 
                         paste(oligo.len, "nt", str, ov, '_N', seq.nb, '_L', seq.len, '_R', repet, sep=''))


  verbose(paste("dir.calib", dir.calib))
  setwd(dir.calib)

  dir.results <- dir.calib
  dir.figures <- file.path(dir.calib, 'figures')
  dir.create.check(dir.figures)
  
  
  ## files
  file.prefix <- paste(organism, "_",oligo.len,"nt_", str, ov, '_n',seq.nb,'_l',seq.len,'_r',repet, sep='')
  file.distrib <- paste(file.prefix,'_distrib.tab', sep="")
  file.stats <- paste(file.prefix,'_stats.tab', sep="")
  file.fitting.poisson <- paste(file.prefix,'_fitting_poisson.tab', sep="")
  file.fitting.negbin <- paste(file.prefix,'_fitting_negbin.tab', sep="")
  
  
  ## pattern count distribution
                                        #distrib <- read.table(file.distrib, sep="\t", comment.char=";",row.names = 1)
                                        #values <- 0:(ncol(distrib)-1)
                                        #names(distrib) <- values
  
  file.fitting <- paste(file.prefix,'_', theor, '.tab', sep="")
  result <- analyze.fitting(file.prefix=file.prefix, theor=theor, dir.input = dir.calib, ...)
  return(result)
}


## Iterate over lengths
one.N.series <- function(N.values,
                         organism="Saccharomyces_cerevisiae",
                         oligo.len=6,
                         repet=1000,
                         seq.len=600,
                         str='-2str',
                         ov='-noov',
                         theor='poisson',
                         ...) { 
  for (theor in c('negbin', 'poisson')) {
    
    series.result <- list()
    summary.table <- matrix(nrow=0,ncol=5)
    i <- 0
    for (N in N.values) {
      i <- i+1
      result <- one.calibration(seq.nb=N,
                                organism=organism,
                                oligo.len=oligo.len,
                                repet=repet,
                                seq.len=seq.len,
                                str=str,
                                ov=ov,
                                theor=theor,
                                ...)
      summary.table <- rbind(summary.table, data.frame(result$file.prefix,
                                                       result$theor,
                                                       N,
                                                       result$rows.fitting,
                                                       result$rows.var.gt.avg
                                                       )
                             )
    }
    names(summary.table) <- c("file","theor","N","rows.fitting", "rows.var.gt.avg")  
    
    series.result$summary.table <- summary.table
    series.prefix <- file.prefix <- paste(organism, "_",oligo.len,"nt_", str, ov, '_l',seq.len,'_r',repet, "_", theor, sep='')
    series.result$prefix <- series.prefix
    
    
    setwd(file.path(dir.main, 'results', organism, 'rand_gene_selections')); export.object(series.result$summary.table, file=paste(series.result$prefix, 'summary', sep='_'),  export.format='table')
  }    
##    return(series.result)
}

