
source("~/research/R-files/config.R")
source(file.path(dir.util, 'util.R'))
source(file.path(dir.util, 'util_test_fitting_poisson.R'))

################################################################
## Calculate P-value, E-value and sig for a fitting file
## and draw some plots
analyze.fitting <- function(file.prefix, ## file resulting from the perl script fit-distribution
			    theor, ## theoretcial distribution; supported: negbin, poisson
			    dir.input='.', ## Input directory (containing the input file)
			    dir.output='.', ## Directory where the fitting files have to be stored)
			    dir.figures=NA ## Directorb where the figrues have to be stored
			    ) {


  if (is.na(dir.figures)) {
      dir.figures <- file.path(dir.output, 'figures')
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
  header <- c('id', 
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
      
  ## count the number of patterns for which the variance is greater than the mean
  verbose(paste(theor, 'fitting; paterns with var > avg'))
  print(table(fitting$var > fitting$avg))
  
  ## plot a graph with variance as a function of avg
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

  
  ## Select the most significant patterns
  selected.patterns <- row.names(fitting.sorted)[1:16]
  
  ## select the less significant patterns
  # selected.patterns <- rev(row.names(fitting.sorted))[1:16]
  
  ## pattern count distribution
  setwd(dir.input); distrib <- read.table(file.prefix, sep="\t", comment.char=";",row.names = 1)
  distrib <- distrib[, c(2:ncol(distrib))]
  values <- 0:(ncol(distrib)-1)
  names(distrib) <- values


  ## Plot the distribution of the selected patterns - one image per word
  for (pattern in selected.patterns) {
      test.fitting.poisson.histo(values,distrib[pattern,],subtitle=paste(toupper(pattern), fitting[pattern, "sig"]))
      lines(values,fitting[pattern,"repet"]*dnegbin(values, m=fitting[pattern,"avg"], v=fitting[pattern,"var"]), col="#008800",lwd=2)
      setwd(dir.figures); export.plot(file.prefix=paste(file.prefix, theor, 'fitting', toupper(pattern), sep='_'), export.formats=export.formats.plots, width=8,height=8)
  }

  ## Plot the distribution of the selected patterns - one big image with 16 frames
  par(mfrow=c(4,4))
  for (pattern in selected.patterns) {
      test.fitting.poisson.histo(values,distrib[pattern,],subtitle=paste(toupper(pattern), fitting[pattern, "sig"]))
      lines(values,fitting[pattern,"repet"]*dnegbin(values, m=fitting[pattern,"avg"], v=fitting[pattern,"var"]), col="#008800",lwd=2)
  }
  par(mfrow=c(1,1))
  setwd(dir.figures); export.plot(file.prefix=paste(file.prefix, 'fitting', 'top16', sep='_'), export.formats=export.formats.plots, width=12,height=9)
}
