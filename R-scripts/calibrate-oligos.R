## ##############################################################
## Analyse pattern count distribution obtained with
## calibrate-oligos.pl
#dir.rsat <- Sys.getenv("RSAT")
dir.rsat <- file.path('~', "rsa-tools")
source(file.path(dir.rsat, 'R-scripts', 'negbin.R'))
source(file.path(dir.rsat, 'R-scripts', 'distrib_fitting.R'))


## ##############################################################
## Initialize variables
# organism <- "Homo_sapiens"
#organism <- "Homo_sapiens"
organisms <- c("Homo_sapiens", "Saccharomyces_cerevisiae", "Plasmodium_falciparum", "Drosophila_melanogaster")

organism <- "Homo_sapiens"
organism <- "Saccharomyces_cerevisiae"
organism <- "Plasmodium_falciparum"

repet <- 1000
seq.nb <- 10
seq.len <- 1000


oligo.len <- 6
str <- '-2str'
ov <- '-noov'
theor <- "poisson"
  
N.values<-c(1,2,3,4,5,10,20,30,50,100,150,200)

## Directories
dir.main <- "~/Documents/biosapiens"

################################################################
## Run calibration for one set of parameters
one.calibration <- function(organism="Saccharomyces_cerevisiae",
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
                         "results", 
                         organism, 
                         "rand_gene_selections", 
                         paste(oligo.len, "nt", str, ov, '_N', seq.nb, '_L', seq.len, '_R', repet, sep=''))

  verbose(dir.calib)
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
  result <- analyze.fitting (file.prefix=file.prefix, theor=theor, dir.input = dir.calib, ...)
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


one.N.series(theor=theor,organism="Plasmodium_falciparum",seq.len=1500,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)
one.N.series(theor=theor,organism="Plasmodium_falciparum",seq.len=1000,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)
one.N.series(theor=theor,organism="Saccharomyces_cerevisiae",seq.len=600,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)
one.N.series(theor=theor,organism="Homo_sapiens_EnsEMBL",seq.len=2000,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)

