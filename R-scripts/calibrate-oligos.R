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
organisms <- c("Homo_sapiens", "Saccharomyces_cerevisiae")

organism <- "Homo_sapiens"
organism <- "Saccharomyces_cerevisiae"

repet <- 10000
seq.nb <- 10
seq.len <- 500


oligo.len <- 6
str <- '-2str'
ov <- '-noov'

## Directories
dir.main <- "~/motif_discovery_competition_2003"
    
dir.calib <- file.path(dir.main, 
		       "results", 
		       organism, 
		       "rand_gene_selections", 
		       paste(oligo.len, "nt", str, ov, '_N', seq.nb, '_L', seq.len, '_R', repet, sep=''))
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
distrib <- read.table(file.distrib, sep="\t", comment.char=";",row.names = 1)
values <- 0:(ncol(distrib)-1)
names(distrib) <- values

theor <- "negbin"
theor <- "poisson"

for (theor in c('poisson', 'negbin')) {
  file.fitting <- paste(file.prefix,'_', theor, '.tab', sep="")
  analyze.fitting (file.prefix=file.prefix, theor=theor, dir.input = dir.calib)

  
}


