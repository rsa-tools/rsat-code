## ##############################################################
## Analyse pattern count distribution obtained with
## calibrate-oligos.pl
dir.rsat <- Sys.getenv("RSAT")
#dir.rsat <- file.path('~', "rsa-tools")
source(file.path(dir.rsat, 'R-scripts', 'negbin.R'))
source(file.path(dir.rsat, 'R-scripts', 'distrib_fitting.R'))

## Directories
dir.main <- "~/Documents/research/biosapiens"

## ##############################################################
## Initialize variables
# organism <- "Homo_sapiens"
#organism <- "Homo_sapiens"
organisms <- c("Homo_sapiens", "Saccharomyces_cerevisiae", "Plasmodium_falciparum", "Drosophila_melanogaster")

#organism <- "Homo_sapiens"
#organism <- "Saccharomyces_cerevisiae"
#organism <- "Plasmodium_falciparum"
organism <- "Drosophila_melanogaster"
repet <- 1000
seq.len <- 1000
seq.nb <- 10


oligo.len <- 6
str <- '-2str'
ov <- '-noov'
theor <- "poisson"
  
N.values<-c(1,2,3,4,5,10,20,30,50,100,150,200)


one.N.series(theor=theor,organism="Saccharomyces_cerevisiae",seq.len=600,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)
one.N.series(theor=theor,organism="Plasmodium_falciparum",seq.len=1500,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)
one.N.series(theor=theor,organism="Plasmodium_falciparum",seq.len=1000,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)
one.N.series(theor=theor,organism="Homo_sapiens_EnsEMBL",seq.len=2000,export.small.plots=F,close.windows=T, N.values=N.values,export.plots=F)


################################################################
## Analysis for Florence Mougel
dir.main <- "~/Documents/research/collaborations/florence_mougel/"
repet <- 99
seq.len <- 5000
seq.nb <- 10

one.N.series(subdir="oligo_calibrations", theor=theor,organism="Drosophila_melanogaster",seq.len=5000,export.small.plots=F,close.windows=T, N.values=10,export.plots=F, repet=repet)
