################################################################
## Analyze the results of upstream calibration

dir.rsat <- file.path('~', "rsa-tools")
source(file.path(dir.rsat, 'R-scripts', 'negbin.R'))
source(file.path(dir.rsat, 'R-scripts', 'distrib_fitting.R'))

org <- 'Saccharomyces_cerevisiae_no_mito'
seq.len <- 1000
oligo.len <- 6
str <- '-1str'
ov <- '-noov'

dir.upstream.calibrations <- file.path('~', 'research', 'upstream_calibrations')

upstream.calibration <- function (org='Homo_sapiens', 
				  seq.len = 1000,
				  oligo.len = 6,
				  str='-1str',
				  ov='-noov') {
  dir.org <- file.path(dir.upstream.calibrations, 'results', org)
  dir.distrib <- file.path(dir.org, 'oligo_distrib')
  setwd(dir.distrib) 
  
  theor <- 'negbin'
  for (theor in c('negbin', 'poisson')) {
    file.prefix <- paste(org, '_allup', seq.len, '_', oligo.len, 'nt', str, ov, sep="")
    analyze.fitting(file.prefix=file.prefix, theor=theor, dir.input=dir.distrib, dir.output=dir.distrib)
  }
}


upstream.calibration(org='Saccharomyces_cerevisiae_no_mito', seq.len=800)
upstream.calibration(org='Homo_sapiens', seq.len=1000)
upstream.calibration(org='Arabidopsis_thaliana', seq.len=1000)
upstream.calibration(org='Escherichia_coli_K12', seq.len=200)
