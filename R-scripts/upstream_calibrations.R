################################################################
## Analyze the results of upstream calibration

dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {dir.rsat <- file.path('~', "rsa-tools")}
source(file.path(dir.rsat, 'R-scripts', 'negbin.R'))
source(file.path(dir.rsat, 'R-scripts', 'distrib_fitting.R'))

## Groups of organisms
fungi <- c('Saccharomyces_cerevisiae','Schizosaccharomyces_pombe','Ashbya_gossypii')
procaryotes <- c('Mycoplasma_genitalium','Escherichia_coli_K12','Bacillus_subtilis','Salmonella_typhimurium_LT2')
other.orgs <- c('Homo_sapiens','Drosophila_melanogaster','Caenorhabditis_elegans','Arabidopsis_thaliana','Plasmodium_falciparum')


## Default parameters
org <- 'Saccharomyces_cerevisiae'
seq.len <- 800
oligo.len <- 6
str <- '-1str'
ov <- '-noov'

result.table <- matrix(ncol=7,nrow=0)
export.plots <- F

dir.upstream.calibrations <- file.path('~', 'research', 'upstream_calibrations')


upstream.calibration <- function (org='Saccharomyces_cerevisiae', 
				  seq.len = 800,
				  oligo.len = 6,
				  str='-1str',
				  ov='-noov',
                                  result.table,
                                  ...) {
  dir.org <- file.path(dir.upstream.calibrations, 'results', org)
  dir.distrib <- file.path(dir.org, 'oligo_distrib')
  setwd(dir.distrib) 
  
  theor <- 'negbin'
  for (theor in c('poisson','negbin')) {
    file.prefix <- paste(org, '_allup', seq.len, '_', oligo.len, 'nt', str, ov, sep="")
    result <- analyze.fitting(file.prefix=file.prefix, theor=theor, dir.input=dir.distrib, dir.output=dir.distrib,...)
    result.table <- rbind(result.table,t(as.matrix(result)))
  }
  return (result.table)
}


for (ol in 1:6) {
  for (org in fungi) {
    result.table <- upstream.calibration(org=org, oligo.len=ol, seq.len=800,export.plots=export.plots,result.table=result.table)
  }
  for (org in procaryotes) {
    result.table <- upstream.calibration(org=org, oligo.len=ol, seq.len=200,export.plots=export.plots,result.table=result.table)
  }
  for (org in other.orgs) {
    for (len in c(200,500,1000,2000)) {
      result.table <- upstream.calibration(org=org, oligo.len=ol, seq.len=len,export.plots=export.plots,result.table=result.table)
    }
  }
}

for (ol in 1:6) {  
# result.table <- upstream.calibration(org='Ashbya_gossypii', seq.len=800,export.plots=export.plots,result.table=result.table)
  result.table <- upstream.calibration(org='Saccharomyces_cerevisiae', oligo.len=ol, seq.len=800,export.plots=export.plots,result.table=result.table)
  result.table <- upstream.calibration(org='Homo_sapiens', oligo.len=ol, seq.len=200,export.plots=export.plots,result.table=result.table)
  result.table <- upstream.calibration(org='Drosophila_melanogaster', oligo.len=ol, seq.len=200,export.plots=export.plots,result.table=result.table)
# result.table <- upstream.calibration(org='Arabidopsis_thaliana', seq.len=1000,export.plots=export.plots,result.table=result.table)
# result.table <- upstream.calibration(org='Escherichia_coli_K12', seq.len=200,export.plots=export.plots,result.table=result.table)
}

setwd(dir.results); export.object(result.table,file.prefix='fitting_tests_summary',export.formats='table')

