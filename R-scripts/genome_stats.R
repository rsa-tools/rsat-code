################################################################
# Plot some results of genome statistics
# This script should be used after the makefile
# ${RSAT}/makefiles/genome_statistics.mk


## source('~/Documents/enseignement/bioinformatics_courses/bioinformatics_introductory_course/R-scripts/genome_stats.R')

source("~/research/R-files/config.R")
setwd(dir.util); source('util.R')

dir.main <- file.path(dir.home,"Documents","enseignement","bioinformatics_courses","bioinformatics_introductory_course")


################################################################
## Read genome statistics
dir.data <- file.path(dir.main, "data","genome_stats")
setwd(dir.data)
dir.results <- dir.data
dir.figures <- file.path(dir.results, "figures")
file.genome.stats <- file.path(dir.data, "stats_all_organisms.tab")
x.ori <- read.table(file.genome.stats, header=T,row.names=1,sep="\t")

## ignore incomplete genomes from WASHU
x.ori <- x.ori[-grep("WASHU",row.names(x.ori)),]

## sort taxonomy by alphabetical order
#o1 <- order(as.vector(x.ori$taxonomy))
#x<- x.ori[o1,]

## Sort by taxonomy, bacteria first (to avoid hiding archea on the plots)
b <- grep("Bacteria",as.vector(x.ori$taxonomy))
a <- grep("Archaea",as.vector(x.ori$taxonomy))
e <- grep("Eukaryota",as.vector(x.ori$taxonomy))
f <- grep("Fungi",as.vector(x.ori$taxonomy))
o2 <- c(b,a,e)
x <- x.ori[o2,]

## Calculate gene density
x$density <- 1000*x$genes/x$size
x$gene.spacing <- (x$size/x$gene)/1000

attach(x)

## Associate colors to the different taxa
tax.palette = c(
  "Bacteria"=rgb(0,0.8,0.8),
  "Archaea"=rgb(0,0,1),
  "Fungi"=rgb(1,0.5,0),
  "Viridiplantae"=rgb(0,0.7,0),
  "Metazoa"=rgb(0.5,0,1),
  "Mammalia"=rgb(1,0,1),
  "Alveolata"=rgb(0.5,0.5,0),
  "Microsporidia"=rgb(1, 0,0)
  )
tax.colors <- rep(rgb(0,0,0), nrow(x))
names(tax.colors) = row.names(x)
for (taxon in names(tax.palette)) {
  tax.colors[grep(taxon, as.vector(taxonomy))] <- tax.palette[taxon]
}

## Report organisms whose taxonomy is not in the palette
print(as.vector(taxonomy[tax.colors==rgb(0,0,0)]))

################################################################
## Plot number of genes as a function of genome size 
plot(size, genes,
     log="",
     xlab="genome size",
     ylab="number of genes",
     xlim=c(0,1.5e8),
     col=tax.colors,
     panel.first=(grid(col='#000000')),
     pch=19
     )
legend("bottomright",legend=names(tax.palette),col=tax.palette,pch=20,bg="white",bty="o")
setwd(dir.figures); export.plot('size_vs_genes',export.formats="pdf",width=8,height=8)

## with log scale
plot(size, genes,log="xy",
     xlab="genome size (log scale)",
     ylab="number of genes (log scale)",
     col=tax.colors,
     panel.first=(grid(col='#000000')),
     pch=19
     )
legend("bottomright",legend=names(tax.palette),col=tax.palette,pch=20,bg="white",bty="o")
setwd(dir.figures); export.plot('size_vs_genes_log',export.formats="pdf",width=8,height=8)

################################################################
## Fraction of intergenic sequences as a function of genome size
plot(size, fract.1,log="x",
     xlab="genome size (log scale)",
     ylab="intergenic fraction",     
     ylim=c(0,1),
     col=tax.colors,
     panel.first=(grid(col='#000000')),
     pch=19
     )
legend("topleft",legend=names(tax.palette),col=tax.palette,pch=20,bg="white",bty="o")
setwd(dir.figures); export.plot('size_vs_intergenic_fraction',export.formats="pdf",width=8,height=8)

################################################################
## Gene spacing as a function of genome size
plot(size, gene.spacing,log="xy",
     xlab="genome size (log scale)",
     ylab="gene spacing (kb/gene) (log)",     
     col=tax.colors,
     panel.first=(grid(col='#000000')),
     pch=19
     )
legend("topleft",legend=names(tax.palette),col=tax.palette,pch=20,bg="white",bty="o")
setwd(dir.figures); export.plot('size_vs_gene_spacing',export.formats="pdf",width=8,height=8)

detach(x)

################################################################
## Read protein statistics
dir.proteins <- file.path(dir.main, "data","protein_stats")
setwd(dir.proteins)
file.protein.stats <- file.path(dir.proteins, "protein_stats_all_organisms.tab")
p.ori <- read.table(file.protein.stats, header=T,sep="\t")
row.names(p.ori) <- p.ori[,"organism_name"]
dir.results <- dir.proteins
dir.figures <- file.path(dir.results, "figures")

p <- p.ori[(p.ori$X.n > 400 & p.ori$mean > 99),]

orgs <- intersect(row.names(p),row.names(x))

plot(x[orgs,"size"],
     p[orgs,"mean"],
     log="x",
     xlab='genome size (log scale)',
     ylab='average protein size',
     ylim=c(0,max(p[orgs,"mean"])),
     col=tax.colors[orgs],
     panel.first=(grid(col='#000000')),
     pch=19
     )
legend("topleft",legend=names(tax.palette),col=tax.palette,pch=20,bg="white",bty="o")
setwd(dir.figures); export.plot('genome_vs_protein_sizes',export.formats="pdf",width=8,height=8)
