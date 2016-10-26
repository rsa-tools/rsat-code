################################################################
## At this stage you can already check some simple RSAT command 

## Test a simple Perl script that does not require for organisms to be
## installed.
which random-seq
random-seq -l 100

## Test a simple python script that does not require organisms to be
## installed.
random-motif -l 10 -c 0.90

################
## Test some external programs

## vmatch (used in purge-sequence)
random-seq -l 100 | purge-sequence

## Check that seqlogo is installed
which seqlogo
seqlogo

## Check that weblogo 3 is installed
which weblogo
weblogo --help

## ghostscript
which gs
gs --version

## Check tat the model genomes have been correctly installed
## A simple and quick test: retrieve all the start codons and count
## oligonucleotide frequencies (most should be ATG).
retrieve-seq -org Saccharomyces_cerevisiae -all -from 0 -to +2 \
    | oligo-analysis -l 3 -1str -return occ,freq -sort

## Check that the Perl graphical libraries are properly installed
XYgraph -help
feature-map -help

## Test the program supported-organisms-server, which relies on Web
## services without stub
supported-organisms-server -url http://${IP}/rsat/ | wc
supported-organisms-server -url http://localhost/rsat/ | wc
supported-organisms-server -url http://rsat-tagc.univ-mrs.fr/ | wc
