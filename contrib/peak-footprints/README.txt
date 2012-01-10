CONFIGURATION
=============

...

INSTALLING MAF FILES
===================

The current version of peak-footpints requires to download the
complete set of multiple alignment files (MAF) from the UCSC Genome
Browser, for a selected set of species and genome versions.

Downloading the MAFs
--------------------

- Open a connection to UCSC server (http://genome.ucsc.edu/)

- Click on "Downloads"

- On the Download page, click on the reference organism (the organism
  from which the ChIP-seq peaks were obtained e.g. Mouse)

- Choose a genome version (e.g. mm8)
   Under the title "Multiple Alignments", follow the link 
    "Multiple alignments of NN vertebrate genomes with XXX"
    (e.g. Multiple alignments of 16 vertebrate genomes with Mouse)

- Download all the files except those containing the suffix
  "random.maf.gz" and uncompress them with the command gunzip.

Indexing MAFs
-------------

- ...
