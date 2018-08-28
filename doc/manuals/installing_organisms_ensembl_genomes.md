---
title: "How to install organisms from Ensembl Genomes"
author: "Bruno Contreras Moreira"
date: 'Last update: 2018/06/19'
output:
  html_document:
    fig_caption: yes
    highlight: zenburn
    theme: cerulean
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    fig_caption: yes
    highlight: zenburn
    toc: yes
    toc_depth: 2
---

# Introduction

This document explains how to install genomes and annotations from Ensembl Genomes in RSAT. 
Please note that while Ensembl covers Vertebrates, Ensembl Genomes includes the rest 
(Protists, Fungi, Plants, Bacteria, Metazoa). These instructions have been tested with Ensembl Plants.
This manual uses *makefiles/ensemblgenomes_FTP_client.mk* for these tasks.

This document does NOT use the following scripts: install-ensembl-genome, download-ensembl-genome,
download-ensembl-features, download-ensembl-variations .


# Downloading a genome from another RSAT web server

Usually the fastest way of getting a genome installed is to fetch it from other RSAT web servers, 
as oligo frequencies are not computed; instead thay are copied over:

```{r, engine='bash', eval=FALSE}
download-organism -server http://rsat.eead.csic.es/plants -org Arabidopsis_thaliana.TAIR10.29
```

# Installing genome sequences and annotations from Ensembl genomes

For this to work we first need to get the current list of supported genomes in the respective GROUP,
'Plants' in the example. This can be done by running the following command, where X is the desired release number:

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X organisms
```

Once this is done then you can install all genomes with one command.

```{r, engine='bash', eval=FALSE}
cd $RSAT

# download FASTA and GTF files
nohup make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X download_all_species

# parse input files, extract genomic features and compute oligo frequencies
nohup make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X install_all_species

# check upstream sequences can be retrieved
nohup make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X check_all_species

# compute descriptive genome stats and produce HTML report
# such as http://rsat.eead.csic.es/plants/data/stats/
make -f makefiles/ensemblgenomes_FTP_client.mk calc_stats
```

Note that this can take a long time, days for current Ensembl Plant releases.

In case you want to install a single genome you can do that by:

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X SPECIES=oryza_longistaminata download_one_species
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X SPECIES=oryza_longistaminata install_one_species
```

If variation data is available for your species of interest at Ensembl Genomes you can download it with:

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=X SPECIES=oryza_sativa variations_one_species
```


You can even install selected genomes from older releases (Y in the example):

```{r, engine='bash', eval=FALSE}
cd $RSAT
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=Y organisms
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=Y SPECIES=oryza_sativa download_one_species
make -f makefiles/ensemblgenomes_FTP_client.mk GROUP=Plants RELEASE=Y SPECIES=oryza_sativa install_one_species
```

# Installing from other sources

The installation procedures described can also be used to install arbitrary genomes from other sources,
provided that 4 input files are available with the following extensions: 

* SPECIES_RSAT_ID.dna.toplevel.fa : raw genomic sequence
* SPECIES_RSAT_ID.dna_rm.genome.fa : repeat-hard-masked genomic sequence 
* SPECIES_RSAT_ID.gtf : annotation file 
* SPECIES_RSAT_ID.pep.all.fa : peptide sequences of CDS features

where SPECIES_RSAT_ID is a string identifying this organism and its annotation in RSAT. 
For instance, for assembly Wm82.a2.v1 of *Glycine max* from <JGI>, we could install with:
 
```{r, engine='bash', eval=FALSE}
cd $RSAT

mkdir -p /var/www/html/rsat/data/genomes/Glycine_max.Wm82.a2.v1.JGI/genome
# put there those 4 files (dna.toplevel.fa,dna_rm.genome.fa,.gtf,.pep.all.fa)

make -f makefiles/ensemblgenomes_FTP_client.mk SPECIES=Glycine_max \
    SPECIES_DIR=/var/www/html/rsat/data/genomes/Glycine_max.Wm82.a2.v1.JGI \
    SPECIES_RSAT_ID=Glycine_max.Wm82.a2.v1.JGI TAXON_ID=3847 GTF_SOURCE=JGI \
    install_from_gtf
```

Note that TAXON_ID can be obtained at https://www.ncbi.nlm.nih.gov/taxonomy
