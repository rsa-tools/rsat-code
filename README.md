# Regulatory Sequence Analysis Tools (RSAT)

[![Biocontainers](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/biocontainers/rsat/tags)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/rsat-core/badges/version.svg)](https://anaconda.org/bioconda/rsat-core)

This folder contains the code required to run an instance of the 
software suite **Regulatory Sequence Analysis Tools** (**RSAT**).

Releases are used to mark development progress and might trigger
updates [Docker containers](https://hub.docker.com/r/biocontainers/rsat/tags) and 
[conda recipes](https://anaconda.org/bioconda/rsat-core). 
Release names have the format 2025.04.04 so that they work also in conda.

**Contact**: Jacques.van-Helden[AT]univ-amu.fr , rsat-contact[AT]list01.bio.ens.psl.eu

**Authors**: Jacques van Helden + the RSAT team.


****************************************************************
# Files and folders

## Folders coming with the release

- **contrib** software tools developed by collaborators, in C,
  requiring compilation. 
  
- **doc** documentation in various formats.

- **perl-scripts** Perl scripts, which constitute most of the RSAT and
NeAT code.

- **public_html** Web interface (cgi + php scripts).

- **python-scripts** Python scripts

- **ws_clients** examples of Web services clients in different
  languages.

## Folders created for management

- **app_sources** download the code source of diverse applications
  (dependencies plus optional complementary software tools).
  
- **bin** executable binaries (including the compiled C programs).

- **downloads** this folder is created at RSAT initialisation. It
  serves to download the genome data from NCBI, Ensembl or other
  databases, which will serve to install organisms on RSAT.

****************************************************************

# Documentation

- [**Installing RSAT**](https://rsa-tools.github.io/installing-RSAT): 
	+ on [Unix systems](https://rsa-tools.github.io/installing-RSAT/unix-install-rsat/installing_RSAT_procedure.html)
	+ with [conda](https://rsa-tools.github.io/installing-RSAT/conda-install-rsat/bioconda-rsat-core.html)
	+ as container with motif collections: [Docker](https://rsa-tools.github.io/installing-RSAT/RSAT-Docker/RSAT-Docker-tuto.html), [Apptainer](https://rsa-tools.github.io/installing-RSAT/RSAT-Docker/RSAT-Apptainer-tuto.html)

- [**Managing RSAT**](https://rsa-tools.github.io/managing-RSAT)

- [**Motif databases**](https://github.com/rsa-tools/motif_databases)
