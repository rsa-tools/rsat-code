# Regulatory Sequence Analysis Tools (RSAT)

This folder contains the code required to run a local version of the
software suite **Regulatory Sequence Analysis Tools** (**RSAT**).


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

- [**Installing RSAT**](https://rsa-tools.github.io/installing-RSAT): on [Unix systems](https://rsa-tools.github.io/installing-RSAT/unix-install-rsat/installing_RSAT_procedure.html), [Docker container](https://rsa-tools.github.io/installing-RSAT/unix-install-rsat/installing_RSAT_procedure.html)

- [**Managing RSAT**](https://rsa-tools.github.io/managing-RSAT)

- [**Motif databases**](https://github.com/rsa-tools/motif_databases)
