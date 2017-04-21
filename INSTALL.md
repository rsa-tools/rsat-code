---
title: "RSAT installation guide"
author: "Jacques van Helden"
date: 'Last update: 2017/04/21'
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

This document explains how to install and configure the Regulatory Sequence Analysis Tools (RSAT). 

# Downloading RSAT

1. Go to the RSAT portal (**<http:rsat.eu/>**)
2. Select any server. 
3. In the left-sided menu,  click on the **Download** link. 
4. Fill in your name and coordinates and accept the license.


You now have access to the download page containing the links to  
     * the latest stable release named  
     `rsat_20XX-XX-XX.tar.gz`  
     (XX-XX-XX are replaced by the release date).
  
     * the previous release archives available in the folder  
       `previous_versions`.
	 
## Downloading the latest RSAT release via the Web installer  

5. Download the tar archive named  

`rsat_20XX-XX-XX.tar.gz`  

where `XX-XX-XX` is the latest release date and put it in your chosen directory.  

6. Uncompress the archive. This will create a directory named `rsat` where we will continue the install procedure.

```
sudo mkdir -p /packages
sudo chown rsat.rsat /packages
cd /packages

## Replace XX-XX-XX by the actual release
tar -xpvzf /home/rsat/Downloads/rsat_20XX-XX-XX.tar.gz
cd rsat
```

## Configuring your local RSAT instance

7. We will start by auto-configuring RSAT in order to deine suitable basic parameters. For this, in the `rsat` directory, type:

```
perl perl-scripts/configure_rsat.pl -auto  \
  rsat_site=rsat-vb-2017-04 \
  rsat_www=http://192.168.56.101/rsat/ \
  rsat_ws=http://192.168.56.101/rsat/ \
  ucsc_tools=1 \
  ensembl_tools=1
```

8. We can now refine the configuration by choosing custom parameter to your RSAT instance (for example the email of the local admin, the instance name, ...).

```
perl perl-scripts/configure_rsat.pl
```


## Installing RSAT

```
sudo bash
cd /packages/rsat
source RSAT_config.bashrc
bash installer/01_ubuntu16.4_packages.bash
bash installer/02_python_packages.bash 
bash installer/03_install_rsat.bash
bash installer/04_perl_packages.bash 
bash installer/06_install_organisms.bash
bash installer/07_R-and-packages.bash 
bash installer/08_apache_config.bash 
bash installer/09_rsat_ws.bash 
bash installer/10_clean_unnecessary_files.bash
```

## Testing the command lines

```
make -f makefiles/install_tests.mk all
```

This makefile runs a series of tests for different components of the *RSAT* suite. Each test result is stored in a separate file in the test directory (`./install_tests` by default). Output file names are printed out after each test. 


## Testing the Web server

****************************************************************
# Supplementary information

## Files

rsat_YYYY-MM-DD.tar.gz

   Stand-alone versions + web servers for RSAT and NeAT (except the
   metabolic tools).

metabolic-tools_YYYYMMDD.tar.gz

   Metabolic pathway analysis tools (supported on some NeAT servers).

## RSAT/NeAT installation

After having uncompressed the archive, you will find the installation
and user guides in the directory

      rsa-tools/doc/manuals/*.pdf

## Regulatory Sequence Analysis Tools (RSAT)

RSAT installation guide:   RSAT_install_guide.pdf
Web configuration guide:   rsat_web_server.pdf
Command-linde user guide:  tutorial_shell_rsat.pdf

## Network Analysis Tools (NeAT)

Web server configuration:  neat_web_server.pdf
Command-line user guide:   neat_tutorial.pdf

****************************************************************