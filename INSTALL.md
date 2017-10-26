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

- the latest stable release named `rsat_20XX-XX-XX.tar.gz` (XX-XX-XX must be replaced by the release date).
  
- the previous release archives available in the folder `previous_versions`.
	 
## Downloading the latest RSAT release via the Web installer  

5. Download the tar archive named  

`rsat_20XX-XX-XX.tar.gz`  

where `XX-XX-XX` is the latest release date and put it in your chosen directory.  

6. Uncompress the archive. This will create a directory named `rsat` where we will continue the install procedure.

```
## By default we install the package in the tool directory.
## This should be adapted according to your local setup. 
export INSTALL_ROOT=/packages
sudo mkdir -p ${INSTALL_ROOT}

## Replace XX-XX-XX by the actual release
tar -C ${INSTALL_ROOT}/ -xpvzf rsat_20XX-XX-XX.tar.gz
cd ${INSTALL_ROOT}/rsat
```

## Configuring your local RSAT instance

7. We will start by auto-configuring RSAT in order to deine suitable basic parameters. 


For VirtualBox instances, replace [your.server.IP] by the actual IP address of your server. 

```
## Get your IP address and check if it is 192.168.56.101
ifconfig | grep inet

## Semi-auto configuration for VirtualBox VM
## (adapt IP address if required)
perl perl-scripts/configure_rsat.pl -auto  \
  rsat_site=rsat-vb-2017-10 \
  rsat_www=http://192.168.56.101/rsat/ \
  rsat_ws=http://192.168.56.101/rsat/ \
  ucsc_tools=1 \
  ensembl_tools=1
```

For the IFB cloud (IP address will change at each instance)

```
perl perl-scripts/configure_rsat.pl -auto  \
  rsat_site=rsatvm-ifb-2017-10 \
  RSAT=${INSTALL_ROOT}/rsat \
  rsat_www=auto \
  rsat_ws=auto \
  phylo_tools=0 \
  compara_tools=0 \
  variations_tools=0 \
  ucsc_tools=0 \
  ensembl_tools=0 \
  SUDO=sudo
  
chmod 755 /root # required for apache user to access the packages
```

8. We can now refine the configuration by choosing custom parameter to your RSAT instance (for example the email of the local admin, the instance name, ...).

```
perl perl-scripts/configure_rsat.pl
```


## Installing RSAT

Before running the installation, it might be worth updating the Linux distribution (`apt-get update`) in order to get the latest versions of the basic packages. 

```
## This requires admin privileges
sudo bash

## Check who you are  (should be root)
whoami


## Go to the RSAT directory
export INSTALL_ROOT=/packages
cd ${INSTALL_ROOT}/rsat

## Read config and run bash installation scripts
source RSAT_config.bashrc && \ 
bash installer/01_ubuntu_packages.bash && \
bash installer/02_python_packages.bash  && \
bash installer/03_install_rsat.bash && \
bash installer/04_perl_packages.bash  && \
bash installer/06_install_organisms.bash && \
bash installer/07_R-and-packages.bash  && \
bash installer/08_apache_config.bash && \
bash installer/09_rsat_ws.bash && \
bash installer/10_clean_unnecessary_files.bash

## Restore rsat as owner of the $RSAT folder
chown -R rsat.rsat $RSAT

## Exit sudo session
exit

## Check who you are (should be back to normal user identity)
whoami
```

## Testing the command lines

```
make -f makefiles/install_tests.mk all
```

This makefile runs a series of tests for different components of the *RSAT* suite. Each test result is stored in a separate file in the test directory (`./install_tests` by default). Output file names are printed out after each test. 


## Testing the Web server

TO BE WRITTEN


****************************************************************
# Supplementary information

## Files

rsat_YYYY-MM-DD.tar.gz

   Stand-alone versions + web servers for RSAT and NeAT (except the
   metabolic tools).

metabolic-tools_YYYYMMDD.tar.gz

   Metabolic pathway analysis tools (supported on some NeAT servers).

## RSAT/NeAT installation and user guides

After having uncompressed the archive, you will find the installation
and user guides in the `doc/manuals` directory

```
ls -1 public_html/release/*.pdf
```


| Guide | File |
|------------------------|---------------------------|
| RSAT installation guide |   RSAT_install_guide.pdf |
| RSAT Web configuration guide |   rsat_web_server.pdf |
| RSAT Command-linde user guide |  tutorial_shell_rsat.pdf |
| NeAT Web server configuration |  neat_web_server.pdf |
| NeAT Command-line user guide |   neat_tutorial.pdf |

****************************************************************

