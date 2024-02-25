#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Install python libraries
##

# ## optional: an utility to measure internet bandwidth
# sudo pip2 install speedtest-cli

#${OS_INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-suds

echo "# Installing Python libraries for MAC OS with brew"; 
brew install numpy
brew install scipy
brew install python-matplotlib
brew install fisher
brew install snakemake
brew install python-requests
brew install PyYAML

# may2023, ubuntu 20.04
# to avoid error with rpy2-3.5.12.tar.gz "note: This error originates from a subprocess, and is likely not a problem with pip",
# not sure this solves "THIS FAILS on the IFB cloud. To be checked"

## 2024-02-24 : I don't find these libraries wih brew for MAC OS
# brew install python-wheel
# brew install rpy2 
# sudo pip3 install suds ## Error on Mac OSX
# brew install suds
# sudo pip3 install pysimplesoap


## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages
# sudo pip3 install suds-jurko ## Inactivated on 2023-02-06 because does not work anymore with Ubuntu 22.04
# soappy is not maintained, and don't install well on python3
# sudo pip3 install soappy

# required for downloading GO terms from biomart

# overide old pyyaml on python2

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
