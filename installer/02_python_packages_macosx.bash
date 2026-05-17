#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Install python libraries
##

# ## optional: an utility to measure internet bandwidth
# sudo pip2 install speedtest-cli

#${OS_INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-sudsecho "# Installing Python libraries for MAC OS with brew"; 
brew install numpy
brew install scipy
brew install python-matplotlib
brew install fisher
brew install snakemake
brew install python-requests
brew install PyYAML
brew install pipx


## 2024-02-24 : I don't find these libraries wih brew for MAC OS
# brew install python-wheel
# brew install rpy2 
# sudo pip3 install suds ## Error on Mac OSX
# brew install suds
# sudo pip3 install pysimplesoap


## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

# Install weblogo system-wide so that it is accessible
# to the Apache user running RSAT services.
echo "Installing weblogo via python2 pip requires sudo rights"
sudo -H python3 -m pip install weblogo

# required for downloading GO terms from biomart

# overide old pyyaml on python2

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
