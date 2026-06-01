#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Install some python libraries with pip
##
## Note: numpy, scipy and matplotlib are supposed to have previously
## been installed with apt-get under Ubuntu. For other OS, they should
## be added to the pip installation.

echo
echo "Installing Python packages for RSAT"
echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

sudo apt-get update

    PACKAGES_REQUIRED="
python3-numpy
python3-scipy
python3-matplotlib 
snakemake 
python3-rpy2 
python3-pysimplesoap 
python3-requests 
python3-yaml
python3-suds
python3-venv
python3-pip
python3-pygraphviz
python3-requests
"

sudo apt-get install -y ${PACKAGES_REQUIRED}

# python3-fisher ## 2025-12-29: unable to install

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
