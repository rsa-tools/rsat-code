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

apt install -y python3-numpy
apt install -y python3-scipy
apt install -y python3-matplotlib 
apt install -y snakemake 
apt install -y python3-rpy2 
apt install -y python3-pysimplesoap 
apt install -y python3-requests 
apt install -y python3-yaml
#apt install -y python3-wsdl
pip3 install wstools
apt install -y python3-suds
apt install -y python3-venv
apt install -y python3-pip
# apt install -ypython3-fisher ## 2025-12-29: unable to install
apt install -y python3-pygraphviz
# repeated # apt install -y python3-requests   # required for downloading GO terms from biomart

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
