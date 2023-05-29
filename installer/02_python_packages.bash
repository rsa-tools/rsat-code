#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Install some python libraries with pip
##
## Note: numpy, scipy and matplotlib are supposed to have previously
## been installed with apt-get under Ubuntu. For other OS, they should
## be added to the pip installation.

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

## pip2 is no longer supported
# sudo pip2 install numpy
# sudo pip2 install scipy
# sudo pip2 install matplotlib
# sudo pip2 install suds
# sudo pip2 install soappy
# sudo pip2 install fisher
# sudo pip2 install httplib2

# ## Some python modules required for the REST web services
# sudo pip2 install flask flask-restplus requests werkzeug jinja2 click itsdangerous
# ## pip install pygraphviz ## OSError: Error locating graphviz.

# ## optional: an utility to measure internet bandwidth
# sudo pip2 install speedtest-cli

#${OS_INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-suds

## Failures: no distributions at all found
# pip3 install wsdl
# pip3 install wstools
sudo pip3 install numpy
sudo pip3 install scipy
sudo pip3 install matplotlib
# sudo pip3 install suds ## Error on Mac OSX
sudo pip3 install fisher
sudo pip3 install snakemake

# may2023, ubuntu 20.04
# to avoid error with rpy2-3.5.12.tar.gz "note: This error originates from a subprocess, and is likely not a problem with pip",
# not sure this solves "THIS FAILS on the IFB cloud. To be checked"
sudo pip install wheel setuptools pip --upgrade
sudo pip3 install wheel setuptools pip --upgrade
sudo pip3 install rpy2 

## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages
sudo pip3 install pysimplesoap
# sudo pip3 install suds-jurko ## Inactivated on 2023-02-06 because does not work anymore with Ubuntu 22.04
# soappy is not maintained, and don't install well on python3
# sudo pip3 install soappy

# required for downloading GO terms from biomart
sudo pip3 install requests

# overide old pyyaml on python2
pip3 install --ignore-installed PyYAML

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
