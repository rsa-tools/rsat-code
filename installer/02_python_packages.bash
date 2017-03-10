################################################################
## Install some python libraries with pip
##
## Note: numpy, scipy and matplotlib are supposed to have previously
## been installed with apt-get under Ubuntu. For other OS, they should
## be added to the pip installation.

source 00_config.bash

pip install soappy
pip install fisher
pip install httplib2
## pip install pygraphviz ## OSError: Error locating graphviz.

## optional: an utility to measure internet bandwidth
pip install speedtest-cli

#${OS_INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-suds

## Failures: no distributions at all found
# pip3 install wsdl
# pip3 install wstools
pip3 install fisher
pip3 install snakemake
pip3 install rpy2  ## THIS FAILS on the IFB cloud. To be checked.
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages
pip3 install suds-jurko
pip3 install pysimplesoap
pip3 install soappy

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
