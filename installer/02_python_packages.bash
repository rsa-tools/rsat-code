################################################################
## Install some python libraries with pip
##
## Note: numpy, scipy and matplotlib are supposed to have previously
## been installed with apt-get under Ubuntu. For other OS, they should
## be added to the pip installation.

# source installer/00_config.bash

source ${RSAT}/RSAT_config.bashrc

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
sudo pip3 install rpy2  ## THIS FAILS on the IFB cloud. To be checked.
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages
sudo pip3 install suds-jurko
sudo pip3 install pysimplesoap
sudo pip3 install soappy

## Check disk usage
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
