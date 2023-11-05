#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Install all the Ubuntu packages required prior to the installation
## of the Regulatory Sequence Analysis Tools (RSAT; http://rsat.eu/).

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo


## PRIOR REQUIREMENT: Mac developer tools and brew (for make and other utilities)

################################################################
## Required brew packages
################################################################
## Required apt-get packages
PACKAGES_REQUIRED="
make
java
cvs
wget
zip
unzip
screen
libgd
apache2
ghostscript
gnuplot
graphviz
openssl
mysql-connector-c
mysql
python3
emacs
ntp
curl
php
"

## 2023-11-05: JvH suppresses this because python2 is not suppported anymore
# python


################################################################
## This variable was left for consistency with the corresponding
## ubuntu script
PACKAGES_PERL=""


## Install the brew libraries
PACKAGES="${PACKAGES_REQUIRED} ${PACKAGES_PERL}"
echo "Packages to be installed with ${OS_INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES}"
for LIB in ${PACKAGES}; do \
    echo "`date '+%Y/%m/%d %H:%M:%S'` installing brew library ${LIB}" ; \
    echo "${OS_INSTALLER} install ${INSTALLER_OPT} ${LIB}" ; \
    ${OS_INSTALLER} install ${INSTALLER_OPT} ${LIB} > ${RSAT}/install_logs/install_${LIB}_log.txt ; \
    df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${LIB}_installed.txt ; \
done
echo "Log files are in folder ${RSAT}/install_logs"
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt

## DONE: installation of Mac OS X packages via brew
################################################################
