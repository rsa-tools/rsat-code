#!/usr/bin/env bash

source $(dirname $0)/00_config.bash
OS_INSTALLER=brew
INSTALLER_OPT=

################################################################
## Install all the Mac OS packages required prior to the installation
## of the Regulatory Sequence Analysis Tools (RSAT; http://rsat.eu/).

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

## PRIOR REQUIREMENT: Mac developer tools and brew (for make and other utilities)

################################################################
## Required brew packages
PACKAGES_REQUIRED="
make
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
emacs
curl
php
pkg-config
python3
bedtools
blast
diamond
"

PACKAGED_NOT_REQUIRED="
mcl
ntp
java
mysql-connector-c
mysql
"

################################################################
## This variable was left for consistency with the corresponding
## ubuntu script, but actually on MAC OS we install all Perl packages
## with cpan.
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
