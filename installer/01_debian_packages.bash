#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Install all the Debian packages required prior to the installation
## of the Regulatory Sequence Analysis Tools (RSAT; http://rsat.eu/).
##
## Note: as of Debian 11 bullseye, there are packages that require enabling the non-free repositores
NONFREE_PACKAGES="libxml-compile-soap-perl"

echo
echo "Debian packages that require enabling non-free repos: ${NONFREE_PACKAGES}"
echo

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

################################################################
## Before anything else, check that the date, time and time zone are
## correctly specified
date

${OS_INSTALLER} install -y openssh-client

# manage local time
if [ ! -f /etc/localtime ]
then
    ${OS_INSTALLER} install ${INSTALLER_OPT} tzdata
    sudo ln -sf /usr/share/zoneinfo/Europe/Paris /etc/localtime
    sudo dpkg-reconfigure --frontend noninteractive tzdata
fi

################################################################
## Required apt-get packages
PACKAGES_REQUIRED="
apt-utils
make
net-tools
ssh
git
cvs
wget
zip
unzip
screen
g++
apache2
libgdbm-dev
libgd-tools
libgd-gd2-perl
ghostscript
gnuplot
graphviz
default-mysql-client
default-jre
python
python-setuptools
python3
python3-pip
python3-setuptools
python3-numpy
python3-scipy
python3-matplotlib
python3-rpy2
emacs
x11-apps
eog
ntp
curl
libcurl4-openssl-dev
libcurl4-gnutls-dev
libxml2-dev
libnet-ssleay-perl
libcrypt-ssleay-perl
libssl-dev
php
libapache2-mod-php
libapache2-mod-wsgi-py3
rsync
"


################################################################
## Packages to be checked by JvH.
## These are useful to me, but I am not sure they are required for RSAT.
PACKAGES_OPT="
ess
yum
php-elisp
libgd2-xpm-dev
libxml2-dev
links
gfortran
libmysqlclient-dev
texlive-latex-base
python-virtualenv
ipython
ipython-notebook
libreadline-gplv2-dev:i386
lib64readline-gplv2-dev:i386
libreadline-gplv2-dev
libx11-dev
libxt-dev
libxml2-dev
tcl8.5-dev
tk8.5-dev
libxss-dev
libpng12-dev
libjpeg62-dev
libcairo2-dev
lib32z1
lib32ncurses5
lib32bz2-1.0
libc6-dev
build-essential
python-dev
python3-dev
libnet-ssleay-perl
libcrypt-ssleay-perl
exfat-fuse
exfat-utils
at
firefox
ncbi-blast+
finger
"

################################################################
## apt-get packages to install Perl modules (not properly speaking
## necessary, could be done with cpan, but ensure consistency with
## ubuntu OS)
PACKAGES_PERL="
bioperl-run
libbio-das-lite-perl
libbio-perl-perl
libclass-std-perl
libdbd-mysql-perl
libdbi-perl
libdigest-md5-file-perl
libemail-sender-perl
libemail-simple-creator-perl
libemail-simple-perl
libgd-perl
libio-all-perl
libjson-perl
liblockfile-simple-perl
liblog-log4perl-perl
libnet-address-ip-local-perl
libnet-smtp-tls-perl
libnet-smtps-perl
libnumber-format-perl
libobject-insideout-perl
libole-storage-lite-perl
libparallel-forkmanager-perl
libpostscript-simple-perl
librest-client-perl
libsoap-lite-perl
libsoap-wsdl-perl
libspreadsheet-xlsx-perl
libstatistics-distributions-perl
libxml-compile-cache-perl
libxml-compile-soap-perl
libxml-compile-wsdl11-perl
libxml-parser-perl
libxml-perl
libxml-simple-perl
libyaml-perl
perl-doc
pmtools
"


## We did not find apt-get packages for some required Perl
## libraries. These will have to be installed with cpan.
PACKAGES_PERL_MISSING="
libalgorithm-cluster-perl
libutil-properties-perl
"


## Install the apt-get libraries
PACKAGES="${PACKAGES_REQUIRED} ${PACKAGES_PERL}"
echo "Packages to be installed with ${OS_INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES}"
for LIB in ${PACKAGES}; do \
    echo "`date '+%Y/%m/%d %H:%M:%S'`  installing apt-get library ${LIB}" ; \
    ${OS_INSTALLER} install ${INSTALLER_OPT} ${LIB} > ${RSAT}/install_logs/install_${LIB}_log.txt ; \
    df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${LIB}_installed.txt ; \
done
echo "Log files are in folder ${RSAT}/install_logs"

################################################################
## To free space, remove apt-get packages that are no longer required.a
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
${OS_INSTALLER} ${INSTALLER_OPT}  autoremove
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_autoremoved.txt
${OS_INSTALLER} ${INSTALLER_OPT}  clean
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned.txt
## This really helps: it saves several hundreds Mb
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt

## DONE: installation of Debian packages
################################################################
