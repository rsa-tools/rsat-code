################################################################
## Install all the Ubuntu packages required prior to the installation
## of the Regulatory Sequence Analysis Tools (RSAT; http://rsat.eu/).

source installer/00_config.bash

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

## Set time zone in non-interactive mode
#echo Europe/Rome > /etc/timezone
# sudo dpkg-reconfigure -f noninteractive tzdata

## If not, set up the time zone, date and time with this command
## (source: https://help.ubuntu.com/community/UbuntuTime).
# dpkg-reconfigure tzdata



## We need to update apt-get, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu

## We can then check the increase of disk usage during the different
## steps of the installation
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt

################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Still required for Ubuntu 16.04 ? TO BE CHECKED
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Fix a problem with rabbitmq in the Ubuntu 14.04 distrib
## wget -O- https://www.rabbitmq.com/rabbitmq-release-signing-key.asc | apt-key add -

## Install aptitude, more efficient than apt-get to treat dependencies
## when installing and uninstalling packages.
## TO SAVE SPACE, I SUPPRESS aptitude
## apt-get install aptitude
apt-get update
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_apt-get_updated.txt

################################################################
## I tried to run dist-upgrade because it 'can "intelligently" handle
## changes in the dependencies system. This includes removing packages
## that are no longer necessary or resolve conflicts between packages
## that arose because of changes in the dependencies.'
## http://askubuntu.com/questions/194651/why-use-apt-get-upgrade-instead-of-apt-get-dist-upgrade
## HOWEVER, THE PERL UPDATE DOES NOT WORK ANYMORE AFTER THAT !!!
${OS_INSTALLER} ${INSTALLER_OPT} upgrade
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${OS_INSTALLER}_upgrade.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt


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
mysql-client
default-jre
python
python-pip
python-setuptools 
python-numpy
python-scipy
python-matplotlib
python-suds
python-rpy2
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
perl-doc
pmtools
libyaml-perl
libemail-simple-perl
libemail-sender-perl
libemail-simple-creator-perl
libpostscript-simple-perl
libstatistics-distributions-perl
libio-all-perl
libobject-insideout-perl
libobject-insideout-perl
libsoap-lite-perl
libsoap-wsdl-perl
libxml-perl
libxml-simple-perl
libxml-compile-cache-perl
libdbi-perl
liblockfile-simple-perl
libobject-insideout-perl
libgd-perl
libdbd-mysql-perl
libjson-perl
libbio-perl-perl
libdigest-md5-file-perl
libnet-address-ip-local-perl
libemail-sender-transport-smtp-tls-perl
"

## We did not find apt-get packages for some required Perl
## libraries. These will have to be installed with cpan.
PACKAGES_PERL_MISSING="
libalgorithm-cluster-perl
digest-md5-file-perl
liblockfile-simple
libutil-properties-perl
librest-client-perl
libxml-compile-soap11-perl
libxml-compile-wsdl11-perl
libxml-compile-transport-soaphttp-perl
libbio-das-perl        
"

## Install the apt-get libraries
echo "Packages to be installed with ${OS_INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES}"
echo "Perl module packages to be installed with ${OS_INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES_PERL}"
for LIB in ${PACKAGES} ${PACKAGES_PERL}; \
do \
   echo "`date '+%Y/%m/%d %H:%M:%S'`  installing apt-get library ${LIB}" ; \
   ${OS_INSTALLER} install ${INSTALLER_OPT} ${LIB} > ${RSAT}/install_logs/${OS_INSTALLER}_install_${LIB}.txt ; \
   df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${LIB}_installed.txt ; \
done
echo "Log files are in folder ${RSAT}/install_logs"
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt


## PROBLEMS WITH Ubuntu 16
#
#apt-get install php5 ## E: Package 'php5' has no installation candidate
#apt-get install libapache2-mod-php5 ## E: Package 'libapache2-mod-php5' has no installation candidate
# Fix: php (7) is part of the Ubuntu distribution
# See: http://askubuntu.com/questions/756879/cant-install-php5-on-ubuntu-16-04
# I NEED TO CHECK IF THE PHP INTERFACES STILL WORK WITH PHP7.
#
# 2016/10/24 07:23:02  installing apt-get library python3-rpy2
# E: Unable to correct problems, you have held broken packages.


## This package has to be installed in an interactive mode (dialog
## box)
#${OS_INSTALLER} install ${INSTALLER_OPT} console-data

################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## DOES NOT WORK ON Ubunutu 16.04.1 anmore
## Not sure it is still required though.
## TO BE CHECKED
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ################################################################
# ## Specific treatment for some Python libraries
# ##
# ## A fix for a problem to install scipy with pip: use ${OS_INSTALLER} build-dep 
# ## taken from here: http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
# ## Note that these dependencies cost 400Mb ! To be checked
# ${OS_INSTALLER} ${INSTALLER_OPT} build-dep python-numpy python-scipy
# df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_numpy-scipy_dependencies_installed.txt
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt

################################################################
## To free space, remove apt-get packages that are no longer required.a
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt
${OS_INSTALLER} ${INSTALLER_OPT}  autoremove
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_autoremoved.txt
${OS_INSTALLER} ${INSTALLER_OPT}  clean
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned.txt
## This really helps: it saves several hundreds Mb
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt

## DONE: installation of Ubuntu packages
################################################################

