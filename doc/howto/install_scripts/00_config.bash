################################################################
## Configuration for the installation of the Regulatory Sequence
## Analysis Tools (RSAT; http://rsat.eu/) on an Ubuntu Linux system.
##

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!   NOTE ABOUT KASPERSKY ANTIVIRUS   !!!!!!!!!!!!!!
## 
## For the installation of Perl package and for third-party Linux
## packages, I need to temporarily inactivate the antivirus software
## Kaspersky.
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Note: for Debian, I must set the following locales manually so I
## show them commented.
#
# export LANGUAGE=en_US.UTF-8
# export LANG=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# locale-gen en_US.UTF-8
# dpkg-reconfigure locales

## Note: you should probably replace the automatic site (hostname) by
## a name of your choice.
export RSAT_SITE=`hostname`

## Path for the local installation
export RSAT_PARENT_PATH=/packages
export RSAT=${RSAT_PARENT_PATH}/rsat

## URL to download the RSAT distribution
export RSAT_RELEASE=2017-02-14 ## Version to be downloaded from the tar distribution
export RSAT_ARCHIVE=rsat_${RSAT_RELEASE}.tar.gz
export RSAT_DISTRIB_URL=http://pedagogix-tagc.univ-mrs.fr/download_rsat/${RSAT_ARCHIVE}

## Configuration for the installation
export INSTALLER=apt-get
export INSTALLER_OPT="--quiet --assume-yes"
## alternative: INSTALLER=aptitude


## Create a separate directory for RSAT, which must be readable by all
## users (in particular by the apache user)
echo "Creating RSAT_PARENT_PATH ${RSAT_PARENT_PATH}"
mkdir -p ${RSAT_PARENT_PATH}
cd ${RSAT_PARENT_PATH}
mkdir -p ${RSAT_PARENT_PATH}/install_logs
chmod 777 ${RSAT_PARENT_PATH}/install_logs
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_start.txt


## Check the installation device 
DEVICE=`df -h | grep '\/$' | perl -pe 's/\/dev\///' | awk '{print $1}'`
echo "Installation device: ${DEVICE}"
## This should give something like sda1 or vda1. If not check the device with df
