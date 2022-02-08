#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

################################################################
## Download RSAT distribution from the tar.gz archive at
## http://download.rsat.eu/

RSAT_DOWNLOAD_SITE=http://pedagogix-tagc.univ-mrs.fr/download_rsat/
RSAT_ARCHIVE=`wget -q -O- ${RSAT_DOWNLOAD_SITE} | grep rsat | grep -v scripts | grep class | awk -F'"' '{print  $6}'`
RSAT_VERSION=`echo $RSAT_ARCHIVE | sed 's/rsat_//g' | sed 's/.tar.gz//g'`
RSAT_DISTRIB_URL=${RSAT_DOWNLOAD_SITE}/${RSAT_ARCHIVE}

echo "Downloading RSAT archive $RSAT_ARCHIVE"
echo "   URL: ${RSAT_DISTRIB_URL}"
echo ${RSAT_VERSION} > ${RSAT_PARENT_PATH}/RSAT_version
cd ${RSAT_PARENT_PATH}
wget --no-clobber ${RSAT_DISTRIB_URL}
tar -xpzf ${RSAT_ARCHIVE}
cd ${RSAT}


## Note: the git distribution requires an account at the ENS git
## server, which is currently only possible for RSAT developing team.
## In the near future, we may use git also for the end-user
## distribution.


## Alternative: get RSAT from the git repository.
## THIS IS RESERVED TO RSAT DEVELOPING TEAM
# cd ${RSAT_PARENT_PATH}
# git config --global user.mail rsat@rsat-vm-${RSAT_RELEASE}
# git config --global user.name "rsat"
# git config --global core.editor emacs
# git config --global merge.tools meld
# git config --list
# git clone git@depot.biologie.ens.fr:rsat


## For users who don't have an account on the RSAT git server, the
## code can be downloaded as a tar archive from the Web site.

## Make a link from home directory to find RSAT home
ln -fs ${RSAT} ${HOME}/rsat

## Metabolic pathway tools installation
##
## TO BE DONE LATER
## wget http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/metabolic-tools_20110408.tar.gz

# ## DONE : obtained RSAT distribution
# ################################################################
