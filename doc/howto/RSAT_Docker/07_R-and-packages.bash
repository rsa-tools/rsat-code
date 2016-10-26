source 00_config.bash

################################################################
##
## Install R
##
## This is done after the other packages becayse I need to declare
## R-cran as source in order to install the latest version of R (3.3.1
## on 2016-10) which is required for some R scripts, but not
## distributed with Ubuntu 14.04 (this Ubuntu release 14.04 comes with
## R version 3.0.2).
##
## I add the row before the original sources.list because there is
## some problem at the end of the update.
grep -v cran.rstudio.com /etc/apt/sources.list > /etc/apt/sources.list.bk
echo "## R-CRAN repository, to install the most recent version of R" > /etc/apt/sources.list.rcran
echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list.rcran
echo "" >> /etc/apt/sources.list.rcran
cat /etc/apt/sources.list.rcran   /etc/apt/sources.list.bk >  /etc/apt/sources.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
add-apt-repository ppa:marutter/rdev
apt-get update
${INSTALLER} install ${INSTALLER_OPT} r-base > ${RSAT_PARENT_PATH}/install_logs/${INSTALLER}_install_r-base.txt

