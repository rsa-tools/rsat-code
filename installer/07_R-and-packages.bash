source installer/00_config.bash


echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables

################################################################
##
## Install R
##
## This is done after the other packages because I need to declare
## R-cran as source in order to install the latest version of R (3.3.1
## on 2016-10) which is required for some R scripts, but not
## distributed with Ubuntu 14.04 (this Ubuntu release 14.04 comes with
## R version 3.0.2).
##
## I add the row before the original sources.list because there is
## some problem at the end of the update.
export CODENAME=`grep CODENAME /etc/lsb-release | cut -c 18-`
export CRAN_DIR=${CODENAME}-cran35
#export CRAN_URL=https://pbil.univ-lyon1.fr/CRAN
export CRAN_URL=http://cran.rstudio.com
# grep -v -i cran /etc/apt/sources.list > /etc/apt/sources.list.bk
# sudo echo "## R-CRAN repository, to install the most recent version of R" > /etc/apt/sources.list.rcran
# sudo echo "deb ${CRAN_URL}/bin/linux/ubuntu ${CRAN_DIR}/" >> /etc/apt/sources.list.rcran
# sudo echo "" >> /etc/apt/sources.list.rcran
# sudo cat /etc/apt/sources.list.rcran   /etc/apt/sources.list.bk >  /etc/apt/sources.list
sudo apt install apt-transport-https software-properties-common
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu ${CRAN_DIR}/"
sudo grep -i "cran" /etc/apt/sources.list # Check the cran repository


# Get code to allow using cran as apt-get source
#sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 
## Not working on the IFB cloud, so I comment it and use the option --force-yes

################################################################
## NOTE: THIS IS FOR UBUNTU ONLY ###

# sudo ${OS_INSTALLER} install ${INSTALLER_OPT} software-properties-common python-software-properties
# sudo add-apt-repository ppa:marutter/rdev
sudo apt-get update
# INSTALLER_OPT=--force-yes
export INSTALLER_OPT="--allow-unauthenticated --allow-downgrades --allow-remove-essential --allow-change-held-packages"
sudo ${OS_INSTALLER} install ${INSTALLER_OPT} r-base:all/${CRAN_DIR} # > ${RSAT}/install_logs/install_r-base_log.txt 
sudo ${OS_INSTALLER} install ${INSTALLER_OPT} r-base-dev:all/${CRAN_DIR} # > ${RSAT}/install_logs/install_r-base-dev_log.txt 



################################################################
## Install selected R librairies, required for some RSAT scripts
################################################################

## Installation of R packages
cd $RSAT

make -f makefiles/install_rsat.mk install_r_packages

# ## The command R CMD INSTALL apparently does not work at this stage.
# ##	root@rsat-tagc:/workspace/rsat# R CMD INSTALL reshape
# ##	Warning: invalid package 'reshape'
# ##	Error: ERROR: no packages specified
# cd $RSAT; make -f makefiles/install_rsat.mk  r_modules_list 
# ### I install them from the R interface. This should be revised to
# ### make it from the bash, but I need to see how to specify the CRAN
# ### server from the command line (for the time being, I run R and the
# ### programm asks me to specify my preferred CRAN repository the first
# ### time I install packages).
# R
# ## At the R prompt, type the following R commands.
# ##
# ## Beware, the first installation of bioconductor may take a while,
# ## because there are many packages to install
# ##
# ## Note: since this is the first time you install R packages on this
# ## VM, you need to choose a RCRAN server nearby to your site.
# install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras", "dendextend"))
# source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
# quit()
# ## At prompt "Save workspace image? [y/n/c]:", answer "n"

## Check remaining disk space
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_R_packages_installed.txt
