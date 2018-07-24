################################################################
## Install all the Ubuntu packages required prior to the installation
## of the Regulatory Sequence Analysis Tools (RSAT; http://rsat.eu/).

source RSAT_config.mk

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo


## PRIOR REQUIREMENT: Mac developer tools (for make and other utilities)

################################################################
## Required apt-get packages
PACKAGES_REQUIRED="
screen
libgd
ghostscript
gnuplot
graphviz
mysql
curl
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
# grep ${DEVICE} ${RSAT}/install_logs/df_*.txt

## Specific brew repository for java.
## NOTE: admin password will be required during installation. 
brew cask install java

## DONE: installation of Mac OS X packages via brew
################################################################

