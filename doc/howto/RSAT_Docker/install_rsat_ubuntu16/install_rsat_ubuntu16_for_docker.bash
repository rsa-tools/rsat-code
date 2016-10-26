################################################################
## Instructions used to install a Virtual Machine on the IFB cloud
## (Institut Francais de Bioinformatique), and on a VirtualBox VM.

## THIS IS NOT REALLY A bahsrc FILE, IT IS A SUCCESSION OF
## INSTRUCTIONS AND COMMENTS, THAT SHOULD BE DONE MANUALLY. I (JvH)
## SHOULD IMPROVE THIS WHEN I CAN.

## TO RUN Virtualbox on Mac OSX
#    Go into your Ubuntu Settings
#    Navigate to Keyboard → Keyboard Layout Settings
#    add English (Macintosh)


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##
## For the installation of Perl package and for third-party Linux
## packages, I need to temporarily inactivate the antivirus software
## Kaspersky.
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## For Debian, I must set the locales manually
# export LANGUAGE=en_US.UTF-8
# export LANG=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# locale-gen en_US.UTF-8
# dpkg-reconfigure locales

## Path for the local installation
export RSAT_PARENT_PATH=/packages
export RSAT_HOME=${RSAT_PARENT_PATH}/rsat

## URL to download the RSAT distribution
export RSAT_RELEASE=2016-10-26 ## Version to be downloaded from the tar distribution
export RSAT_ARCHIVE=rsat_${RSAT_RELEASE}.tar.gz
export RSAT_DISTRIB_URL=http://pedagogix-tagc.univ-mrs.fr/download_rsat/${RSAT_ARCHIVE}

## Configuration for the installation
export INSTALLER=apt-get
export INSTALLER_OPT="--quiet --assume-yes"
## alternative: INSTALLER=aptitude

################################################################
## Must be executed as root. If you are non-root but sudoer user, you
## can become it withn "sudo bash"

################################################################
## Before anything else, check that the date, time and time zone are
## correctly specified
date

## If not, set up the time zone, date and time with this command
## (source: https://help.ubuntu.com/community/UbuntuTime).
# dpkg-reconfigure tzdata


## We need to update apt-get, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu

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

## We can then check the increase of disk usage during the different
## steps of the installation
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt


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
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_apt-get_updated.txt

################################################################
## I tried to run dist-upgrade because it 'can "intelligently" handle
## changes in the dependencies system. This includes removing packages
## that are no longer necessary or resolve conflicts between packages
## that arose because of changes in the dependencies.'
## http://askubuntu.com/questions/194651/why-use-apt-get-upgrade-instead-of-apt-get-dist-upgrade
## HOWEVER, THE PERL UPDATE DOES NOT WORK ANYMORE AFTER THAT !!!
${INSTALLER} ${INSTALLER_OPT} upgrade
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${INSTALLER}_upgrade.txt
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt



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
"

PACKAGES="
apt-utils
ssh
git
cvs
wget
zip
unzip
screen
make
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
finger
"

################################################################
## apt-get packages to install Perl modules (not properly speaking
## necessary, could be done with cpan, but ensure consistency with
## ubuntu OS)

PACKAGES_PERL="perl-doc
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
echo "Packages to be installed with ${INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES}"
echo "Perl module packages to be installed with ${INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES_PERL}"
for LIB in ${PACKAGES} ${PACKAGES_PERL}; \
do \
   echo "`date '+%Y/%m/%d %H:%M:%S'`  installing apt-get library ${LIB}" ; \
   ${INSTALLER} install ${INSTALLER_OPT} ${LIB} > ${RSAT_PARENT_PATH}/install_logs/${INSTALLER}_install_${LIB}.txt ; \
   df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${LIB}_installed.txt ; \
done
echo "Log files are in folder ${RSAT_PARENT_PATH}/install_logs"
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt


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
${INSTALLER} install ${INSTALLER_OPT} console-data

################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## DOES NOT WORK ON Ubunutu 16.04.1 anmore
## Not sure it is still required though.
## TO BE CHECKED
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ################################################################
# ## Specific treatment for some Python libraries
# ##
# ## A fix for a problem to install scipy with pip: use ${INSTALLER} build-dep 
# ## taken from here: http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
# ## Note that these dependencies cost 400Mb ! To be checked
# ${INSTALLER} ${INSTALLER_OPT} build-dep python-numpy python-scipy
# df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_numpy-scipy_dependencies_installed.txt
# grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt

################################################################
## To free space, remove apt-get packages that are no longer required.a
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt
${INSTALLER} ${INSTALLER_OPT}  autoremove
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_autoremoved.txt
${INSTALLER} ${INSTALLER_OPT}  clean
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned.txt
## This really helps: it saves several hundreds Mb
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt

## DONE: installation of Ubuntu packages
################################################################


################################################################
## Install some python libraries with pip
##
## Note: at this stage, numpy, scipy and matplotlib have already been
## installed with apt-get under Ubuntu. For other OS, they should be
## added to the pip installation
pip install soappy
pip install fisher
pip install httplib2
## pip install pygraphviz ## OSError: Error locating graphviz.

## optional: an utility to measure internet bandwidth
pip install speedtest-cli

#${INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-suds

## Failures: no distributions at all found
# pip3 install wsdl
# pip3 install wstools
pip3 install fisher
pip3 install snakemake
pip3 install rpy2  ## THIS FAILS on the IFB cloud. To be checked.
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages
pip3 install suds-jurko
pip3 install pysimplesoap
pip3 install soappy

## Check disk usage
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt


# ################################################################
# ## TO BE CHECKED: TO WE STILL NEED TO DO ALL THE TRICKY STUFF BELOW ?
# ## NOT SURE: on 2014/08/15 Jv installed a VM on IFB cloud, and the
# ## SOAP/WSDL Web services seem to work without this.
#
# ## The installation of SOAP:WSDL under cpan is particularly tricky. 
# ## In Ubuntu, there is a way to install it with ${INSTALLER}. 
# ## http://www.installion.co.uk/ubuntu/trusty/universe/l/libsoap-wsdl-perl/fr/install.html
# emacs -nw /etc/apt/sources.list
#
# ## Ensure that the following line is set to "universe"
# deb http://us.archive.ubuntu.com/ubuntu trusty main universe
# ## You can now quit emacs

apt-get update
apt-get --quiet --assume-yes install libmodule-build-perl
apt-get --quiet --assume-yes install libsoap-wsdl-perl

# ## Note: this is still not sufficient to get SOAP::WSDL to run the two
# ## following targets
# ##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stubb
# ##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stubb_test
#
# ## We first need to fix some problem with CPAN : on Ubuntu, I cannot
# ## install the SOAP::WSDL module, which is required for several
# ## functionalities. More precisely, after fiddling around for a few
# ## hours, the server is able to answer to web services requests, but I
# ## cannot run clients on it. Since NeAT relies on WSDL clients, it is
# ## impossible to have neat running on and Ubuntu server. The
# ## installation is however possible, since the stubb can be generated
# ## on rsat-tagc.univ-mrs.fr.  I have no idea how we did to install
# ## SOAP::WSDL there. In any case, the 
# ##
# ## Solution proposed here: http://stackoverflow.com/questions/3489642/dependency-problem-of-perl-cpan -modules
# ## Not sure it works by its own, but cannot harm.
# cpan
# ## At the cpan prompt, type the following
# install Module::Build
# install Module::Build::Compat
# install CPAN ## This takesa HUGE time. I answer all questions by the default answer (simply type the Enter key)
# upgrade ## Takes a HUGE time, since all packages are apparently re-tested
# quit
##
## DONE: tricky stuff
################################################################


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

################################################################
## Activate the Apache Web server


## Make a backup of original apache configuration
zip -ry apacheconf_backup_`date '+%Y-%m-%d_%H%M%S'`.zip /etc/apache2/

## We need to replace some lines in the Apache configuration files
# nano /etc/apache2/sites-available/000-default.conf
## Uncomment the following line:
# Include conf-available/serve-cgi-bin.conf
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#Include conf-available/serve-cgi-bin.conf|Include conf-available/serve-cgi-bin.conf|' /etc/apache2/sites-available/000-default.conf

## To avoid puzzling warning at apache start, set ServerName globally.
echo "" >> /etc/apache2/apache2.conf
echo "ServerName localhost" >> /etc/apache2/apache2.conf

## In the file /etc/apache2/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi
##
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddHandler cgi-script .cgi|AddHandler cgi-script .cgi|' /etc/apache2/mods-available/mime.conf
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddEncoding x-|AddEncoding x-|' /etc/apache2/mods-available/mime.conf

# echo "" >>  /etc/apache2/mods-available/mime.conf
# echo "        ## RSAT options" >>  /etc/apache2/mods-available/mime.conf
# echo "        AddType text/plain .fasta" >>  /etc/apache2/mods-available/mime.conf
# echo "        AddType text/plain .bed" >>  /etc/apache2/mods-available/mime.conf

## Optional : also associate a plain/text mime type to extensions for
## some classical bioinformatics files.
##   AddType text/plain .fasta
##   AddType text/plain .bed
## I also uncomment the following, for convenience
##        AddEncoding x-compress .Z
##        AddEncoding x-gzip .gz .tgz
##        AddEncoding x-bzip2 .bz2



## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Adapt the PHP parameters
## TO DO: adapt for Ubuntu 16.04.1, where php5 has been replaced by php7. 
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# emacs -nw /etc/php5/apache2/php.ini
## Modify the following parameters
##      post_max_size = 100M
##      upload_max_filesize=200M
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## The following lines are required to activate cgi scripts.  Found at
## http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/
chmod 755 /usr/lib/cgi-bin
chown root.root /usr/lib/cgi-bin
a2enmod cgi ## this is apparently required to enable cgi
a2dismod mpm_event
a2enmod mpm_prefork

## THIS DOES NOT BELONG TO THE INSTALLER BUT TO THE RUN / TEST
#service apache2 restart

## DONE: apache server configured and started
## You can check it by opening a Web connection to 
## http://[IP]
################################################################

################################################################
################       RSAT installation        ################
################################################################


################################################################
## Download RSAT distribution

## Note: the git distribution requires an account at the ENS git
## server, which is currently only possible for RSAT developing team.
## In the near future, we may use git also for the end-user
## distribution.

cd ${RSAT_PARENT_PATH}
echo "Downloading RSAT from  ${RSAT_DISTRIB_URL}"
wget --no-clobber ${RSAT_DISTRIB_URL}
tar -xpzf ${RSAT_ARCHIVE}
cd ${RSAT_HOME}

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
ln -fs ${RSAT_HOME} ${HOME}/rsat

## Metabolic pathway tools installation
##
## TO BE DONE LATER
## wget http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/metabolic-tools_20110408.tar.gz

# ## DONE : obtained RSAT distribution
# ################################################################

## Run the configuration script, to specify the environment variables.
cd ${RSAT_HOME}
perl perl-scripts/configure_rsat.pl --auto rsat_site=rsat-vb-2016-10 \
    ucsc_tools=0 ensembl_tools=0 phylo_tools=0 compara_tools=0 \
    rsat_server_admin=RSAT_admin

## Parameters to change
##   rsat_site   rsat-vm-2016-03
##   rsat_server_admin    I don't specify it, because I don't want to receive notifications from all the VMs
## I activate the optional tools ucsc_tools and ensembl_tools, but not the other ones because they require many genomes (phylo tools) or big genomes (compara_tools, variation_tools).

## Load the (updated) RSAT environment variables
source RSAT_config.bashrc

## Check that the RSAT environment variable has been properly configured
echo ${RSAT}

## Initialise RSAT folders
make -f makefiles/init_rsat.mk init

################################################################
## Previous way to specify bashrc parameters, via
## /etc/bash_completion.d/. I change it (2014-09-23) because it does
## not allow to run remote commands via ssh (/etc/bash_completion.d is
## apparently only loaded in interactive mode).
## 
## Link the RSAT bash configuration file to a directory where files
## are loaded by each user at each login. Each user will then
## automatically load the RSAT configuration file when opening a bash
## session.
rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/
## ln -fs ${RSAT_HOME}/RSAT_config.bashrc /etc/bash_completion.d/

#emacs -nw /etc/bash.bashrc




################################################################
## Installation of Perl modules required for RSAT
################################################################
## Notes
##
## 1) limxml2-dev is required to compile the Perl module XML::LibXML
# apt-get install limxml2-dev 
## 2) For some modules, installation failed until I used "force"
##	 force install SOAP::WSDL
##	 force install SOAP::Lite
##
## 3) Problem of dependency when installint REST::Client, even with "force".
##	 MCRAWFOR/REST-Client-271.tar.gz              : make_test FAILED but failure ignored because 'force' in effect
##	 NANIS/Crypt-SSLeay-0.64.tar.gz               : make NO
## To solve it I found this
## 	apt-get install libnet-ssleay-perl
##	apt-get install libcrypt-ssleay-perl
##	apt-get install libssl-dev


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## NOTE: at this stage I had to inactivate the Kaspersky antivirus because it
## prevents CPAN from downloading and installing Perl modules !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cpan YAML ## Answer "yes" !!!
cpan CPAN ## Update CPAN

## Set working directory to RSAT
cd ${RSAT}

## Get the list of Perl modules to be installed
make -f makefiles/install_rsat.mk  perl_modules_list

## Check which Perl modules are already installed
make -f makefiles/install_rsat.mk perl_modules_check
## The result file name will be displayed at the end of the tests

## Check the result of perl modules test
more check_perl_modules_eval.txt
## On Ubuntu 14.04, Object::InsideOut has status "Fail" but there is
## apparently no problem

## Check missing modules
grep Fail  check_perl_modules_eval.txt

## Identify Perl modules that were not OK after the ubuntu package installation
grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;'
MISSING_PERL_MODULES=`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | cut -f 2 | xargs`
echo "Missing Perl modules:     ${MISSING_PERL_MODULES}"

## Beware: the _noprompt suffix is optional. It has the advantage to
## avoid for the admin to confirm each installation step, but the
## counterpart is that errors may be overlooked.
make -f makefiles/install_rsat.mk perl_modules_install PERL_MODULES="${MISSING_PERL_MODULES}"

## Check if all required Perl modules have now been correctly installed
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt
## Note: Object::InsideOut should be ignored for this test, because it
## always display "Fail", whereas it is OK during installation.

## Note: I had to force installation for the some modules, because
## there seem to be some circular dependencies.
grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v "Object::InsideOut"
make -f makefiles/install_rsat.mk perl_modules_install_by_force PERL_MODULES_TO_FORCE="`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v Object::InsideOut| cut -f 2 | xargs`"


## Last check for Perl modules. 
## If some of them still fail (except Object::InsideOut), manual intervention will be required.
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt

## Measure remaining disk space
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_perl_modules_installed.txt
grep ${DEVICE} ${RSAT_PARENT_PATH}/install_logs/df_*.txt


################################################################
## Configure RSAT web server

## Edit the file to replace [RSAT_PARENT_FOLDER] byt the parent directory
## of the rsat directory.
cd ${RSAT}
rsync -ruptvl RSAT_config.conf /etc/apache2/sites-enabled/rsat.conf

## OPTIONAL: since I am using this to install a virtual machine whose
## only function will be to host the RSAT server, I replace the normal
## default web folder by RSAT web folder. 
##
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|DocumentRoot /var/www/html|DocumentRoot /packages/rsat/public_html|' /etc/apache2/sites-available/000-default.conf

apache2ctl restart
## The server will now immediately display RSAT home page when you
## type its IP address.

## You should now test the access to the RSAT Web server, whose URL is
## in the environment variable RSAT_WWW
echo $RSAT_WWW


################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## THIS DOES NOT WORK ON VirtualBox VM 16.04.1
## CHECK IF I STILL NEED IT
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ## If the value is "auto", get the URL as follows
# export IP=`ifconfig eth0 | awk '/inet /{print $2}' | cut -f2 -d':'`
# # export IP=192.168.56.101
# echo ${IP}
# export RSAT_WWW=http://${IP}/rsat/
# echo $RSAT_WWW

################################################################
## Next steps require to be done as rsat administrator user

## compile RSAT programs written in C
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_app_compiled.txt

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!  BUG    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!! I HAVE A PROBLEM TO COMPILE KWALKS. SHOULD BE CHECKED !!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae \
 -org Escherichia_coli_K_12_substr__MG1655_uid57779

## Optionally, install some pluricellular model organisms
# download-organism -v 1 -org Drosophila_melanogaster
# download-organism -v 1 -org Caenorhabditis_elegans
# download-organism -v 1 -org Arabidopsis_thaliana

## Get the list of organisms supported on your computer.
supported-organisms

df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_organism_installed.txt

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
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_R_packages_installed.txt

################################################################
## Install some third-party programs required by some RSAT scripts.
cd ${RSAT}
make -f makefiles/install_software.mk
make -f makefiles/install_software.mk install_ext_apps
df -m > ${RSAT_PARENT_PATH}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_extapp_installed.txt

################################################################
## Configure the SOAP/WSDL Web services

## Check the URL of the web services (RSAT_WS). By default, the server
## addresses the WS requests to itself (http://localhost/rsat) because
## web services are used for multi-tierd architecture of some Web
## tools (retrieve-ensembl-seq, NeAT).
cd $RSAT

## Set a fixed IP address for the web services
perl ${RSAT}/perl-scripts/configure_rsat.pl auto RSAT_WS=auto RSAT_WWWW=auto SUDO=''
# rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/ ## HERE OR BELOW ???

################################################################
## Install the Web services
cd ${RSAT}
make -f makefiles/init_rsat.mk ws_param  ## Check parameters to generate the stub for the Web services
make -f makefiles/init_rsat.mk ws_init  ## Initialize the stub for the Web services
make -f makefiles/init_rsat.mk ws_stub  ## Generate the stub for the Web services
make -f makefiles/init_rsat.mk ws_stub_test  ## Test the stub for the Web services
make -f makefiles/init_rsat.mk ws_nostub_test  ## Test Web services with no stub


################################################################
## tests on the Web site

## Run the demo of the following tools
##
## - retrieve-seq to check the access to local genomes (at least
##   Saccharomyces cerevisiae)
##
## - feature-map to check the GD library
##
## - retrieve-ensembl-seq to check the interface to Ensembl
##
## - fetch-sequences to ceck the interface to UCSC
##
## - some NeAT tools (they rely on web services)
##
## - peak-motifs because it mobilises half of the RSAT tools -> a good
##   control for the overall installation.
##
## - footprint-discovery to check the tools depending on homology
##   tables (blast tables).




################################################################
########### Remove unnecessary files to clean space  ###########
###########   BEFORE DELIVERY for VirtualBox         ###########
################################################################

cd ${RSAT_PARENT_PATH}; rm -f  ${RSAT_ARCHIVE} ## Free space
cd ${RSAT}; make -f makefiles/server.mk clean_tmp ## Clean temporary directory
rm -rf ${RSAT}/app_sources
rm -rf ${RSAT}/public_html/tmp/serialized_genomes
find   ${RSAT_HOME}/public_html/tmp/  -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} \;
rm -rf public_html/tmp/www-data ## Clean Apache user temporary directory
rm -rf public_html/tmp/serialized_genomes ## Clean serialized organisms
rm -rf  ~/.ssh/ ## Make sure that the installer ssh config is not included in the distribution


## THE INSTALLATION OF THE RSAT SERVER IS HOW DONE. THE REST IS OPTIONNAL
################################################################

