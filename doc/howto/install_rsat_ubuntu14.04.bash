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


################################################################
## Must be executed as root. If you are non-root but sudoer user, you
## can become it withn "sudo bash"


## Configuration for the installation
export INSTALLER=apt-get
export INSTALLER_OPT="--quiet --assume-yes"
## alternative: INSTALLER=aptitude
export INSTALL_ROOT_DIR=/bio/
export RSAT_HOME=${INSTALL_ROOT_DIR}/rsat
#export RSAT_DISTRIB=rsat_2014-08-22.tar.gz
#export RSAT_DISTRIB_URL=http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/${RSAT_DISTRIB}

## We need to update apt-get, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu

## Create a separate directory for RSAT, which must be readable by all
## users (in particular by the apache user)
mkdir -p ${INSTALL_ROOT_DIR}
cd ${INSTALL_ROOT_DIR}
mkdir -p ${INSTALL_ROOT_DIR}/install_logs
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_start.txt
apt-get update
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_apt-get_updated.txt

## Install aptitude, more efficient than apt-get to treat dependencies
## when installing and uninstalling packages.
## TO SAVE SPACE, I SUPPRESS aptitude
## apt-get install aptitude

${INSTALLER} ${INSTALLER_OPT} upgrade
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${INSTALLER}_upgraded.txt

## Packages to be checked: to I really need this ?
PACKAGES_OPT="
	curl
	yum
	php-elisp
	libgdbm-dev
	libgd2-xpm-dev
	libxml2-dev
	links
	gfortran
	libmysqlclient-dev
	texlive-latex-base
	python-virtualenv
	ipython
	ipython-notebook
	libssl-dev
	libreadline-gplv2-dev:i386
	lib64readline-gplv2-dev:i386
	libreadline-gplv2-dev
	libx11-dev
	libxt-dev
	libcurl4-openssl-dev
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
	libnet-ssleay-perl \
	libcrypt-ssleay-perl \
"

PACKAGES="
	ssh
	git
	cvs
	wget
	zip
	unzip
	finger
	screen
	make
	g++
	apache2
	php5
	libapache2-mod-php5
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
	python3
	python3-pip
	python3-setuptools 
	python3-numpy
	python3-scipy
	python3-matplotlib
	r-base-core
	emacs
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
"


PACKAGES_PERL_PROBLEM="
	libalgorithm-cluster-perl
	digest-md5-file-perl
	liblockfile-simple
	libutil-properties-perl
"
#E: Unable to locate package libalgorithm-cluster-perl
#E: Unable to locate package digest-md5-file-perl
#E: Unable to locate package liblockfile-simple
#E: Unable to locate package libutil-properties-perl

## Install the apt-get libraries
echo "Packages to be installed with ${INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES}"
echo "Perl module packages to be installed with ${INSTALLER} ${INSTALLER_OPT}"
echo "${PACKAGES_PERL}"
for LIB in ${PACKAGES} ${PACKAGES_PERL}; \
do \
   echo "`date '+%Y/%m/%d %H:%M:%S'`  installing apt-get library ${LIB}" ; \
   ${INSTALLER} install ${INSTALLER_OPT} ${LIB} > ${INSTALL_ROOT_DIR}/install_logs/${INSTALLER}_install_${LIB}.txt ; \
   df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_${LIB}_installed.txt ; \
done
echo "Log files are in folder ${INSTALL_ROOT_DIR}/install_logs"

################################################################
## Specific treatment for some Python libraries
##
## A fix for a problem to install scipy with pip: use ${INSTALLER} build-dep 
## taken from here: http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
## Note that these dependencies cost 400Mb ! To be checked
${INSTALLER} ${INSTALLER_OPT} build-dep python-numpy python-scipy
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_numpy-scipy_dependencies_installed.txt

################################################################
## To free space, remove apt-get packages that are no longer required.a
${INSTALLER} ${INSTALLER_OPT}  autoremove
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_autoremoved.txt
${INSTALLER} ${INSTALLER_OPT}  clean
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_cleaned.txt
## This really helps: it saves several hundreds Mb

## DONE: installation of Ubuntu packages
################################################################



################################################################
## Activate the Apache Web server
##
## !!!!!!!! SOME MANUAL INTERVENTION IS REQUIRED HERE  !!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
emacs -nw /etc/apache2/sites-available/000-default.conf

## Uncomment the following line:
# Include conf-available/serve-cgi-bin.conf

emacs -nw /etc/apache2/mods-available/mime.conf
## In the file /etc/apache2/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi

## The following lines are required to activate cgi scripts.  Found at
## http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/
chmod 755 /usr/lib/cgi-bin
chown root.root /usr/lib/cgi-bin
a2enmod cgi ## this is apparently required to enable cgi

## Restart the apache server to take the new config into account
service apache2 restart

## DONE: apache server configured and started
## You can check it by opening a Web connection to 
## http://[IP]
################################################################

################################################################
## Install some python libraries with pip
##
## Note: at this stage, numpy, scipy and matplotlib have already been
## installed with apt-get under Ubuntu. For other OS, they should be
## added to the pip installation
pip install soappy
pip install fisher
## pip install pygraphviz ## OSError: Error locating graphviz.

#${INSTALLER} install python3-suds
## PROBLEM : No distributions at all found for python-suds
## pip3 install python-suds

## Failures: no distributions at all found
# pip3 install wsdl
# pip3 install wstools
pip3 install fisher
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## PROBLEM: soappy seems to be discontnued for python3 !
#      pip3 install soappy
## Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/wstools
## Storing debug log for failure in /home/rsat/.pip/pip.log
##
## I also tried with easy_install3
#       easy_install3 soappy
## SAME ERROR: ImportError: No module named 'WSDLTools'

## I should test one of the following SOAP packages
pip3 install suds-jurko
pip3 install pysimplesoap

df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_pip_libraries_installed.txt

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
#
# apt-get update
#
# apt-get --quiet --assume-yes install libmodule-build-perl
# apt-get --quiet --assume-yes install libsoap-wsdl-perl
#
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
# ## Solution proposed here: http://stackoverflow.com/questions/3489642/dependency-problem-of-perl-cpan-modules
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
################       RSAT installation        ################
################################################################

## Simplified protocol: get RSAT as tar archive rather than git
mkdir -p ${RSAT_HOME}
cd ~; ln -fs ${RSAT_HOME} rsat
cd ${INSTALL_ROOT_DIR}


################################################################
## Download RSAT distribution

## Note: the git distribution requires an account at the ENS git
## server, which is currently only possible for RSAT developing team.
## In the near future, I envisage to use git also for the end-user
## distribution.
git clone git@depot.biologie.ens.fr:rsat

## For users who don't have an account on the RSAT git server, the
## code can be downloaded as a tar archive from the Web site.
# wget ${RSAT_DISTRIB_URL}
# tar -xpzf ${RSAT_DISTRIB}
# rm -f   ${RSAT_DISTRIB} ## To free space
# df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_downloaded.txt

## Metabolic pathway tools installation
##
## TO BE DONE LATER
## wget http://rsat.ulb.ac.be/~jvanheld/rsat_distrib/metabolic-tools_20110408.tar.gz

# ## DONE : obtained RSAT distribution
# ################################################################

## Run the configuration script, to specify the environment variables.
cd ${RSAT_HOME}
perl perl-scripts/configure_rsat.pl

## Load the (updated) RSAT environment variables
source RSAT_config.bashrc

## Check that the RSAT environment variable has been properly configured
echo ${RSAT}

## Initialise RSAT folders
make -f makefiles/init_rsat.mk init
    
################################################################
## For the next operations, we need to be su

################################################################
## Link the RSAT bash configuration file to a directory where files
## are loaded by each user at each login. Each user will then
## automatically load the RSAT configuration file when opening a bash
## session.
rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/

################################################################
## Installation of Perl modules required for RSAT
################################################################
## Notes
##
## 1) limxml2-dev is required to compile the Perl module XML::LibXML
##        sudo apt-get install limxml2-dev 
##
## 2) For some modules, installation failed until I used "force"
##	 force install SOAP::WSDL
##	 force install SOAP::Lite
##
## 3) Problem of dependency when installint REST::Client, even with "force".
##	 MCRAWFOR/REST-Client-271.tar.gz              : make_test FAILED but failure ignored because 'force' in effect
##	 NANIS/Crypt-SSLeay-0.64.tar.gz               : make NO
## To solve it I found this
## 	sudo apt-get install libnet-ssleay-perl
##	sudo apt-get install libcrypt-ssleay-perl
##	sudo apt-get install libssl-dev


## Set working directory to RSAT
cd $RSAT

## Get the list of Perl modules to be installed
make -f makefiles/install_rsat.mk  perl_modules_list

## Check which Perl modules are already installed
make -f makefiles/install_rsat.mk perl_modules_check
## The result file name will be displayed at the end of the tests

## Check the result of perl modules test
more check_perl_modules_eval.txt
## On Ubuntu 14.04, Object::InsideOut has status "Fail" but there is
## apparently no problem

## Identify Perl modules that were not OK after the ubuntu package installation
grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;'
MISSING_PERL_MODULES=`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | cut -f 2 | xargs`
echo ${MISSING_PERL_MODULES}

## Beware: the _noprompt suffix is
## optional. It has the advantage to avoid for the admin to confirm
## each installation step, but the counterpart is that errors may be
## overlooked.
make -f makefiles/install_rsat.mk perl_modules_install PERL_MODULES="${MISSING_PERL_MODULES}"

## Check if all required Perl modules have now been correctly installed
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt

## Note: I had to force installation for the some modules, because
## there seem to be some circular dependencies.
grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v "Object::InsideOut"
make -f makefiles/install_rsat.mk perl_modules_install_by_force PERL_MODULES_TO_FORCE="`grep -v '^OK'  check_perl_modules_eval.txt | grep -v '^;' | grep -v Object::InsideOut| cut -f 2 | xargs`"


## Last check for Perl modules. 
## If some of them still fail (except Object::InsideOut), manual intervention will be required.
make -f makefiles/install_rsat.mk perl_modules_check
more check_perl_modules_eval.txt

## Measure remaining disk space
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_perl_modules_installed.txt


################################################################
## Install selected R librairies, required for some RSAT scripts
################################################################

## Installation of R packages

## The command R CMD INSTALL apparently does not work at this stage.
##	root@rsat-tagc:/workspace/rsat# R CMD INSTALL reshape
##	Warning: invalid package 'reshape'
##	Error: ERROR: no packages specified

cd $RSAT; make -f makefiles/install_rsat.mk  r_modules_list 

### I install them from the R interface. This should be revised to
### make it from the bash, but I need to see how to specify the CRAN
### server from the command line (for the time being, I run R and the
### programm asks me to specify my preferred CRAN repository the first
### time I install packages).
R

## At the R prompt, type the following R commands.
## Beware, the first installation of bioconductor may take a while, because there are many packages to install
install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras", "dendextend"))
source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
quit()
## At prompt "Save workspace image? [y/n/c]:", answer "n"


## Check remaining disk space
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_R_packages_installed.txt

################################################################
## Configure RSAT web server

## Edit the file to replace [RSAT_PARENT_FOLDER] byt the parent directory
## of the rsat directory.
cd ${RSAT}
rsync -ruptvl RSAT_config.conf /etc/apache2/sites-enabled/rsat.conf
apache2ctl restart

## You should now test the access to the RSAT Web server, whose URL is
## in the environment variable RSAT_WWW
echo $RSAT_WWW

## If the value is "auto", get the URL as follows
# export IP=`ifconfig eth0 | awk '/inet /{print $2}' | cut -f2 -d':'
# export RSAT_WWW=http://${IP}/rsat/
# echo $RSAT_WWW

################################################################
## Next steps require to be done as rsat administrator user

## compile RSAT programs written in C
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_app_compiled.txt

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!  BUG    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!! I HAVE A PROBLEM TO COMPILE KWALKS. SHOULD BE CHECKED !!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Install some third-party programs required by some RSAT scripts.
make -f makefiles/install_software.mk install_ext_apps
df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_extapp_installed.txt

## ONLY FOR THE IFB CLOUD: 
##
## replace the data directory by a link to
## a separate disk containing all RSAT data.
export RSAT_DATA_DIR=/root/mydisk/rsat_data
cd ${RSAT}/public_html
mv data/* ${RSAT_DATA_DIR}/
mv data/.htaccess ${RSAT_DATA_DIR}/
rmdir data
ln -s ${RSAT_DATA_DIR} data
cd $RSAT

## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae
download-organism -v 1 -org Escherichia_coli_K_12_substr__MG1655_uid57779

## Optionally, install some pluricellular model organisms
# download-organism -v 1 -org Drosophila_melanogaster
# download-organism -v 1 -org Caenorhabditis_elegans
# download-organism -v 1 -org Arabidopsis_thaliana

## Get the list of organisms supported on your computer.
supported-organisms

df -m > ${INSTALL_ROOT_DIR}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_organism_installed.txt

################################################################
## At this stage you can already check some simple RSAT command 

## Test a simple Perl script that does not require for organisms to be
## installed.
which random-seq
random-seq -l 100

## Test a simple python script that does not require organisms to be
## installed.
random-motif -l 10 -c 0.90

################
## Test some external programs

## vmatch (used in purge-sequence)
random-seq -l 100 | purge-sequence

## get the help for seqlogo
which seqlogo
seqlogo

## ghostscript
which gs
gs --version


################################################################
## Configure the SOAP/WSDL Web services


## Adapt the URL to your local configuration.
cd $RSAT
make -f makefiles/init_rsat.mk ws_init

## After this, you should re-generate the web services stubb, with the
## following command.
make -f makefiles/init_rsat.mk ws_stub

## Test the local web services
make -f makefiles/init_rsat.mk ws_stub_test

## Test RSAT Web services (local and remote) without using the
## SOAP/WSDL stubb (direct parsing of the remote WSDL file)
make -f makefiles/init_rsat.mk ws_nostub_test

## Test the program supported-organisms-server, which relies on Web
## services without stub
supported-organisms-server -url ${RSAT_WS}
supported-organisms-server -url ${RSAT_WS} | wc
supported-organisms-server -url http://localhost/rsat/ | wc
supported-organisms-server -url http://rsat.eu/ | wc

################################################################
## Install the cluster management system (torque, qsub, ...)

## Check the number of core (processors)
grep ^processor /proc/cpuinfo

## Check RAM
grep MemTotal /proc/meminfo
