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

## Must be executed as root
sudo bash

## We need to update apt-ge, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu
apt-get update

## Install aptitude, more efficient than apt-get to treat dependencies
## when installing and uninstalling packages.
apt-get install aptitude

aptitude --quiet --assume-yes upgrade

## Enable shell connection from the host
aptitude --quiet --assume-yes install ssh

## Concurrent versioning systems
aptitude --quiet --assume-yes install git
aptitude --quiet --assume-yes install cvs

## Web aspirators
aptitude --quiet --assume-yes install wget
aptitude --quiet --assume-yes install curl

## Utilities
aptitude --quiet --assume-yes install zip
aptitude --quiet --assume-yes install unzip
aptitude --quiet --assume-yes install finger
aptitude --quiet --assume-yes install screen
aptitude --quiet --assume-yes install make
aptitude --quiet --assume-yes install g++
aptitude --quiet --assume-yes install yum

## Cannot survive without emacs
aptitude --quiet --assume-yes install emacs

## Apache and utilities
aptitude --quiet --assume-yes install apache2
aptitude --quiet --assume-yes install php5
aptitude --quiet --assume-yes install libapache2-mod-php5
aptitude --quiet --assume-yes install php-elisp

## Graphic libraries and software tools
aptitude --quiet --assume-yes install libgdbm-dev
aptitude --quiet --assume-yes install libgd-tools
aptitude --quiet --assume-yes install libgd-gd2-perl
aptitude --quiet --assume-yes install libgd2-xpm-dev
aptitude --quiet --assume-yes install libxml2-dev

aptitude --quiet --assume-yes install libnet-ssleay-perl
aptitude --quiet --assume-yes install libcrypt-ssleay-perl
aptitude --quiet --assume-yes install libssl-dev
aptitude --quiet --assume-yes install ghostscript
aptitude --quiet --assume-yes install gnuplot
aptitude --quiet --assume-yes install graphviz
aptitude --quiet --assume-yes install lib32z1
aptitude --quiet --assume-yes install lib32ncurses5
aptitude --quiet --assume-yes install lib32bz2-1.0

## Text-mode Web browser, used by some packages
aptitude --quiet --assume-yes install links

## Some linux packages required for R BioConductor
aptitude --quiet --assume-yes install libc6-dev
aptitude --quiet --assume-yes install gfortran
aptitude --quiet --assume-yes install build-essential
aptitude --quiet --assume-yes install libreadline-gplv2-dev:i386
aptitude --quiet --assume-yes install lib64readline-gplv2-dev:i386
aptitude --quiet --assume-yes install libreadline-gplv2-dev
aptitude --quiet --assume-yes install libx11-dev
aptitude --quiet --assume-yes install libxt-dev
aptitude --quiet --assume-yes install libcurl4-openssl-dev
aptitude --quiet --assume-yes install libxml2-dev
## BEWARE: texlive-full occupies a lot of disk space. I should check if this is really required (for R ?)
## aptitude install texlive-full
aptitude --quiet --assume-yes install tcl8.5-dev
aptitude --quiet --assume-yes install tk8.5-dev
aptitude --quiet --assume-yes install libxss-dev
aptitude --quiet --assume-yes install libpng12-dev
aptitude --quiet --assume-yes install libjpeg62-dev
aptitude --quiet --assume-yes install libcairo2-dev

## mysql client is required for ensembl client scripts
aptitude --quiet --assume-yes install mysql-client
aptitude --quiet --assume-yes install libmysqlclient-dev

## Java 
## seems to be required for SOAP::WSDL Perl module
aptitude --quiet --assume-yes install default-jre
## aptitude --quiet --assume-yes install default-jdk

## Latex is required for RSAT doc + other applications (e.g. R). Note
## that it takes a some time to install
aptitude --quiet --assume-yes install texlive-latex-base

################################################################
## Python and modules
aptitude --quiet --assume-yes install python
aptitude --quiet --assume-yes install python-setuptools 
aptitude --quiet --assume-yes install python-virtualenv
aptitude --quiet --assume-yes install python-pip
aptitude --quiet --assume-yes install python-dev
aptitude --quiet --assume-yes install python-suds

aptitude --quiet --assume-yes install ipython
aptitude --quiet --assume-yes install ipython-notebook

## A fix for a problem to install scipy with pip: use aptitude build-dep 
## taken from here: http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
aptitude --quiet --assume-yes build-dep python-numpy python-scipy

pip install numpy
## Note: the installation of scipy and matplotlib takes some time and issues
## a lot of warning messages, but finally it works
pip install scipy 
pip install matplotlib
pip install soappy
pip install fisher
## pip install pygraphviz ## OSError: Error locating graphviz.


## We need both python2.7 and python3 (for different scripts)
aptitude --quiet --assume-yes install python3
aptitude --quiet --assume-yes install python3-pip
aptitude --quiet --assume-yes install python3-setuptools 
aptitude --quiet --assume-yes install python3
aptitude --quiet --assume-yes install python3-pip 
aptitude --quiet --assume-yes install python3-dev
aptitude --quiet --assume-yes install python3-numpy
aptitude --quiet --assume-yes install python3-scipy
aptitude --quiet --assume-yes install python3-matplotlib
#aptitude install python3-suds
## For pip3 also, scipy and matplotlib return a lot of verbosity, but the installation finally works

## Problem : No distributions at all found for python-suds
## pip3 install python-suds

## Problems: 
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

################################################################
## Perl modules
APTGET_PERLMOD="perl-doc \
        pmtools \
	libyaml-perl \
	libemail-simple-perl \
	libemail-sender-perl \
	libemail-simple-creator-perl \
	libpostscript-simple-perl \
	libstatistics-distributions-perl \
	libalgorithm-cluster-perl \
	digest-md5-file-perl \
	libio-all-perl \
	liblockfile-simple \
	libobject-insideout-perl \
	libutil-properties-perl \
	libobject-insideout-perl \
	libsoap-lite-perl \
	libsoap-wsdl-perl \
	libxml-perl \
	libxml-simple-perl \
	libxml-compile-cache-perl \
	libdbi-perl \
	liblockfile-simple-perl \
	libobject-insideout-perl \
	libgd-perl \
	libdbd-mysql-perl \
	libjson-perl \
	libbio-perl-perl \
	libdigest-md5-file-perl
"
## Install the apt-get libraries
mkdir -p install_logs
echo "Installing apt-get libraries"
for LIB in $APTGET_LIBRARIES $APTGET_PERLMOD
do
   echo "`date`        installing apt-get library $LIB"
   aptitude --quiet --assume-yes install $LIB > install_logs/aptitude_install_$LIB
done
echo "Log files are in folder install_logs"


################################################################
## TO BE CHECKED: TO WE STILL NEED TO DO ALL THE TRICKY STUFF BELOW ?

## The installation of SOAP:WSDL under cpan is particularly tricky. 
## In Ubuntu, there is a way to install it with aptitude. 
## http://www.installion.co.uk/ubuntu/trusty/universe/l/libsoap-wsdl-perl/fr/install.html
emacs -nw /etc/apt/sources.list

## Ensure that the following line is set to "universe"
deb http://us.archive.ubuntu.com/ubuntu trusty main universe
## You can now quit emacs

apt-get update

apt-get --quiet --assume-yes install libmodule-build-perl
apt-get --quiet --assume-yes install libsoap-wsdl-perl

## Note: this is still not sufficient to get SOAP::WSDL to run the two
## following targets
##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stubb
##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stubb_test

## We first need to fix some problem with CPAN : on Ubuntu, I cannot
## install the SOAP::WSDL module, which is required for several
## functionalities. More precisely, after fiddling around for a few
## hours, the server is able to answer to web services requests, but I
## cannot run clients on it. Since NeAT relies on WSDL clients, it is
## impossible to have neat running on and Ubuntu server. The
## installation is however possible, since the stubb can be generated
## on rsat-tagc.univ-mrs.fr.  I have no idea how we did to install
## SOAP::WSDL there. In any case, the 
##
## Solution proposed here: http://stackoverflow.com/questions/3489642/dependency-problem-of-perl-cpan-modules
## Not sure it works by its own, but cannot harm.
cpan
## At the cpan prompt, type the following
install Module::Build
install Module::Build::Compat
install CPAN ## This takesa HUGE time. I answer all questions by the default answer (simply type the Enter key)
upgrade ## Takes a HUGE time, since all packages are apparently re-tested
quit



################################################################
## R installation + some libraries
################################################################

sudo bash

## As sudo, I edited the file /etc/apt/sources.list 
## and added the following line 
## (see instructions on http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu/)
##   deb http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu trusty/
## I then updated the apt-get packages
aptitude update

## .. and installed the R base package
aptitude --quiet --assume-yes install r-base
#aptitude --quiet --assume-yes installr-base-dev

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

## At the R prompt
install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras"))
source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
quit()


################################################################
## To free space, remove apt-get packages that are no longer required.
apt-get autoremove
apt-get clean

################################################################
################       RSAT installation        ################
################################################################

## Create a specific user for RSAT. The user is named rsat
sudo adduser rsat
## Full Name: Regulatory Sequence Analysis Tools admin

## Grant sudoer privileges to the rsat user (will be more convenient for
## installing Perl modules, software tools, etc)
visudo
## then add the following line below "User privilege specification"
# rsat    ALL=(ALL:ALL) ALL

## The installation is done under the rsat login
su - rsat

## I first recuperatemy .ssh folder from some other server (by ssh).
## Then I start an agent to manage my passphrase for ssh
## transfers. This is more convenient, I only provide my password
## once.
ssh-agent > agent
source agent
ssh-add

## Note: this is only possible for regular RSAT admin. It requires to
## have specified te RSAT ssh key and sent it to the git server at
## ENS.
git clone git@depot.biologie.ens.fr:rsat

## Run the configuration script, to specify the environment variables.
cd rsat
perl perl-scripts/configure_rsat.pl 

## Load the (updated) RSAT environment variables
source RSAT_config.bashrc

## Initialise RSAT folders
make -f makefiles/init_rsat.mk init

## Exit from the rsat session (and become root again)
exit

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



## Check that RSAT path has been defined
echo $RSAT

## Set working directory to RSAT
cd $RSAT

## Get the list of Perl modules to be installed
make -f makefiles/install_rsat.mk  perl_modules_list

## Check which Perl modules are already installed
make -f makefiles/install_rsat.mk perl_modules_check
## The result file name will be displayed at the end of the tests

## Install these modules Beware: the _noprompt suffix is optional. It
## has the advantage to avoid for the admin to confirm each
## installation step, but the counterpart is that errors may be
## overlooked.
make -f makefiles/install_rsat.mk perl_modules_install_noprompt

## Note: I had to force installation for the some modules, because
## there seem to be some circular dependencies.
make -f makefiles/install_rsat.mk perl_modules_install_by_force

## Check if all required Perl modules have been correctly installed
make -f makefiles/install_rsat.mk perl_modules_check

## Ensure that all files belong to rsat user
chown -R rsat.rsat .


################################################################
## Activate the Apache Web server and RSAT configuration
sudo emacs -nw /etc/apache2/sites-available/000-default.conf

## Activate the following line:
# Include conf-available/serve-cgi-bin.conf

## In the file /etc/apache2/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi

## From http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/
sudo chmod 755 /usr/lib/cgi-bin
sudo chown root.root /usr/lib/cgi-bin
sudo a2enmod cgi ## this is apparently required to enable cgi
apache2ctl restart

################################################################
## Configure RSAT web server

## Edit the file to replace [RSAT_PARENT_FOLDER] byt the parent directory
## of the rsat directory.
cd ${RSAT}
rsync -ruptvl RSAT_config.conf /etc/apache2/sites-enabled/rsat.conf
apache2ctl restart

################################################################
## Next steps require to be done as rsat administrator user

su - rsat

whoami 
## Should return "rsat"

## compile RSAT programs written in C
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all

## Install some third-party programs required by some RSAT scripts.
make -f makefiles/install_software.mk install_ext_apps

## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae
download-organism -v 1 -org Escherichia_coli_K_12_substr__MG1655_uid57779

## Optionally, install some pluricellular model organisms
download-organism -v 1 -org Drosophila_melanogaster
download-organism -v 1 -org Caenorhabditis_elegans
download-organism -v 1 -org Arabidopsis_thaliana

## Get the list of organisms supported on your computer.
supported-organisms

################################################################
## IMPORTANT: request a vmatch license from http://www.vmatch.de/, and
## place the license file (vmatch.lic) in the bin folder $RSAT/bin
make -f makefiles/install_software.mk install_vmatch

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
make -f makefiles/init_rsat.mk ws_stubb

## Test the local web services
make -f makefiles/init_rsat.mk ws_stubb_test

## Test RSAT Web services (local and remote) without using the
## SOAP/WSDL stubb (direct parsing of the remote WSDL file)
make -f makefiles/init_rsat.mk ws_nostubb_test

################################################################
## Install the cluster management system (torque, qsub, ...)

## Check the number of core (processors)
grep ^processor /proc/cpuinfo

## Check RAM
grep MemTotal /proc/meminfo

################################################################
##
## ATTENTION: for the IDB cloud, to ensure persistence, you
## imperatively have to run the following command in the VM before any
## shutdown onthe Web site
##
################################################################

## cd; bash cleaner.sh ; history -c && history -w && logout
