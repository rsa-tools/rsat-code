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
echo "" | apt-get upgrade
echo "" | apt-get install ssh

## Concurrent versioning systems
echo "" | apt-get install git
echo "" | apt-get install cvs

## Web aspirators
echo "" | apt-get install wget
echo "" | apt-get install curl

## Utilities
echo "" | apt-get install zip
echo "" | apt-get install unzip
echo "" | apt-get install finger
echo "" | apt-get install screen
echo "" | apt-get install make
echo "" | apt-get install g++
echo "" | apt-get install yum

## Cannot survive without emacs
echo "" | apt-get install emacs

## Perl packages
echo "" | apt-get install perl-doc
echo "" | apt-get install pmtools

## Apache and utilities
echo "" | apt-get install apache2
echo "" | apt-get install php5
echo "" | apt-get install libapache2-mod-php5
echo "" | apt-get install php-elisp

## Graphic libraries and software tools
echo "" | apt-get install libgdbm-dev
echo "" | apt-get install libgd-tools
echo "" | apt-get install libgd-gd2-perl
echo "" | apt-get install libgd2-xpm-dev
echo "" | apt-get install libxml2-dev

echo "" | apt-get install libnet-ssleay-perl
echo "" | apt-get install libcrypt-ssleay-perl
echo "" | apt-get install libssl-dev
echo "" | apt-get install ghostscript
echo "" | apt-get install gnuplot
echo "" | apt-get install graphviz
echo "" | apt-get install lib32z1
echo "" | apt-get install lib32ncurses5
echo "" | apt-get install lib32bz2-1.0

## Text-mode Web browser, used by some packages
echo "" | apt-get install links

## Some linux packages required for R BioConductor
echo "" | apt-get install libc6-dev
echo "" | apt-get install gfortran
echo "" | apt-get install build-essential
echo "" | apt-get install libreadline-gplv2-dev:i386
echo "" | apt-get install lib64readline-gplv2-dev:i386
echo "" | apt-get install libreadline-gplv2-dev
echo "" | apt-get install libx11-dev
echo "" | apt-get install libxt-dev
echo "" | apt-get install libcurl4-openssl-dev
echo "" | apt-get install libxml2-dev
## BEWARE: texlive-full occupies a lot of disk space. I should check if this is really required (for R ?)
## apt-get install texlive-full
echo "" | apt-get install tcl8.5-dev
echo "" | apt-get install tk8.5-dev
echo "" | apt-get install libxss-dev
echo "" | apt-get install libpng12-dev
echo "" | apt-get install libjpeg62-dev
echo "" | apt-get install libcairo2-dev

## mysql client is required for ensembl client scripts
echo "" | apt-get install mysql-client
echo "" | apt-get install libmysqlclient-dev

## Java 
## seems to be required for SOAP::WSDL Perl module
echo "" | apt-get install default-jre
## echo "" | apt-get install default-jdk

## Latex is required for RSAT doc + other applications (e.g. R). Note
## that it takes a some time to install
echo "" | apt-get install texlive-latex-base


################################################################
## Python and modules
echo "" | apt-get install python
echo "" | apt-get install python-setuptools 
echo "" | apt-get install python-virtualenv
echo "" | apt-get install python-pip
echo "" | apt-get install python-dev
echo "" | apt-get install python-suds

echo "" | apt-get install ipython
echo "" | apt-get install ipython-notebook

echo "" | apt-get install python3
echo "" | apt-get install python3-pip

## A fix for a problem to install scipy with pip: use apt-get build-dep 
## taken from here: http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
echo "" | apt-get build-dep python-numpy python-scipy

pip install numpy
## Note: the installation of scipy and matplotlib takes some time and issues
## a lot of warning messages, but finally it works
pip install scipy 
pip install matplotlib
pip install soappy
pip install fisher
## pip install pygraphviz ## OSError: Error locating graphviz.

## We need both python2.7 and python3 (for different scripts)
echo "" | apt-get install python3-setuptools 
echo "" | apt-get install python3
echo "" | apt-get install python3-pip 
echo "" | apt-get install python3-dev
#apt-get install python3-suds
pip3 install numpy
## For pip3 also, scipy and matplotlib return a lot of verbosity, but the installation finally works
pip3 install scipy
pip3 install matplotlib

## Problem : No distributions at all found for python-suds
## pip3 install python-suds

## Problems: 
# pip3 install wsdl
# pip3 install wstools
pip3 install fisher
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## soappy seems to be discontnued for python3 !
# pip3 install soappy
## I should test one of the following SOAP packages
pip3 install suds-jurko
pip3 install pysimplesoap


################################################################
## Perl modules

## The installation of SOAP:WSDL under cpan is particularly tricky. 
## In Ubuntu, there is a way to install it with apt-get. 
## http://www.installion.co.uk/ubuntu/trusty/universe/l/libsoap-wsdl-perl/fr/install.html
emacs -nw /etc/apt/sources.list

## Ensure that the following line is set to "universe"
deb http://us.archive.ubuntu.com/ubuntu trusty main universe
## You can now quit emacs

apt-get update

echo "" | apt-get install libmodule-build-perl
echo "" | apt-get install libsoap-wsdl-perl

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
## To free space, remove apt-get packages that are no longer required.
apt-get autoremove

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

emacs -nw ${RSAT}/public_html/web_services/RSATWS.wsdl

## At the bottom of the file, locate the following line.
##  <soap:address location="http://rsat.ulb.ac.be/rsat/web_services/RSATWS.cgi"/>

## Adapt the URL to your local configuration.

## After this, you should re-generate the web services stubb, with the
## following command.
cd $RSAT
make -f makefiles/init_rsat.mk ws_stubb

## Test the local web services
make -f makefiles/init_rsat.mk ws_stubb_test

## Test RSAT Web services (local and remote) without using the SOAP/WSDL stubb
## (direct parsing of the remote WSDL file)
make -f makefiles/init_rsat.mk ws_nostubb_test

################################################################
## R installation

## As sudo, I edited the file /etc/apt/sources.list 
## and added the following line 
## (see instructions on http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu/)
##   deb http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu trusty/
## I then updated the apt-get packages
sudo apt-get update

## .. and installed the R base package
sudo apt-get -y install r-base
sudo apt-get install -y r-base-dev

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
sudo R

## At the R prompt
install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras"))
source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
quit()


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
