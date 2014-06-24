################################################################
## Instructions used to install a Virtual Machine on the IFB cloud
## (Institut Francais de Bioinformatique).
##
##

## Must be executed as root
sudo bash

## We need to update apt-ge, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu
echo "" | apt-get update
echo "" | apt-get upgrade
echo "" | apt-get install python-setuptools 
echo "" | apt-get install python
echo "" | apt-get install python-virtualenv
echo "" | apt-get install python-pip
echo "" | apt-get install python-dev

echo "" | apt-get install python3
echo "" | apt-get install python3-pip

echo "" | apt-get install ipython
echo "" | apt-get install ipython-notebook

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

## Cannot survive without emacs
echo "" | apt-get install emacs

echo "" | apt-get install make
echo "" | apt-get install g++
echo "" | apt-get install yum

## Problem : unable to locate package perldoc
## apt-get install perldoc 

## Apache and utilities
echo "" | apt-get install apache2
echo "" | apt-get install php5
echo "" | apt-get install libapache2-mod-php5
echo "" | apt-get install php-elisp

## Latex is required for other packages. Note that it takes a some time to install
echo "" | apt-get install texlive-latex-base

## Graphic libraries and software tools
echo "" | apt-get install libgd2-xpm-dev
echo "" | apt-get install libgd-gd2-perl
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

## Text-mode Web browser
echo "" | apt-get install links

## Some linux packages required for R BioConductor
echo "" | apt-get install -y make
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

################################################################
## Python modules

echo "" | apt-get install python
echo "" | apt-get install python-pip
echo "" | apt-get install python-dev
echo "" | apt-get install python-virtualenv
echo "" | apt-get install python-pip
echo "" | apt-get install python-suds

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
pip3 install python-suds
pip3 install wsdl
pip3 install wstools
pip3 install fisher
## pip3 install pygraphviz ## This fails ! Command python setup.py egg_info failed with error code 1 in /tmp/pip_build_root/pygraphviz

## soappy seems to be discontnued for python3 !
# pip3 install soappy
## I should test one of the following SOAP packages
pip3 install suds-jurko
pip3 install pysimplesoap

################################################################
################       RSAT installation        ################
################################################################

## Create a specific user for RSAT. The user is named rsat
sudo adduser rsat
## Name: Regulatory Sequence Analysis Tools user

## Grant sudo privileges to the rsat user (will be more convenient for
## installing Perl modules, software tools, etc)
visudo
## then add the following line below "User privilege specification"
# rsat    ALL=(ALL:ALL) ALL


## The installation is done under the rsat login
su - rsat

## Get the RSAT package
git clone git@depot.biologie.ens.fr:rsat

## Run the configuration script, to specify the environment variables.
cd rsat
perl perl-scripts/configure_rsat.pl 


## Come back to the rsat identify
source RSAT_config.bash

## Initialise RSAT folders
make -f makefiles/init_rsat.mk init

################################################################
## Link the RSAT bash configuration file to a directory where files
## are loaded by each user at each login. Each user will then
## automatically load the RSAT configuration file when opening a bash
## session.
sudo (cd /home/rsat/rsat/;  rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/)

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

## The following commands should be executed with admin rights.
sudo bash

## Check that RSAT path has been defined
echo $RSAT

## Set working directory to RSAT
cd $RSAT

## Install perl-doc package
echo "" | apt-get install perl-doc

## Get the list of Perl modules to be installed
make -f makefiles/install_rsat.mk  perl_modules_list

## Check which Perl modules are already installed
make -f makefiles/install_rsat.mk perl_modules_check
## The locations of installed modules are stored in perl_modules_check.txt

## Install these modules Beware: the _noprompt suffix is optional. It
## has the advantage to avoid for the admin to confirm each
## installation step, but the counterpart is that errors may be
## overlooked.
make -f makefiles/install_rsat.mk perl_modules_install_noprompt

## Check if all required Perl modules have been correctly installed
make -f makefiles/install_rsat.mk perl_modules_check

## Some modules are not installed with the installation procedure,
## due to problems with their cpan declaration.

## Exit from the root shell, and become rsat user again
exit




## compile RSAT programs written in C
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all


## Install some third-party programs required by some RSAT scripts.
make -f makefiles/install_software.mk install_ext_apps

## Install two model organisms, required for some of the Web tools.
download-organism -v 1 -org Saccharomyces_cerevisiae
download-organism -v 1 -org Escherichia_coli_K_12_substr__MG1655_uid57779


################################################################
## At this stage you can already check some simple RSAT command 

## Test a simple Perl script that does not require for organisms to be
## installed.
which random-seq
random-seq -l 100

## Test a simple python script that does not require organisms to be
## installed.
random-motif -l 10 -c 0.90

## Test some external programs

## get the help for seqlogo
which seqlogo
seqlogo

## ghostscript
which gs
gs --version

## Get the list of organisms supported on your computer.
supported-organisms


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
sudo a2enmod cgi ## this is supposed to enable cgi , but apparently not sufficient
apache2ctl restart


################################################################
## Configure the Web services


emacs -nw ${RSAT}/public_html/web_services/RSATWS.wsdl

## At the bottom of the file, locate the following line.
##  <soap:address location="http://rsat.ulb.ac.be/rsat/web_services/RSATWS.cgi"/>

Adapt the URL to your local configuration.

After this, you should re-generate the web services stubb, with the
following command.

  \begin{lstlisting}
cd $RSAT; 
make -f makefiles/init_rsat.mk ws_stubb
  \end{lstlisting}






################################################################
## R installation

## I edited the file /etc/apt/sources.list 
## and added the following line 
## (see instructions on http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu/)
##   deb http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu trusty

## I then updated the apt-get packages
sudo apt-get update

## .. and installed the R base package
sudo apt-get install r-base
sudo apt-get install r-base-dev

## Installation of R packages

## The command R CMD INSTALL apparently does not work at this stage.
##	root@rsat-tagc:/workspace/rsat# R CMD INSTALL reshape
##	Warning: invalid package 'reshape'
##	Error: ERROR: no packages specified

cd $RSAT; make -f makefiles/install_rsat.mk  r_modules_list 

## In install them from the R interface
sudo R

## At the R prompt
install.packages(c("reshape", "RJSONIO", "plyr", "dendroextras"))
source('http://bioconductor.org/biocLite.R'); biocLite("ctc")
quit()



################################################################
##
## ATTENTION: to ensure persistence, you imperatively have to run the
## following command in the VM before any shutdown onthe Web site
##
################################################################

cd; bash cleaner.sh ; history -c && history -w && logout
