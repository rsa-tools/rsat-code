################################################################
## Instructions used to install a Virtual Machine on the IFB cloud
## (Institut Francais de Bioinformatique), and on a VirtualBox VM.

################################################################
## Althu-ough the last stable version of Ubuntu is 14.04, I test a
## previous release (12.04) because with 14.04 I have problemw with
## the Perl SOAP/WSDL library, whcih is indispensable to RSAT.

## THIS IS NOT REALLY A bahsrc FILE, IT IS A SUCCESSION OF
## INSTRUCTIONS AND COMMENTS, THAT SHOULD BE DONE MANUALLY. I (JvH)
## SHOULD IMPROVE THIS WHEN I CAN.


################################################################
## Preparing to install Ubuntu
##   x Download updates while installing
##   x Install third-party software
##   o Erase disk and install Ubuntu
## Keyboard layout
##   French > French (Macintosh)
## Who are you
##  Your name: Regulatory Sequence Analysis Tools admin
##  Computers name: RSAT-Ub12.04
##  Username: rsat
##  o Require my password to log in

## TO RUN Virtualbox on Mac OSX
#    Go into your Ubuntu Settings
#    Navigate to Keyboard → Keyboard Layout Settings
#    add English (Macintosh)


################################################################
## Must be executed as root. If you are non-root but sudoer user, you
## can become it withn "sudo bash"

## Aptitude is more convenient than apt-get to treat dependencies when
## installing and uninstalling packages.
apt-get --quiet --assume-yes install aptitude

## I install openssh, in order to do the rest via an external ssh
## connection (more convenient for keyboard shortcuts, copy-paste etc).
aptitude --quiet --assume-yes install openssh-server

## Done network specification
reboot

## I run the next installation steps from the host terminal, because I
## feel more confortable in my usual environment (in particular the keyboard).
ssh rsat@192.168.56.112
sudo bash
## Enter rsat password

## We need to update apt-ge, to avoid trouble with python
## See http://askubuntu.com/questions/350312/i-am-not-able-to-install-easy-install-in-my-ubuntu
aptitude update
aptitude --quiet --assume-yes upgrade ## This takes a while

APTGET_LIBRARIES="emacs23
        ia32-libs
	ssh
	git
	cvs
	wget
	curl
	g++
	make
	zip
	unzip
	finger
	screen
	yum
	perl-doc
	pmtools
	apache2
	php5
	libapache2-mod-php5
	php-elisp
	lib32z1
	lib32ncurses5
	lib32bz2-1.0
	libgdbm-dev
	libgd-tools
	libgd-gd2-perl
	libgd2-xpm-dev
	libxml2-dev
	libnet-ssleay-perl
	libcrypt-ssleay-perl
	libssl-dev
	ghostscript
	gnuplot
	graphviz
	links
	libc6-dev
	gfortran
	build-essential
	lib64readline-gplv2-dev
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
	mysql-client
	libmysqlclient-dev
	default-jre
	texlive-latex-base
	python
	python-setuptools
	python-virtualenv
	python-pip
	python-dev
	python-suds
	ipython
	ipython-notebook
	python-numpy
	python-scipy
	python-matplotlib
	python-soappy
	python-pygraphviz
	python3
	python3-setuptools
	python3-dev
	python3-numpy
	python3-scipy
	libmodule-build-perl
	libsoap-wsdl-perl
	libsoap-lite-perl
"


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
	libio-all-perl \
	libobject-insideout-perl \
	libobject-insideout-perl \
	libsoap-lite-perl \
	libsoap-wsdl-perl \
	libxml-perl \
	libxml-simple-perl \
	libdbi-perl \
	liblockfile-simple-perl \
	libobject-insideout-perl \
	libgd-perl \
	libdbd-mysql-perl \
	libjson-perl \
	libbio-perl-perl \
	libdigest-md5-file-perl
"

## The following modules exist in Ubuntu 14.04 but not in Ubuntu 12.04
##	digest-md5-file-perl \
## 	liblockfile-simple \
##	libutil-properties-perl \
##	libxml-compile-cache-perl \


## This library is not found in Ubuntun 12.04
##     libgvc6
## I should ceck if it is OK (required for graphviz)

## Install the apt-get libraries
mkdir -p install_logs
echo "Installing apt-get libraries"
for LIB in $APTGET_LIBRARIES $APTGET_PERLMOD
do
   echo "`date`        installing apt-get library $LIB"
   aptitude --quiet --assume-yes install ${LIB} > install_logs/aptitude_install_${LIB}.txt
   df -h > install_logs/aptitude_install_${LIB}_df.txt
done
echo "Log files are in folder install_logs"


# ## 32-bit compatibility libraries are required for some packages
# aptitude --quiet --assume-yes install ia32-libs

# ## Enable incoming ssh
# aptitude --quiet --assume-yes install ssh

# ## Concurrent versioning systems
# aptitude --quiet --assume-yes install git
# aptitude --quiet --assume-yes install cvs

# ## Web aspirators
# aptitude --quiet --assume-yes install wget
# aptitude --quiet --assume-yes install curl

# ## Compilation tools
# aptitude --quiet --assume-yes install g++
# aptitude --quiet --assume-yes install make

# ## Utilities
# aptitude --quiet --assume-yes install zip
# aptitude --quiet --assume-yes install unzip
# aptitude --quiet --assume-yes install finger
# aptitude --quiet --assume-yes install screen
# aptitude --quiet --assume-yes install yum

# ## Perl packages
# aptitude --quiet --assume-yes install perl-doc
# aptitude --quiet --assume-yes install pmtools

# ## Apache and utilities
# aptitude --quiet --assume-yes install apache2
# aptitude --quiet --assume-yes install php5
# aptitude --quiet --assume-yes install libapache2-mod-php5
# aptitude --quiet --assume-yes install php-elisp

# ## Graphic libraries
# aptitude --quiet --assume-yes install lib32z1
# aptitude --quiet --assume-yes install lib32ncurses5
# aptitude --quiet --assume-yes install lib32bz2-1.0

# aptitude --quiet --assume-yes install libgdbm-dev
# aptitude --quiet --assume-yes install libgd-tools

# aptitude --quiet --assume-yes install libgd-gd2-perl
# aptitude --quiet --assume-yes install libgd2-xpm-dev
# aptitude --quiet --assume-yes install libxml2-dev

# aptitude --quiet --assume-yes install libnet-ssleay-perl
# aptitude --quiet --assume-yes install libcrypt-ssleay-perl
# aptitude --quiet --assume-yes install libssl-dev

# ## Graphic libraries and software tools
# aptitude --quiet --assume-yes install ghostscript
# aptitude --quiet --assume-yes install gnuplot
# aptitude --quiet --assume-yes install graphviz

# ## Text-mode Web browser, used by some packages
# aptitude --quiet --assume-yes install links

# ## Some linux packages required for R BioConductor
# aptitude --quiet --assume-yes install libc6-dev
# aptitude --quiet --assume-yes install gfortran
# aptitude --quiet --assume-yes install build-essential

# aptitude --quiet --assume-yes install libreadline-gplv2-dev:i386
# aptitude --quiet --assume-yes install lib64readline-gplv2-dev:i386
# aptitude --quiet --assume-yes install libreadline-gplv2-dev

# aptitude --quiet --assume-yes install libx11-dev
# aptitude --quiet --assume-yes install libxt-dev
# aptitude --quiet --assume-yes install libcurl4-openssl-dev

# aptitude --quiet --assume-yes install libxml2-dev
# ## BEWARE: texlive-full occupies a lot of disk space. I should check if this is really required (for R ?)
# ## aptitude install texlive-full
# aptitude --quiet --assume-yes install tcl8.5-dev

# aptitude --quiet --assume-yes install tk8.5-dev
# aptitude --quiet --assume-yes install libxss-dev
# aptitude --quiet --assume-yes install libpng12-dev

# aptitude --quiet --assume-yes install libjpeg62-dev
# aptitude --quiet --assume-yes install libcairo2-dev

# ## mysql client is required for ensembl client scripts
# aptitude --quiet --assume-yes install mysql-client
# aptitude --quiet --assume-yes install libmysqlclient-dev

# ## Java
# ## seems to be required for SOAP::WSDL Perl module
# aptitude --quiet --assume-yes install default-jre
# ## aptitude --quiet --assume-yes install default-jdk

# ## Latex is required for RSAT doc + other applications (e.g. R). Note
# ## that it takes a some time to install
# aptitude --quiet --assume-yes install texlive-latex-base


# ################################################################
# ## Python and modules
# aptitude --quiet --assume-yes install python
# aptitude --quiet --assume-yes install python-setuptools
# aptitude --quiet --assume-yes install python-virtualenv
# aptitude --quiet --assume-yes install python-pip
# aptitude --quiet --assume-yes install python-dev
# aptitude --quiet --assume-yes install python-suds

# ## iPython
# aptitude --quiet --assume-yes install ipython
# aptitude --quiet --assume-yes install ipython-notebook


## Problem: "pip install matplotlib ;" ## Does not work. matplotlib
## can be installed with easy_install, but for Ubuntu it is probably
## better to use apt-get. I install what I can with
##
## A fix for a problem to install scipy (does not work with pip): in
## Ubuntu, we can use apt-get build-dep taken from here:
## http://stackoverflow.com/questions/11863775/python-scipy-install-on-ubuntu
# apt-get --quiet --assume-yes build-dep python-numpy
# apt-get --quiet --assume-yes build-dep python-scipy
# apt-get --quiet --assume-yes build-dep python-matplotlib
# apt-get --quiet --assume-yes build-dep python-soappy ; ## For web services
# apt-get --quiet --assume-yes install libgvc6 ; ## Required for pygraphviz
# apt-get --quiet --assume-yes build-dep python-pygraphviz

## On other systems, the install may work with pip (to be
## checked). Alternatively, can the python librarires can be installed
## wit easy_install (commented lines below).
# easy_install -U distribute numpy
# easy_install -U distribute scipy
# easy_install -U distribute matplotlib
#easy_install -U distribute soappy
#easy_install -U distribute pygraphviz

## For some python packages there is no apt-get package
easy_install -U distribute fisher

################
## We need both python2.7 and python3 (for different scripts)
##
## Attention, I use "build-dep" rather than "install".
##
## aptitude --quiet --assume-yes build-dep python3-pip ## Did no exist yet in Ubuntu 12.04, has to be installed with easy_install (see below)
aptitude --quiet --assume-yes build-dep python3-setuptools
aptitude --quiet --assume-yes build-dep python3-dev
aptitude --quiet --assume-yes build-dep python3-numpy
aptitude --quiet --assume-yes build-dep python3-scipy
aptitude --quiet --assume-yes build-dep python3-matplotlib ; ## Picking 'matplotlib' as source package instead of 'python3-matplotlib'
## aptitude --quiet --assume-yes build-dep python3-pygraphviz ; ## E: Unable to find a source package for python3-pygraphviz
## aptitude --quiet --assume-yes install python3-soappy ; ## E: Unable to locate package python3-soappy
## aptitude --quiet --assume-yes install python3-suds ; ## E: Unable to locate package python3-suds

## This differs from Ubuntun 14.04: pip3 must be installed with easy_install ?
easy_install3 pip

## Problem : No distributions at all found for python-suds
## pip3 install python-suds
## easy_install -U distribute python-suds ; ## error: Could not find suitable distribution for Requirement.parse('python-suds')

easy_install3 -U distribute pygraphviz ## SyntaxError: invalid syntax
easy_install3 -U distribute fisher

## Problems:
# pip3 install wsdl
# pip3 install wstools
# pip3 install fisher
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
##     deb http://us.archive.ubuntu.com/ubuntu trusty main universe
## You can now quit emacs

apt-get update

# apt-get --quiet --assume-yes install libmodule-build-perl
# apt-get --quiet --assume-yes install libsoap-wsdl-perl
# apt-get --quiet --assume-yes install libsoap-lite-perl

## Note: this is still not sufficient to get SOAP::WSDL to run the two
## following targets
##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stub
##     make -f ${RSAT}/makefiles/init_rsat.mk ws_stub_test

## We first need to fix some problem with CPAN : on Ubuntu, I cannot
## install the SOAP::WSDL module, which is required for several
## functionalities. More precisely, after fiddling around for a few
## hours, the server is able to answer to web services requests, but I
## cannot run clients on it. Since NeAT relies on WSDL clients, it is
## impossible to have neat running on and Ubuntu server. The
## installation is however possible, since the stub can be generated
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
## To free space, remove apt-get packages that are no longer required,
## and clean the installation cache.
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
cd /home/rsat/rsat
source RSAT_config.bashrc
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
cat check_perl_modules_eval.txt

## Install these modules Beware: the _noprompt suffix is optional. It
## has the advantage to avoid for the admin to confirm each
## installation step, but the counterpart is that errors may be
## overlooked.
make -f makefiles/install_rsat.mk perl_modules_install_noprompt

## Note: I had to force installation for the some modules, because
## there seem to be some circular dependencies.
#make -f makefiles/install_rsat.mk perl_modules_install_by_force

## Check if all required Perl modules have been correctly installed
make -f makefiles/install_rsat.mk perl_modules_check
## Check if all modules are OK
cat check_perl_modules_eval.txt

## I still have problems with InsideOut, SOAP::Transport


################################################################
## Activate the Apache Web server and RSAT configuration


## sudo emacs -nw /etc/apache2/sites-available/000-default.conf
## Activate the following line:
# Include conf-available/serve-cgi-bin.conf
emacs -nw /etc/apache2/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi

## From http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/
# sudo chmod 755 /usr/lib/cgi-bin
# sudo chown root.root /usr/lib/cgi-bin
# sudo a2enmod cgi ## this is apparently required to enable cgi
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

## Optionally, install some pluricellular model organisms (requires
## more space)
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

## After this, you should re-generate the web services stub, with the
## following command.
make -f makefiles/init_rsat.mk ws_stub

## Test the local web services
make -f makefiles/init_rsat.mk ws_stub_test

## Test RSAT Web services (local and remote) without using the SOAP/WSDL stub
## (direct parsing of the remote WSDL file)
make -f makefiles/init_rsat.mk ws_nostub_test

################################################################
## R installation

## As sudo, I edited the file /etc/apt/sources.list
## and added the following line
## (see instructions on http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu/)
##   deb http://mirror.ibcp.fr/pub/CRAN/bin/linux/ubuntu trusty/
## I then updated the apt-get packages
sudo apt-get update

## .. and installed the R base package
sudo apt-get --quiet --assume-yes install r-base
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



################################################################
## TO DO

## test SOAPUI
##   http://www.upubuntu.com/2012/04/how-to-install-soapui-web-service.html
## GUI to test SOAP WS
