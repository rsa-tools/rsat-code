
################################################################
##########   Ubuntu-specific installation
################################################################

################################################################
## install some modules required for proper functioning of RSAT o
## Ubuntu 12
ubuntu12_addons:
	sudo apt-get install zlib1g-dev
	sudo apt-get install libgd2-xpm-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install libmysqlclient15-dev
	sudo apt-get install libdb-dev
	sudo apt-get install libberkeleydb-perl
	sudo apt-get install ia32-libs

## Note Berkeley db package should be isntalled with synaptic
ubuntu12_install_logwatch:
	sudo apt-get install logwatch
	sudo logwatch.conf

################################################################
## Install tomcat7 for Ubuntu server
ubuntu12_tomcat:
	sudo apt-get install libtomcat7-java
	sudo apt-get install tomcat7-common
	sudo apt-get install tomcat7
	sudo apt-get install tomcat7-docs

################################################################
## Packages specifically required for Ubuntu 14
install_ubuntu14_packages:
	sudo apt-get install python
	sudo apt-get install python-pip
	sudo apt-get install python-dev
	sudo apt-get install ipython
	sudo apt-get install ipython-notebook
	sudo apt-get install git
	sudo apt-get install emacs
	sudo apt-get install make
	sudo apt-get install g++
	sudo apt-get install yum
	sudo apt-get install apache2
	sudo apt-get install php5
	sudo apt-get install libapache2-mod-php5
	sudo apt-get install php-elisp
	sudo apt-get install perldoc
	sudo apt-get install texlive-latex-base
	sudo apt-get install cvs
	sudo apt-get install wget
	sudo apt-get install finger
	sudo apt-get install zip
	sudo apt-get install unzip
	sudo apt-get install libgd2-xpm-dev
	sudo apt-get install libgd-gd2-perl
	sudo apt-get install libxml2-dev
	sudo apt-get install libnet-ssleay-perl
	sudo apt-get install libcrypt-ssleay-perl
	sudo apt-get install libssl-dev
	sudo apt-get install ghostscript
	sudo apt-get install gnuplot
	sudo apt-get install graphviz
	sudo apt-get install links
	sudo apt-get install lib32z1 lib32ncurses5 lib32bz2-1.0

## Some linux packages required for R BioConductor
	sudo apt-get install -y make
	sudo apt-get install libc6-dev
	sudo apt-get install gfortran
	sudo apt-get install build-essential
	sudo apt-get install libreadline-gplv2-dev:i386 lib64readline-gplv2-dev:i386 libreadline-gplv2-dev
	sudo apt-get install libx11-dev
	sudo apt-get install libxt-dev
	sudo apt-get install libcurl4-openssl-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install texlive-full
	sudo apt-get install tcl8.5-dev
	sudo apt-get install tk8.5-dev
	sudo apt-get install libxss-dev
	sudo apt-get install libpng12-dev
	sudo apt-get install libjpeg62-dev
	sudo apt-get install libcairo2-dev

## Unable to locate package. WHy do I have it in the list ? Was this
## package reomved after Ubuntu 12 ?
# 	sudo apt-get install apache2-util

## The following packages "have no candidates". Maybe they were required for Ubuntun 12 ?
##	gfortran-4.3 \
##	gcj
##
##  .. to be completed

################################################################
## Python modules
	sudo apt-get install python python-pip python-dev
	sudo pip install numpy
	sudo pip install scipy
	sudo pip install matplotlib
	sudo pip install soappy

## We need both python2.7 and python3 (for different scripts)
	sudo apt-get install python3 python3-pip 
	sudo pip3 install numpy
	sudo pip3 install scipy
	sudo pip3 install matplotlib
	sudo pip3 install soappy

## This line does not work, although I found it in some docs
## sudo apt-get install pyton3-dev
