#!/usr/bin/env bash

################################################################
## This script installs RSAT on a Linux Debian operating system

## Author: François-Xavier Theodule and Jacques van Helden <Jacques.van-Helden@univ-amu.fr>

## Update et upgrade apt-get
apt-get update
apt-get -y upgrade
apt install  -y openssh-client

## Time setup
echo Europe/Paris > /etc/timezone
dpkg-reconfigure -f noninteractive tzdata

## Install debian packages
apt install -y git cvs r-base-core net-tools python python-numpy libssl-dev libexpat1-dev libxml2-dev apache2 apt-utils libapache2-mod-perl2 libapache2-mod-perl2-dev php5 perl-doc libdb5.3-dev cpanminus wget sudo

## Install perl modules
cpanm --force Net::SSLeay
##cpanm --force Net::SMTPSXML::Compile::SOAP11 ## --> could not find JvH PROBLEM: we need SOAP11 for Web services, to be revised
cpanm --force Email::Sender::Transport::SMTPS
cpanm --force XML::Parser::Expat
cpanm --force SOAP::WSDL
cpanm --force XML::LibXML
cpanm --force XML::LibXML::Simple
cpanm --force XML::Compile
cpanm --force XML::Compile::Cache
cpanm --force XML::Compile::SOAP11
cpanm --force XML::Compile::WSDL11
cpanm --force YAML
cpanm --force Module::Build::Compat
cpanm --force CGI
cpanm --force Email::Sender
cpanm --force Email::Simple
cpanm --force Email::Simple::Creator
cpanm --force PostScript::Simple
cpanm --force Statistics::Distributions
cpanm --force Math::CDF
cpanm --force Algorithm::Cluster
cpanm --force File::Spec
cpanm --force POSIX
cpanm --force Data::Dumper
cpanm --force Digest::MD5::File
cpanm --force IO::All
cpanm --force LockFile::Simple
cpanm --force Object::InsideOut
cpanm --force Util::Properties
cpanm --force Class::Std::Fast
cpanm --force GD
cpanm --force DBI
cpanm --force DBD::mysql
cpanm --force DB_File
cpanm --force LWP::Simple
cpanm --force REST::Client
cpanm --force JSON
cpanm --force HTTP::Tiny
cpanm --force XML::Compile::Transport::SOAPHTTP
cpanm --force SOAP::Lite
cpanm --force SOAP::Packager
cpanm --force SOAP::Transport::HTTP
cpanm --force SOAP::WSDL
cpanm --force Bio::Perl
cpanm --force Bio::Das
cpanm --force XML::DOM
cpanm --force Spreadsheet::WriteExcel::Big
cpanm --force Spreadsheet::WriteExcel
cpanm --force Log::Log4perl
cpanm --force Number::Format
cpanm --force OLE::Storage_Lite
cpanm --force Template::Plugin::Number::Format
cpanm --force Readonly

## Validate Perl cgi et php under apache2
a2enmod cgi
a2enmod perl
a2enmod php5

## Download and install rsat tar.gz in /rsat
cd /
wget --no-clobber http://pedagogix-tagc.univ-mrs.fr/download_rsat/rsat_2016-05-24.tar.gz
mv rsat_2016-05-24.tar.gz rsat.tar.gz
gunzip rsat.tar.gz
tar -xvf rsat.tar

## RSAT configuration and installation
##
## JvH: THIS SHOULD BE REVISED using configure-rsat.pl, but after I
## add an automatic parameter specification

## On ne va pas jour le scriptperl perl-scripts/configure_rsat.pl on modifie les ficheirs directement
## Trois fichiers concernes :
##  RSAT_config.props
##    config file read by RSAT programs in various languages: Perl,
##    python, java
## cp  RSAT_config_default.props  RSAT_config.props
##
## RSAT_config.mk
##    environment variables for the makefiles
##cp RSAT_config_default.mk RSAT_config.mk
##
## RSAT_config.conf
##    RSAT configuration for the Apache web server.
## cp  RSAT_config_default.conf  RSAT_config.conf
## mise a jour RSAT_config.bashrc
## path /rsat

## Mise a jour fichier RSAT_config.bashrc a partir de RSAT_config_default.bashrc
cd /rsat
sed 's/\[RSAT_PARENT_PATH\]//g' ./RSAT_config_default.bashrc > ./RSAT_config.bashrc
echo -e "export RSAT_PARENT_PATH=\n$(cat  RSAT_config.bashrc)" > RSAT_config.bashrc
source RSAT_config.bashrc
echo "source /rsat/RSAT_config.bashrc" >>  /etc/bash.bashrc
echo "service apache2 start" >>  /etc/bash.bashrc 

## Update RSAT_config.props
i=`ifconfig eth0 | awk '/inet adr:/{print $2}' | awk -F ':' '{print $2}'` | sed -e "s/your_server_name/$i/g"  ./RSAT_config_default.props | sed 's/your.mail@your.mail.server/postmaster@rsat.com/g' | sed 's/\[RSAT_PARENT_PATH\]//g' > ./RSAT_config.props

## Update RSAT_config.mk apt-get for debian
## RSAT_config.mk apt-get
i=`ifconfig eth0 | awk '/inet adr:/{print $2}' | awk -F ':' '{print $2}'` | sed -e "s/your_server_name/$i/g"  ./RSAT_config_default.mk | sed 's/your.mail@your.mail.server/postmaster@rsat.com/g' | sed 's/\[RSAT_PARENT_PATH\]//g' | sed 's/SUDO=/SUDO=sudo/g' > ./RSAT_config.mk

## Update RSAT_config.mk YUM for REDHAT
## RSAT_config.mk yum
## i=`ifconfig eth0 | awk '/inet adr:/{print $2}' | awk -F ':' '{print $2}'` | sed -e "s/your_server_name/$i/g"  ./RSAT_config_default.mk | sed 's/your.mail@your.mail.server/postmaster@rsat.com/g' | sed 's/\[RSAT_PARENT_PATH\]//g' | sed 's/SUDO=/SUDO=sudo/g' | sed 's/apt-get/yum/g' > ./RSAT_config.mk

## Make directories
mkdir /rsat/public_html/logs
mkdir /DockerIn
mkdir /DockerOut
mkdir /DockerTodo
make -f makefiles/init_rsat.mk init

## Compile
make -f makefiles/init_rsat.mk compile_all

##vmatch licence qui est dans /somewhere/vmatch_RSATVM-IFB_2015-07-06.lic
cp /RsatInstall/vmatch_RSATVM-IFB_*.lic /rsat/bin/vmatch.lic

## Install thrid party

################ JvH: to revise. ################
## PROBLEME install blast et version de vmatch (2.2.4-->2.2.5)
mv makefiles/install_software.mk makefiles/install_software_default.mk
sed -e 's/install_blast//g' makefiles/install_software_default.mk | sed -e 's/2.2.4/2.2.5/g' > makefiles/install_software.mk
make -f makefiles/install_software.mk install_ext_apps

## Install blast
cd /rsat/bin
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz
gunzip ncbi-blast-2.4.0+-x64-linux.tar.gz
tar -xvf ncbi-blast-2.4.0+-x64-linux.tar
mv ncbi-blast-2.4.0+/bin/* .
rm ncbi-blast-2.4.0+-x64-linux.tar
rm -rf ncbi-blast-2.4.0+-
cd /rsat

## Install two small-genome model organisms
download-organism -v 1 -org Saccharomyces_cerevisiae
download-organism -v 1 -org Escherichia_coli_K_12_substr__MG1655_uid57779

## Modify Apache2 configuration from the template file  /rsat/RSAT_config_default.conf
sed -e 's/\[RSAT_PARENT_PATH\]//g' RSAT_config_default.conf > /etc/apache2/sites-enabled/rsat.conf

## Adapt Apache2 config for Perl execution
## affected files: /etc/apache2/apache2.conf RSAT_config_default.conf (/etc/apache2/sites-enabled/rsat.conf) 
echo 'AddHandler cgi-script .cgi .pl' >> /etc/apache2/apache2.conf

## Restart apache2 service
service apache2 restart

