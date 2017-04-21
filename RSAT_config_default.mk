################################################################
## Site-specific configuration for the RSAT makefiles
## 
## This file is automatically loade by the script
## ${RSAT}/makefiles/util.mk, which is itself loaded by all the RSAT
## makefile scripts.

## Name of the site (which will appear in the log file)
RSAT_SITE=your_server_name

## Group specificity of this server
GROUP=None

## Default species for genome installation, which should be adapted to the group specificity
SPECIES=Saccharomyces_cerevisiae

## Web server URL
# RSAT_WWW=http://localhost/rsat/

## Email address of RSAT administrator
RSAT_SERVER_ADMIN=your.mail@your.mail.server

## Operating system
## Supported: linux  | macosx
##
## This information is required for installing some pre-compiled
## software tools.
OS=linux
ARCHITECTURE=x64

################################################################
## Manager for installing Unix packages
##
## Examples:
##   PACKAGE_MANAGER=brew
##   PACKAGE_MANAGER=yum
##   PACKAGE_MANAGER=apt-get
##   PACKAGE_MANAGER=aptitude
PACKAGE_MANAGER=apt-get


################################################################
## Operating system to download UCSC userApps from
## http://hgdownload.cse.ucsc.edu/admin/exe/
##
## Supported 
##  UCSC_OS=macOSX.i386
##  UCSC_OS=linux.x86_64
UCSC_OS=linux.x86_64

## Directory to store the downloaded software (sources, before
## compilation)
SRC_DIR=${RSAT}/app_sources

################################################################
## Your sudoer status
##
## The sudo command is required for installing some third-party
## software in directories accessible to all users (e.g. /usr/bin).
##
## Supported: sudo (or empty if you are not administrator)
##
## If you are not sudoer of your server, you should let this variable
## empty, and make sure that the RSAT_BIN defined below is set toa
## directory of your ownership. 
##
## Your non-admin status might cause problems for the "install" step
## of some third-party software, but it will at least allow you to
## have a functional RSAT disribution locally.
##
#SUDO=sudo
SUDO=

## Default installation dir for binaries. By default, will beinstalled
## in /usr/local/bin so they will be accessible to all users. This
## however requires admin rights. If you don't dispose of admin
## rights, an alternative is to set this directory in the RSAT fodler.
##
## Typical values:
##   RSAT_BIN=/usr/local/bin
##   RSAT_BIN=\${RSAT}/bin
RSAT_BIN=${RSAT}/bin


################################################################
## Web services
##
## This variable defines the address of the Web services, which will
## be used to generate a stub. This stub will be used (under others)
## by the Web interface of retrieve-ensembl-seq, and by the NeAT
## tools. It is thus important to specify it correctly and test it.
##
## By default, the web services are running on the host machine
## itself. We thus set it to http://localhost/rsat/. This means that
## the Web site and Web services will run on the same
## machine. However, in some cases it might be useful to use a static
## IP address. In principle, this would allow to redirect the Web
## services towards a separate RSAT server. However, this is somewhat
## tricky and might have side effects that we did not test. 
##
## We recommend to leave this value to its default (localhost).
RSAT_WS=http://localhost/rsat/

################################################################
##
## Configuration of a PC cluster
##
## The parameter QUEUE_MANAGER is required for running tasks on a PC
## cluster. The program used to send jobs to the queue (qsub) has
## different parameters depending on the queue manager. Thus, a qsub
## command for torque will not be understood by sge and
## reciprocally. To circumvent this problem, we specify the queue
## manager in this property file, so RSAT can adapt the options
## consequently.
##
## supported: batch | torque | sge
QUEUE_MANAGER=batch

## Name of the queue where the jobs have to be sent
CLUSTER_QUEUE=

################################################################
## Ensembl and EnsemblGenomes releases. 
##
## These variables should in principle be updated according to the
## relases of the two databases:
##
##  Ensembl: http://www.ensembl.org/
##  EnsemblGenomes: http://ensemblgenomes.org/
ENSEMBL_RELEASE=88
ENSEMBLGENOMES_RELEASE=35
