################################################################
## Site-specific configuration for the RSAT makefiles
## 
## This file is automatically loade by the script
## ${RSAT}/makefiles/util.mk, which is itself loaded by all the RSAT
## makefile scripts.

## Name of the site (which will appear in the log file)
RSAT_SITE=your_server_name

## Web server URL
RSAT_WWW=http://localhost/rsat/

## Email address of RSAT administrator
RSAT_ADMIN_EMAIL=your.mail@your.mail.server

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
## Options:
##   PACKAGE_MANAGER=brew install
##   PACKAGE_MANAGER=yum install
##   PACKAGE_MANAGER=get-apt
PACKAGE_MANAGER=

#UCSC_OS=macOSX.i386
UCSC_OS=linux.x86_64

## Directory to store the downloaded software (sources, before
## compilation)
SRC_DIR=${RSAT}/app_sources

## Directory to store libraries
PERLLIB_DIR=${RSAT}/lib

################################################################
## Your sudoer status
##
## The sudo command is required for installing some third-party
## software in directories accessible to all users (e.g. /usr/bin).
##
## Supported: sudo (or empty if you are not administrator)
##
## If you are not sudoer of your server, you should let this variable
## empty, and make sure that the BIN_DIR defined below is set toa
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
BIN_DIR=${RSAT}/bin
#BIN_DIR=/usr/local/bin

################################################################
## Web services
##
## This may differ from the RSAT site, because some tasks may be
## delegated to remote RSAT web services.
RSAT_WS=http://localhost/rsat

## URL used by the RSAT Web services to give Web access to temporary
## files
RSAT_WS_TMP=http://localhost/rsat/tmp

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
## supported: torque | sge
QUEUE_MANAGER=torque

## Name of the queue where the jobs have to be sent
QUEUE=

