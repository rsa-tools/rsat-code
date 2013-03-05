################################################################
## Site-specific configuration for the RSAT makefiles
## 
## This file is automatically loade by the script
## ${RSAT}/makefiles/util.mk, which is itself loaded by all the RSAT
## makefile scripts.

## Name of the site (which will appear in the log file)
RSAT_SITE=your_server_name


## Email address of RSAT administrator
RSAT_ADMIN_EMAIL=your.address@your.domain

## Operating system
## Supported: linux  | macosx
OS=linux
ARCHITECTURE=x64

## sudo command is required for installing some software
SUDO=sudo

#UCSC_OS=macOSX.i386
UCSC_OS=linux.x86_64

## Main directory for installing third-party software)
SOFT_DIR=${PWD}

## Directory to store the downloaded software (sources, before compilation)
SRC_DIR=${SOFT_DIR}/src

## Default installation dir for binaries
BIN_DIR=/usr/local/bin


################################################################
##
## Configuration of a PC cluster
##
## The parameter QUEUE_MANAGER is required for running tasks on a PC
## cluster. The same program (qsub) has different parameters depending on the
## queue manager. Thus, a qsub command for torque will not be understood by
## sge and reciprocally;
##
## supported: torque | sge
QUEUE_MANAGER=sge

## Name of the queue where the jobs have to be sent
QUEUE=

