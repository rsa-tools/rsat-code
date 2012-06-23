################################################################
## Site-specific configuration for the RSAT makefiles
## 
## This file is automatically loade by the script
## ${RSAT}/makefiles/util.mk, which is itself loaded by all the RSAT
## makefile scripts.

## Operating system
## Supported: linux  | macosx
#OS=linux
OS=linux
ARCHITECTURE=x64

#UCSC_OS=macOSX.i386
UCSC_OS=linux.x86_64

## Default installation dir for binaries
BIN_DIR=/usr/local/bin

## Main directory for software installation
SOFT_DIR=${PWD}

## Directory to store the downloaded software (sources, before compilation)
SRC_DIR=${PWD}/src


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

## Email address of RSAT administrator
RSAT_ADMIN_EMAIL=your.address@your.domain

