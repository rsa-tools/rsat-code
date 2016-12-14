################################################################
## Create an image for the command-line use of RSAT.
## A Web server will be created in a separate Dockerfile.
################################################################

## Hot line for JvH: Olivier Sallou <olivier.sallou@irisa.fr>

## TO DO
## - remove the apache start + functional tests from the ubuntu installer (run != installer)
## - create a separate tester (tester != run != install)

FROM ubuntu:16.04
MAINTAINER Jacques van Helden <Jacquese.van-Helden@univ-amu.fr>

## bioshadock labels
LABEL    rsat.description="Regulatory Sequence Analysis Tools (RSAT)" \ 
    rsat.homepage="http://rsat.eu/" \ 
    rsat.version="2016-10-29" \ 
    rsat.copyright="© 2008-10 Genome Research Limited. All Rights Reserved." \ 
    rsat.license="GPL-2+" \ 
    rsat.binaries="/packages/rsat/bin, /packages/rsat/perl-scripts, /packages/rsat/python-scripts" \ 
    rsat.topics="['Cis-regulatory elements']" \ 
    rsat.tags="implemented-in::perl, role::cis-regulatory elements,motif discovery" \ 
    bioshadock.Vendor="RSAT"

## Avoid the red messages at build by setting non-interactive mode.
## See: https://github.com/phusion/baseimage-docker/issues/58
ENV DEBIAN_FRONTEND noninteractive

## We do not define any entrypoint since this is a software suite. 
## Users invoke the container to run a command-line
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y wget

## FOR LATER (once the installation scripts will be finished and tested locally)
## RUN wget http://download.rsat.eu/install_rsat_ubuntu16/

## Copy the install scripts at the root of the container
COPY . install_rsat_ubuntu16

## Start the installation.
## Each time a RUN succeeds, it is cahced -> does not need to be
## re-executed for next build.
## To force a re-build, try: 
##      docker build --no-cache .
WORKDIR /install_rsat_ubuntu16
RUN bash 00_config.bash
RUN bash 01_ubuntu16.4_packages.bash
RUN bash 02_python_packages.bash
RUN bash 04_download_rsat.bash
RUN bash 05_install_rsat.bash
RUN bash 06_perl_packages.bash
RUN bash 07_R-and-packages.bash
RUN bash 08_apache_config.bash
RUN bash 09_rsat_ws.bash
RUN bash 10_clean_unnecessary_files.bash

## Restore the interactive mode, which will be required to run the container
ENV DEBIAN_FRONTEND teletype

WORKDIR /
