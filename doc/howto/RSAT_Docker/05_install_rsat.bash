source 00_config.bash

################################################################
################       RSAT installation        ################
################################################################

## Run the configuration script, to specify the environment variables.
cd ${RSAT_HOME}
perl perl-scripts/configure_rsat.pl --auto rsat_site=rsat-vb-2016-10 \
    ucsc_tools=0 ensembl_tools=0 phylo_tools=0 compara_tools=0 \
    rsat_server_admin=RSAT_admin

## Parameters to change
##   rsat_site   rsat-vm-2016-03
##   rsat_server_admin    I don't specify it, because I don't want to receive notifications from all the VMs
## I activate the optional tools ucsc_tools and ensembl_tools, but not the other ones because they require many genomes (phylo tools) or big genomes (compara_tools, variation_tools).

## Load the (updated) RSAT environment variables
source RSAT_config.bashrc

## Check that the RSAT environment variable has been properly configured
echo ${RSAT}

## Initialise RSAT folders
make -f makefiles/init_rsat.mk init

################################################################
## Previous way to specify bashrc parameters, via
## /etc/bash_completion.d/. I change it (2014-09-23) because it does
## not allow to run remote commands via ssh (/etc/bash_completion.d is
## apparently only loaded in interactive mode).
## 
## Link the RSAT bash configuration file to a directory where files
## are loaded by each user at each login. Each user will then
## automatically load the RSAT configuration file when opening a bash
## session.
rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/
## ln -fs ${RSAT_HOME}/RSAT_config.bashrc /etc/bash_completion.d/

#emacs -nw /etc/bash.bashrc


