source installer/00_config.bash

################################################################
################       RSAT installation        ################
################################################################

## Run the configuration script, to specify the environment variables.
cd ${RSAT}
perl perl-scripts/configure_rsat.pl --auto rsat_site=${RSAT_SERVER_NAME} \
    ucsc_tools=0 ensembl_tools=0 phylo_tools=0 compara_tools=0 \
    rsat_server_admin=RSAT_admin


## Parameters to change
##   rsat_site   rsat-vm-2016-03
##   rsat_server_admin    I don't specify it, because I don't want to receive notifications from all the VMs
## I activate the optional tools ucsc_tools and ensembl_tools, but not the other ones because they require many genomes (phylo tools) or big genomes (compara_tools, variation_tools).

## Load the (updated) RSAT environment variables
cd ${RSAT}; source RSAT_config.bashrc

## Check that the RSAT environment variable has been properly
## configured. Note: I also define it in the beginning of the script
## because I will beed it for the different installation chunks.
echo "RSAT path: ${RSAT}"
cd ${RSAT}

## Initialise RSAT folders
make -f makefiles/init_rsat.mk init

################################################################
## Next steps require to be done as rsat administrator user

## compile RSAT programs written in C
cd ${RSAT}
make -f makefiles/init_rsat.mk compile_all
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_app_compiled.txt

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!  BUG    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!! I HAVE A PROBLEM TO COMPILE KWALKS. SHOULD BE CHECKED !!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


################################################################
## Install some third-party programs required by some RSAT scripts.
cd ${RSAT}
make -f makefiles/install_software.mk list_ext_apps
make -f makefiles/install_software.mk install_ext_apps
df -m > ${RSAT}/install_logs/df_$(date +%Y-%m-%d_%H-%M-%S)_rsat_extapp_installed.txt


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
echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo
echo "Bash configuration for all users"
echo "    /etc/bash_completion.d/RSAT_config.bashrc"
sudo rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/

## THIS WAS CREATING THE PROBLEM WITH have command not found !
## THE COMPLETION FILES MUST BE EXECUTED WITH bash_completion
# echo '' >> /etc/bash.bashrc
# echo '## Custom bash completion'  >> /etc/bash.bashrc
# echo 'for file in /etc/bash_completion.d/* ; do'  >> /etc/bash.bashrc
# echo '    source "$file"'  >> /etc/bash.bashrc
# echo 'done'  >> /etc/bash.bashrc
# echo ''  >> /etc/bash.bashrc

echo "WARNING"
echo "In order to support bash completion, you must uncomment the following lines in /etc/bash.bashrc"
echo "# # enable bash completion in interactive shells"
echo "# if ! shopt -oq posix; then"
echo "#   if [ -f /usr/share/bash-completion/bash_completion ]; then"
echo "#     . /usr/share/bash-completion/bash_completion"
echo "#   elif [ -f /etc/bash_completion ]; then"
echo "#     . /etc/bash_completion"
echo "#   fi"
echo "# fi"

## ln -fs ${RSAT}/RSAT_config.bashrc /etc/bash_completion.d/

ln -fs ${RSAT} ~/rsat

#emacs -nw /etc/bash.bashrc

