################################################################
## Optional elements for RSAT installation on specific servers.
##
## BEWARE: this script should not be run as such.  Each chunk of code
## should be read by a Human being in order to evaluate if it has to
## run for the local RSAT instance.

################################################################
##
## Web server with Hard-wired IP. Only for fixed servers, not for
## Virtual machines !
##
## Check the URL of the web services (RSAT_WS). By default, the server
## addresses the WS requests to itself (http://localhost/rsat) because
## web services are used for multi-tierd architecture of some Web
## tools (retrieve-ensembl-seq, NeAT).
cd $RSAT

## Get the current IP address
export IP="$(ifconfig | grep 'inet addr' |head -1 | cut -d ':' -f 2| cut -d ' ' -f 1)"
# export IP=192.168.56.101
echo ${IP}

## Set a fixed IP address for the web services
perl ${RSAT}/perl-scripts/configure_rsat.pl auto RSAT_WS=http://${IP}/rsat RSAT_WWWW=http://${IP}/rsat SUDO=''
#perl ${RSAT}/perl-scripts/configure_rsat.pl auto RSAT_WS=auto RSAT_WWWW=auto
# rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/


################################################################
## Last step before delivery: reset the passowrd of the RSAT
## administrator user (rsat), and define a user (vmuser).

## Create a user for the virtual machine
##
## This VM user is separate from the rsat user, which only serves to
## manage the RSAT software suite and related packages.
##
## For the sake of security, we force this user to change password at
## first login

## First delete this user (in case it was previously defined)
##  sudo userdel --remove vmuser

## Then create vmuser
sudo useradd --password `openssl passwd -1 -salt xyz tochng`\
    --home /home/vmuser \
    --create-home \
    --shell /bin/bash \
    --comment "VM user" \
    vmuser

## Force vmuser to change password at first login
sudo chage -d 0 vmuser

## Force rsat user to change password at first login
usermod --password `openssl passwd -1 -salt xyz tochng` rsat
sudo chage -d 0 rsat

## Add sudoer rights to vmuser and rsat users
sudo chmod 644 /etc/sudoers
sudo emacs -nw /etc/sudoers
## Find the following line
##     # User privilege specification
##     root    ALL=(ALL:ALL) ALL
## Below it, add the following line:
##     rsat    ALL=(ALL:ALL) ALL
##     vmuser  ALL=(ALL:ALL) ALL


################################################################
## We now change the owner of the RSAT package to the rsat user.
cd ${RSAT_PARENT_PATH}
chown -R rsat.rsat ${RSAT_HOME}


## New (2016-03-25) : for the IFB cloud I suppress the RSAT user, and
## install everything as root.

## Note (2016-10-17) : I could actually always do the whole
## installation as root, and if required create RSAT user only at the
## very end, and chown the rsat directory then.

# ## Create a specific user for RSAT. The user is named rsat
# sudo adduser rsat
# ## Full Name: Regulatory Sequence Analysis Tools admin

# ## Grant sudoer privileges to the rsat user (will be more convenient for
# ## installing Perl modules, software tools, etc)
# visudo
# ## then add the following line below "User privilege specification"
# # rsat    ALL=(ALL:ALL) ALL

################################################################
## Install some software tools for NGS analysis
################################################################

cd ${RSAT}

## TO BE DONE
##
## I should evaluate if I can use conda for this, it would greatly
## facilitate the portability.



################################################################
## Install gene-regulation package

export GR_PARENT_PATH=/packages
mkdir -p ${GR_PARENT_PATH}
cd ${GR_PARENT_PATH}

################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!    ONLY FOR THE VM ON THE IFB CLOUD    !!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## Adapt RSAT icon for IFB cloud (should be adapted depending on VM
## type).
cp ${RSAT}/public_html/images/ifb-logo-s.jpg   ${RSAT}/public_html/images/RSAT_icon.jpg

## replace the data directory by a link to
## a separate disk containing all RSAT data.
export RSAT_DATA_DIR=/root/mydisk/rsat_data
cd ${RSAT}/public_html
mv data/* ${RSAT_DATA_DIR}/
mv data/.htaccess ${RSAT_DATA_DIR}/
rmdir data
ln -s ${RSAT_DATA_DIR} data
cd $RSAT

################################################################
################ Install Sun Grid Engine (SGE) job scheduler
################################################################

## Check the number of core (processors)
grep ^processor /proc/cpuinfo

## Check RAM
grep MemTotal /proc/meminfo

## Beware, before installing the grid engine we need to modify
## manually tjhe file /etc/hosts
emacs -nw /etc/hosts
## Initial config (problematic) 
##    127.0.0.1       localhost       rsat-vm-2016-03
##    127.0.1.1      rsat-vm-2016-03
## Config to obtain: 
##    127.0.0.1       localhost       rsat-vm-2016-03
##    #127.0.1.1      rsat-vm-2016-03
apt-get install --quiet --assume-yes gridengine-client
apt-get install --quiet --assume-yes gridengine-exec
apt-get install --quiet --assume-yes gridengine-master
apt-get install --quiet --assume-yes gridengine-qmon 

qconf -aq default  ## aggregate a new queue called "default"
qconf -mq default  ## modify the queue "default"
qconf -as localhost ## aggregate the localhost tho the list of submitters
## -> set the following values
## hostlist              localhost

## Take all default parameters BUT For the SGE master parameter, type
## localhost (it must be the hostname)

## Test that jobs can be sent to the job scheduler



################################################################
## Ganglia: tool to monitor a cluster (or single machine)
## https://www.digitalocean.com/community/tutorials/introduction-to-ganglia-on-ubuntu-14-04
sudo apt-get install -y ganglia-monitor rrdtool gmetad ganglia-webfrontend
sudo cp /etc/ganglia-webfrontend/apache.conf /etc/apache2/sites-enabled/ganglia.conf
sudo apachectl restart
