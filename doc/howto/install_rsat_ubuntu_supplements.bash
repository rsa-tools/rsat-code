

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

################################################################
## Install some software tools for NGS analysis
################################################################

cd ${RSAT}
## TO BE DONE




################################################################
## Install gene-regulation package

export GR_PARENT_PATH=/packages
mkdir -p ${GR_PARENT_PATH}
cd ${GR_PARENT_PATH}
