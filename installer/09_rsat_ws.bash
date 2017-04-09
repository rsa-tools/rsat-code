source installer/00_config.bash

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables


################################################################
## Configure the SOAP/WSDL Web services

## Check the URL of the web services (RSAT_WS). By default, the server
## addresses the WS requests to itself (http://localhost/rsat) because
## web services are used for multi-tierd architecture of some Web
## tools (retrieve-ensembl-seq, NeAT).
cd $RSAT



################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## THIS DOES NOT WORK ON VirtualBox VM 16.04.1
## CHECK IF I STILL NEED IT
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ## If the value is "auto", get the URL as follows
# export IP=`ifconfig  | awk '/inet /{print $2}' | cut -f2 -d':' | head -1`
# # export IP=192.168.56.101
# echo ${IP}
# export RSAT_WWW=http://${IP}/rsat/
# echo $RSAT_WWW
# export RSAT_WS=http://${IP}/rsat/
# echo $RSAT_WS

## Set a fixed IP address for the web services
# perl ${RSAT}/perl-scripts/configure_rsat.pl auto RSAT_WS=$RSAT_WWW RSAT_WWWW=$RSAT_WS SUDO=''

## Set localhost address for the web services
perl ${RSAT}/perl-scripts/configure_rsat.pl auto RSAT_WS=http://localhost/rsat RSAT_WWWW=http://localhost/rsat SUDO=''

## Set auto IP address for the web services
# perl ${RSAT}/perl-scripts/configure_rsat.pl auto RSAT_WS=auto RSAT_WWWW=auto SUDO=''

# rsync -ruptvl RSAT_config.bashrc /etc/bash_completion.d/ ## I DON'T THINK THIS IS REQUIRED HERE

################################################################
## Install the Web services
cd ${RSAT}
service apache2 restart ## Make sure the Apache server is running because it is required to generate the WS stub
make -f makefiles/init_rsat.mk ws_param  ## Check parameters to generate the stub for the Web services
make -f makefiles/init_rsat.mk ws_init  ## Initialize the stub for the Web services
make -f makefiles/init_rsat.mk ws_stub  ## Generate the stub for the Web services
make -f makefiles/init_rsat.mk ws_stub_test  ## Test the stub for the Web services
make -f makefiles/init_rsat.mk ws_nostub_test  ## Test Web services with no stub
service apache2 stop
