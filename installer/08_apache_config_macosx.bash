#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

cd ${RSAT};
# source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables
# have been reloaded by the source on 00_config.bash

################################################################
## Configure the Apache Web server

export APACHE_CONFIG_FOLDER=/etc/apache2/
# on Centos: APACHE_CONFIG_FOLDER=/etc/httpd/ but the rest is different as well -> best solution is to edit manually

## Make a backup of original apache configuration
zip -ry apacheconf_backup_`date '+%Y-%m-%d_%H%M%S'`.zip ${APACHE_CONFIG_FOLDER}

## Check the existence of the config file
export APACHE_CONF_MAIN=${APACHE_CONFIG_FOLDER}/httpd.conf
echo " APACHE_CONF_MAIN ${APACHE_CONF_MAIN}"
ls -ltr ${APACHE_CONF_MAIN}
tail $APACHE_CONF_MAIN

## Define the future location of the RSAT Apache configuration file
export CONF_DIR_USER=/private/etc/apache2/user/
export CONF_DIR_OTHER=/private/etc/apache2/other/
echo " CONF_DIR_USER ${CONF_DIR_USER}"
echo " CONF_DIR_OTHER ${CONF_DIR_OTHER}"
find ${APACHE_CONFIG_FOLDER} -name '*.conf' | xargs grep -i rsat | cut -d ':' -f 1 | sort -u

export APACHE_CONF_RSAT=${OTHER}rsat.conf
echo " APACHE_CONF_RSAT ${APACHE_CONF_RSAT}"


## Copy the RSAT-speciifc Apache config from RSAT folder to Apache config folder
export APACHE_CONF_RSAT_SRC=${RSAT}/RSAT_config.conf
echo " APACHE_CONFIG_RSAT_SRC ${APACHE_CONF_RSAT_SRC}"
ls -ltr ${APACHE_CONF_RSAT_SRC}
grep -v '^#'   ${APACHE_CONF_RSAT_SRC}

export SYNC_CMD="rsync -ruptvl ${APACHE_CONF_RSAT_SRC} ${APACHE_CONF_RSAT}"
echo ${SYNC_CMD}
${SYNC_CMD}

## Check that the config file is at the right position
ls -ltr ${APACHE_CONF_RSAT}
grep -v '^#' ${APACHE_CONF_RSAT}

## Add apache config file to the list of files to load
echo "" >>  ${APACHE_CONF_MAIN}
echo "## Regulatory Sequence Analysis Tools (rsat)" >>  ${APACHE_CONF_MAIN}
echo Include /private/etc/apache2/other/rsat.conf >>  ${APACHE_CONF_MAIN}
echo "" >>  ${APACHE_CONF_MAIN}
tail   ${APACHE_CONF_MAIN}

## Check the current status of Apache
apachectl status

## Check the Apache configuration (including our recent addition of rsat.conf)
apachectl configtest

## Stop Apache server
apachectl stop

## Check that the server is well stopped
## Try to open this URL in your browser:
##    http://localhost
## Should display an error message, since the server has been stopped.

## Start Apache server
apachectl start

## Check that the server is well restarted
## Try to open this URL in your browser:
##    http://localhost
## Should display the welcome page ("It works !")

## Try to open the RSAT Web server at this address
##    http://localhost/rsat


## Activate CGI
## Note: we keep a copy of the file in the state if was before modifying it
grep cgi ${APACHE_CONF_MAIN}
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddHandler cgi-script .cgi|AddHandler cgi-script .cgi|' ${APACHE_CONF_MAIN}
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#LoadModule cgi_module libexec/apache2/mod_cgi.so|LoadModule cgi_module libexec/apache2/mod_cgi.so|' ${APACHE_CONF_MAIN}
apachectl restart

echo("TO BE CONTINUED")
