#!/usr/bin/env bash

source $(dirname $0)/00_config.bash

echo
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!     BEWARE: INSTALLATION REQUIRES SUDO RIGHTS       !!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables

################################################################
## Configure the Apache Web server

export APACHE_CONFIG_FOLDER=/etc/apache2/
export APACHE_CONFIG_MAIN=${APACHE_CONFIG_FOLDER}/sites-available/000-default.conf
# on Centos: APACHE_CONFIG_FOLDER=/etc/httpd/ but the rest is different as well -> best solution is to edit manually

## Make a backup of original apache configuration
zip -ry apacheconf_backup_`date '+%Y-%m-%d_%H%M%S'`.zip ${APACHE_CONFIG_FOLDER}

## We need to replace some lines in the Apache configuration files
# nano ${APACHE_CONFIG_MAIN}
## Uncomment the following line:
# Include conf-available/serve-cgi-bin.conf
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#Include conf-available/serve-cgi-bin.conf|Include conf-available/serve-cgi-bin.conf|' \
     ${APACHE_CONFIG_MAIN}


## To avoid puzzling warning at apache start, set ServerName globally.
echo "" >> ${APACHE_CONFIG_FOLDER}/apache2.conf
echo "ServerName localhost" >> ${APACHE_CONFIG_FOLDER}/apache2.conf

## In the file ${APACHE_CONFIG_FOLDER}/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi
##
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddHandler cgi-script .cgi|AddHandler cgi-script .cgi|' ${APACHE_CONFIG_FOLDER}/mods-available/mime.conf
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddEncoding x-|AddEncoding x-|' ${APACHE_CONFIG_FOLDER}/mods-available/mime.conf

echo
echo "EDIT THE FILE ${APACHE_CONFIG_FOLDER}/mods-available/mime.conf"
echo "AND ADD THE FOLLOWINT LINES before the </IfModule> tag"
echo
cat installer/apache_mime_additions.conf
echo


## I also uncomment the following, for convenience
##        AddEncoding x-compress .Z
##        AddEncoding x-gzip .gz .tgz
##        AddEncoding x-bzip2 .bz2



## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Adapt the PHP parameters
## TO DO: adapt for Ubuntu 16.04.1, where php5 has been replaced by php7.
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# emacs -nw /etc/php5/apache2/php.ini
## Modify the following parameters
##      post_max_size = 100M
##      upload_max_filesize=200M
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


################################################################
## Configure RSAT web server

## Edit the file to replace [RSAT_PARENT_FOLDER] byt the parent directory
## of the rsat directory.
cd ${RSAT}
rsync -ruptvl RSAT_config.conf ${APACHE_CONFIG_FOLDER}/sites-enabled/rsat.conf

## OPTIONAL: since I am using this to install a virtual machine whose
## only function will be to host the RSAT server, I replace the normal
## default web folder by RSAT web folder.
##
## TO DO: CHECK IF THIS DOES NOT CREATE PROBLEMS
## perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|DocumentRoot /var/www/html|DocumentRoot ${RSAT}/public_html|' ${APACHE_CONFIG_FOLDER}/sites-available/000-default.conf

# apache2ctl restart
## The server will now immediately display RSAT home page when you
## type its IP address.

## You should now test the access to the RSAT Web server, whose URL is
## in the environment variable RSAT_WWW
echo $RSAT_WWW


## The following lines are required to activate cgi scripts.  Found at
## http://www.techrepublic.com/blog/diy-it-guy/diy-enable-cgi-on-your-apache-server/
chmod 755 /usr/lib/cgi-bin
chown root.root /usr/lib/cgi-bin
a2enmod cgi ## this is apparently required to enable cgi
a2dismod mpm_event
a2enmod mpm_prefork
a2enmod wsgi   ## Required for REST web services


## THIS DOES NOT BELONG TO THE INSTALLER BUT TO THE RUN / TEST
# service apache2 restart

## DONE: apache server configured and started
## You can check it by opening a Web connection to
## http://[IP]
################################################################
