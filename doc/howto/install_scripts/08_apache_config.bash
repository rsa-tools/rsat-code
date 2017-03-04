source 00_config.bash

cd ${RSAT}; source RSAT_config.bashrc ## Reload the (updated) RSAT environment variables

################################################################
## Configure the Apache Web server


## Make a backup of original apache configuration
zip -ry apacheconf_backup_`date '+%Y-%m-%d_%H%M%S'`.zip /etc/apache2/

## We need to replace some lines in the Apache configuration files
# nano /etc/apache2/sites-available/000-default.conf
## Uncomment the following line:
# Include conf-available/serve-cgi-bin.conf
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#Include conf-available/serve-cgi-bin.conf|Include conf-available/serve-cgi-bin.conf|' /etc/apache2/sites-available/000-default.conf

## To avoid puzzling warning at apache start, set ServerName globally.
echo "" >> /etc/apache2/apache2.conf
echo "ServerName localhost" >> /etc/apache2/apache2.conf

## In the file /etc/apache2/mods-available/mime.conf
## uncomment the line
##  AddHandler cgi-script .cgi
##
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddHandler cgi-script .cgi|AddHandler cgi-script .cgi|' /etc/apache2/mods-available/mime.conf
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|\#AddEncoding x-|AddEncoding x-|' /etc/apache2/mods-available/mime.conf

# echo "" >>  /etc/apache2/mods-available/mime.conf
# echo "        ## RSAT options" >>  /etc/apache2/mods-available/mime.conf
# echo "        AddType text/plain .fasta" >>  /etc/apache2/mods-available/mime.conf
# echo "        AddType text/plain .bed" >>  /etc/apache2/mods-available/mime.conf

## Optional : also associate a plain/text mime type to extensions for
## some classical bioinformatics files.
##   AddType text/plain .fasta
##   AddType text/plain .bed
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
rsync -ruptvl RSAT_config.conf /etc/apache2/sites-enabled/rsat.conf

## OPTIONAL: since I am using this to install a virtual machine whose
## only function will be to host the RSAT server, I replace the normal
## default web folder by RSAT web folder. 
##
perl -pi.`date '+%Y-%m-%d_%H%M%S'`.back -e 's|DocumentRoot /var/www/html|DocumentRoot /packages/rsat/public_html|' /etc/apache2/sites-available/000-default.conf

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

## THIS DOES NOT BELONG TO THE INSTALLER BUT TO THE RUN / TEST
#service apache2 restart

## DONE: apache server configured and started
## You can check it by opening a Web connection to 
## http://[IP]
################################################################
