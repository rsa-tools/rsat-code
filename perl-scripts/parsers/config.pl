#!/usr/bin/perl
############################################################
#
# $Id: config.pl,v 1.1 2003/10/28 22:40:38 jvanheld Exp $
#
# Time-stamp: <2003-10-02 18:06:22 jvanheld>
#
############################################################


################################################################
#
# configuration for the aMAZE parsers
#
################################################################


################################################################
# check that the program runs on a machine with the appropriate 
# configuration
$hostname = `hostname`;
$domainname = `domainname`;
chomp $hostname;

################################################################
# SQL export format
$supported_dbms{oracle} = 1;
$supported_dbms{postgresql} = 1;
$supported_dbms{mysql} = 1;

$default{dbms} = "mysql";
$default{host} = "localhost.localdomain";
$default{schema} = "amaze";
$default{user} = "amaze";
$default{password} = "amaze";


################################################################
# machine-specific configuration
#
# In principle, it is sufficient to change these variables 
# according to your local configuration, and the programs should run. 
if (($hostname =~ /^brol/i) ||
    ($hostname =~ /^gwilh/i))
{
    ### input directory. This directory must contain the mirror of the databases to be p
    $RSA= "/home/jvanheld/rsa-tools/";
#    $Databases = "/win/databases/downloads"; 
    $Databases = "/lin/genomics/downloads"; 
    ### main output directory: the result is exported to subdirectories of  this directory
#    $dir{skeletons} = "/win/amaze/amaze_data/pathway_skeletons/";
    $parsed_data = "/win/amaze/amaze_data/parsed_data"; 
    
    #### pathway skeletons
    $dir{skeletons} = "/win/amaze/amaze_team/didier_gonze/corrected_data";


} elsif ($hostname =~ /rsat/) {
    
    ### input directory. This directory must contain the mirror of the databases to be parsed
    $RSA= "/ext/jvanheld/rsa-tools/";
    $Databases = "/ext/jvanheld/downloads";

    ### main output directory: the result is exported to subdirectories of  this directory
    $parsed_data = "/rubens/dsk4/jvanheld/amaze/parsed_data";

    #### pathway skeletons
    $dir{skeletons} = "/rubens/dsk2/dgonze/annotation/corrected_data";

} elsif (($hostname =~ /merlin/) || ($hostname =~ /^n\d+$/)) { ### cluster at SCMBB
    
    ### input directory. This directory must contain the mirror of the databases to be parsed
    $RSA= "/home/jvanheld/rsa-tools/";
    $Databases = "/home/jvanheld/downloads";

    ### main output directory: the result is exported to subdirectories of  this directory
    $parsed_data = "/home/jvanheld/parsed_data";

    #### pathway skeletons
#    $dir{skeletons} = "/rubens/dsk2/dgonze/annotation/corrected_data";

} elsif (($domainname =~ /ucmb/) ||
	 ($domainname =~ /scmbb/)
	 ) {
    ### input directory. This directory must contain the mirror of the databases to be parsed
    $RSA= "/rubens/dks2/jvanheld/rsa-tools/";

    ### main output directory: the result is exported to subdirectories of  this directory
    $parsed_data = "/rubens/dsk4/jvanheld/amaze/parsed_data";

    #### pathway skeletons
    $dir{skeletons} = "/rubens/dsk2/dgonze/annotation/corrected_data";

} elsif (($domainname =~ /cifn/) || ($hostname eq "itzamna")) {
    ### input directory. This directory must contain the mirror of the databases to be parsed
    $Databases = "/home/jvanheld/amaze/downloads/";
    ### main output directory: the result is exported to subdirectories of  this directory
    $parsed_data = "/home/jvanheld/amaze/parsed_data/";

} else {
    warn "Warning: no configuration for host '$hostname'\n";
}

################################################################
# initialization
$start_time = `date +%Y-%m-%d.%H%M%S`;
$delivery_date = `date +%Y%m%d`;
chomp $delivery_date;

return 1;
