#!/usr/bin/perl
############################################################
#
# $Id: load_classes.pl,v 1.4 2010/11/12 22:08:38 rsat Exp $
#
# Time-stamp: <2003-07-10 12:08:29 jvanheld>
#
# perl classes for the AMAZE project
# developed by Jacques van Helden (jvanheld@ucmb.ulb.ac.be)
# creation date : Dec 1999
#
############################################################

require "lib/util.pl";
#require "config.pl";
use strict;
no strict "refs";
no strict "vars";


#### load all class files
$class_dir = $ENV{RSAT}."/perl-scripts/parsers/classes/*.pl";
my @class_files = glob($class_dir);

if ((defined($main::verbose)) && ($main::verbose >= 3)){
    warn ("loading classes", "\n\t",  
	  join ("\n\t", 
		$main::RSA, 
		$class_dir,
		@class_files), 
	  "\n");
}
foreach my $file (@class_files) {
    require $file;
}




return 1;
