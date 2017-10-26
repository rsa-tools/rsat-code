#!/usr/bin/perl
############################################################
#
# $Id: RSAT_home.cgi,v 1.80 2013/10/09 07:04:55 jvanheld Exp $
#
# Time-stamp: <2003-10-22 11:53:22 jvanheld>
#
############################################################
#### this cgi script fills the HTML form for the program dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

print "Content-type:txt/html\n\n";
print "<html>";
print "<body>";

@orgs = &RSAT::OrganismManager::get_supported_organisms_web();
print scalar @orgs;

print "</body></html>";
exit(0);
