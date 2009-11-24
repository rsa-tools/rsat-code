#!/usr/bin/perl
#### this cgi script fills the HTML form for the program convert-matrix
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

### Read the CGI query
$query = new CGI;


### print the form ###
&RSA_header("Implant sites", "form");
&ListParameters() if ($ENV{rsat_echo} >=2);

### default values for filling the form
$default{sites_espp} = 0.001;
$default{sites} = "";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

# if (!($query->param("sites_filename") eq ""))
# {
#     $filename = $query->param("sites_filename");
#     $default{sites} = `cat $filename`;
# }

### head
print "<center>";
print "Randomly implant given sites into given sequences.\n";
print "<br>Reference: <a target='_blank' href=\"http://www.ncbi.nlm.nih.gov/pubmed/19689955\">Defrance & van Helden, Bioinformatics 2009</a>.</a>";
print "</center>";

print $query->start_multipart_form(-action=>"implant-sites.cgi");

#print "<font face='Helvetica'>";

#### Input sequences
&DisplaySequenceChoice;
print "<hr>";

#### Input sites
print "<B>Sites to implant</B> (in fasta format)<br />";
print $query->textarea(-name=>'sites',
                     -default=>$default{sites},
                     -rows=>10,
                     -columns=>40);
print "<hr>";

# #### Options
# print "<B>Options</B>";
# print "<br />";
# 
# #### Expected number of sites per positions
# print " <a href='help.implant-sites.html#sites_espp'>Excpected number of sites per position</a> ";
# print $query->textfield(-name=>'sites_espp',
#           -default=>$default{sites_espp},
#           -size=>3);
# print "<hr>";

### send results by email or display on the browser
print "<P>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.implant-sites.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_implant-sites.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

