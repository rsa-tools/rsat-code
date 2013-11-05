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
use RSAT::matrix;
use RSAT::MatrixReader;
require "patser.lib.pl";

local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

### Read the CGI query
$query = new CGI;


### default values for filling the form
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "tab";
$default{sites_nb} = 10;
#$default{output_format} = "tab";
$default{pseudo_counts} = 1;
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";

&ReadMatrixFromFile();

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("Random sites", "form");

### head
print "<center>";
print "Generate random sites given a motif (PSSM).\n";
print "<br>Reference: <a target='_blank' href=\"http://www.ncbi.nlm.nih.gov/pubmed/19689955\">Defrance & van Helden, Bioinformatics 2009</a>.</a>";
print "</center>";

print $query->start_multipart_form(-action=>"random-sites.cgi");

#print "<font face='Helvetica'>";

#### Input matrix
#print "<hr>";
&GetMatrix();
print "<hr>";


#### Options
print "<B>Options</B>";
print "<br />";

#### Number of sites
print " <a href='help.random-sites.html#sites_nb'>Number of sites</a> ";
print $query->textfield(-name=>'sites_nb',
			-default=>$default{sites_nb},
			-size=>3);
print "<hr>";

### send results by email or display on the browser
print "<P>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.random-sites.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_random-sites.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

