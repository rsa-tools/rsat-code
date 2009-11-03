#!/usr/bin/perl
#### this cgi script fills the HTML form for the program genome-scale-patser
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "$ENV{RSAT}/public_html/genome-scale.lib.pl";

### Read the CGI query
$query = new CGI;

$default{sequence_format} = "wc"; ### Important ! the format supported by patser


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
    if ($query->param($key)) {
	$default{$key} = $query->param($key);
    }
} 

### if a matrix file is specified in the query,
### read matrix from this file
if (($matrix_file = $query->param("matrix_file")) &&
     (-e $matrix_file)) {
  open PAT, $matrix_file;
  while (<PAT>) {
    $default{matrix} .= $_;
  }
  close PAT;
}

### print the form ###
&RSA_header("genome-scale patser", "form");

### head
print "<CENTER>";
print "Scan all upstream or downstream regions with a profile matrix<BR>\n";
print "Program developed by <A HREF='mailto:hertz\@colorado.edu (Jerry Hertz)'>Jerry Hertz</A>. \n";
print "Web interface by <A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>Jacques van Helden</A>.<BR>";
print "The stand-alone version of <i>patser</i> is available at <a target=_blank href=http://ural.wustl.edu/software.html>http://ural.wustl.edu/software.html</a><P>";
print "</CENTER><hr>";

#&ListParameters;


print $query->start_multipart_form(-action=>"genome-scale-patser.cgi");

################################################################
#
# retrieve-seq options
#

&DisplayRetrieveSeqOptions();

&DisplayPatserOptions();

### send results by email or display on the browser
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"genome-scale-patser_form.cgi");

$demo_matrix = "A    |    1    3    2    0    8    0    0    0    0    0    1    2
C    |    2    2    3    8    0    8    0    0    0    2    0    2
G    |    1    2    3    0    0    0    8    0    5    4    5    2
T    |    4    1    0    0    0    0    0    8    3    2    2    2";

print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'lthreshold_method',-default=>'weight');

print $query->hidden(-name=>'lthreshold',-default=>9);
print $query->hidden(-name=>'alphabet',-default=>"a:t 0.325 c:g 0.175");
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->hidden(-name=>'organism',-default=>'Saccharomyces_cerevisiae');
print $query->hidden(-name=>'set_name',-default=>'upstream sequences from the yeast PHO genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.patser.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_genome-scale-patser.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





