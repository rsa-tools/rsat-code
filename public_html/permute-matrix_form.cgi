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
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
use RSAT::matrix;
use RSAT::MatrixReader;

### Read the CGI query
$query = new CGI;

local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

################################################################
### default values for filling the form
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "transfac";
$default{output_format} = "transfac";
$default{perm} = 1;


&ReadMatrixFromFile();

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 


################################################################
### print the form ###


################################################################
### header
&RSA_header("permute-matrix", "form");
print "<CENTER>";
print "Randomize a set of input matrices by permuting their columns.<br/> The
    resulting motifs have the same nucleotide composition and information
    content as the original ones.<p/>
";
#print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"permute-matrix.cgi");

#print "<FONT FACE='Helvetica'>";

################################################################
#### Matrix specification
print "<hr>";

## Input matrix
&GetMatrix( 'nowhere'=>1,'no_pseudo'=>1);

print "<p/><hr>";


### Output matrix format
print "<br>";
print "<b><a href='help.convert-matrix.html#output_format'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
			 -Values=>[@supported_output_formats],
			 -default=>$default{output_format});
print "<BR>\n";


#### permutations
print "<BR>\n";
print "<B><A>Number of permutations</A></b>\n";
print $query->textfield(-name=>'perm',
			-default=>$default{perm},
			-size=>2);

print "<hr>";

################################################################
### send results by email or display on the browser
print "<p>\n";
&SelectOutput("display");

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"permute-matrix_form.cgi");
my $demo_matrix=`cat demo_files/compare-matrices_demo.tf`;
print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'input_format',-default=>'transfac');
print $query->hidden(-name=>'output_format',-default=>'transfac');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;



print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

