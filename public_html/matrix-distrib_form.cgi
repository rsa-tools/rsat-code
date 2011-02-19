#!/usr/bin/perl
#### this cgi script fills the HTML form for the program convert-matrix
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
use RSAT::matrix;
use RSAT::MatrixReader;
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
use RSAT::MatrixReader;

### Read the CGI query
$query = new CGI;

local @supported_input_formats = sort(keys( %RSAT::MatrixReader::supported_input_format));
local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

################################################################
### default values for filling the form
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "tab";
$default{pseudo_weight}=1;
$default{decimals}=1;
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{bg_method}="bgfile";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = "1";
$default{organism} = "Saccharomyces cerevisiae";


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
&RSA_header("matrix-distrib", "form");
print "<CENTER>";
print "Computes the theoretical distribution of score probabilities of a given PSSM.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"matrix-distrib.cgi");

#print "<FONT FACE='Helvetica'>";

################################################################
#### Matrix specification
print "<hr>";
&GetMatrix();
print "<hr>";

my %bg_params =("markov" => 1,
		"markov_message" => 1
				);
&GetBackgroundModel(%bg_params);

print "<hr>";

print "<br/>";
print "<A HREF='help.convert-matrix.html#decimals'><B>score decimals</B></A>\n";
print $query->popup_menu(-name=>'decimals',
			 -Values=>['0',
				   '1','2'],
			 -default=>$default{decimals});

################################################################
### send results by email or display on the browser
print "<p>\n";
&SelectOutput("server");

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"matrix-distrib_form.cgi");
$demo_matrix = "
; MET4 matrix, from Gonze et al. (2005). Bioinformatics 21, 3490-500.
A |   7   9   0   0  16   0   1   0   0  11   6   9   6   1   8
C |   5   1   4  16   0  15   0   0   0   3   5   5   0   2   0
G |   4   4   1   0   0   0  15   0  16   0   3   0   0   2   0
T |   0   2  11   0   0   1   0  16   0   2   2   2  10  11   8";
print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'input_format',-default=>'tab');
print $query->hidden(-name=>'organism',-default=>'Saccharomyces_cerevisiae');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'markov_order',-default=>'0');


print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.matrix-distrib.html'>MANUAL</A></B></TD>\n";
#print "<TD><B>TUTORIAL</B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

