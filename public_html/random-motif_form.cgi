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

local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

### Read the CGI query
$query = new CGI;


### default values for filling the form
$default{motif_nb} = 1;
$default{width} = 12;
$default{conservation} = 0.80;
$default{multiplier} = 15;
$default{round} = "checked";
$default{output_format} = "tab";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("Random motif", "form");

### head
print "<center>";
print "Generate random motifs with a given level of conservation in each column.\n";
print "<br>Reference: <a target='_blank' href=\"http://www.ncbi.nlm.nih.gov/pubmed/19689955\">Defrance & van Helden, Bioinformatics 2009</a>.</a>";
print "</center>";


print $query->start_multipart_form(-action=>"random-motif.cgi");

#print "<font face='Helvetica'>";


print "<b>Parameters<b/><br />";
#### Motif width
print " <a href='help.random-motif.html#width'>Motif width</a> ";
print $query->textfield(-name=>'width',
			-default=>$default{width},
			-size=>3);

#### Conservation
print " <a href='help.random-motif.html#conservation'>Conservation</a> ";
print $query->textfield(-name=>'conservation',
			-default=>$default{conservation},
			-size=>3);

print "<hr />";
print "<b>Options<b/><br />";

#### Multiplier
print "<br>";
print " <a href='help.random-motif.html#multiplier'>Multiplier</A> ";
print $query->textfield(-name=>'multiplier',
			-default=>$default{multiplier},
			-size=>3);

### Round the values
print $query->checkbox(-name=>'round',
		       -checked=>$default{round},
		       -label=>'');
print "<a href='help.random-motif.html#round'>round frequencies</a>\n";


#### Number of motifs
print "<br>";
print " <a href='help.random-motif.html#motif_nb'>Number of motifs</A> ";
print $query->textfield(-name=>'motif_nb',
			-default=>$default{motif_nb},
			-size=>3);

### Output matrix format
print "<br>";
print "<b><a href='help.convert-matrix.html#output_format'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
			 -Values=>[@supported_output_formats],
			 -default=>$default{output_format});
print "<BR>\n";

print "<hr />";

### send results by email or display on the browser
print "<P>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.random-motif.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_random-motif.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

