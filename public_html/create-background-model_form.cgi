#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program convert-matrix
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
use RSAT::MarkovModel;
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "matrix_web_forms.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

local @supported_input_formats = sort(keys( %RSAT::MarkovModel::supported_input_formats));
#local @supported_output_formats = sort(keys( %RSAT::MarkovModel::supported_output_format));
local @supported_output_formats = ("tab","transitions","tables","patser","oligo-analysis", "meme", "MotifSampler");

################################################################
### default values for filling the form
$default{output}="server";
$default{output_format} = "transitions";
$default{markov_order} = "2";
$default{noov} = "CHECKED";
$default{sequence} = "";
$default{sequence_url} = "";

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
&RSA_header("create-background-model", "form");
print "<CENTER>";
print "Estimate the parameters of a Markov model from sequences.<p>\n";
#print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

# &ListDefaultParameters() if ($ENV{rsat_echo} >= 0);

print $query->start_multipart_form(-action=>"create-background-model.cgi", -id=>"form");

################################################################
#### sequence
print "<fieldset>
<legend><b><a class='iframe' href='help.formats.html'>Sequences </a></b></legend>";
#print &SequenceChoice();
&MultiSequenceChoice("Background sequences",1);
print "</fieldset><p/>";


################################################################
#### Background specification


print "<fieldset>
<legend><b><a class='iframe' href='help.formats.html'>Background specifications </a></b></legend>";

## markov order
print ("<b><a class='iframe' href=help.matrix-scan.html#markov_order>Markov order</a></b> &nbsp;");
print $query->popup_menu(-name=>'markov_order',
			       -Values=>[0..7],
			       -default=>$default{markov_order});

print "&nbsp;"x2, "<i>Markov order =  k-mer size of your subsequent analysis - 1. ie: markov order 5 for 6-mers</i>";
print "<br/>";
print "<p/>";

### Output bg format
print "<B><A class='iframe' HREF='help.convert-background-model.html#item__2dto_output_format'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
			 -Values=>[@supported_output_formats],
			 -default=>$default{output_format});
print "<BR>\n";


## overlap
# print ($query->checkbox(-name=>'noov',-checked=>$default{noov},-label=>''));
# print "<B><A class='iframe' HREF='help.convert-background-model.html#item__2dnoov'>prevent overlapping matches (noov)</A></b>\n";
# print "<br/>";

print "</fieldset><p/>";

print "<br/>";


################################################################
### send results by email or display on the browser
print "<p>\n";
&SelectOutput("display");

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
## Demo button
$demo_url1= $ENV{rsat_www}."/demo_files/Oct4_peaks_top1000.fa";

print '<script>
function setDemo(demo_url){
    $("#reset").trigger("click");
    sequence_url1.value = demo_url;
}
</script>';
print "<TD><b>";
print '<button type="button" onclick="setDemo('. "'$demo_url1'" .')">DEMO</button>';
print "</B></TD>\n";

print "<TD><B><A class='iframe' HREF='help.convert-background-model.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

