#!/usr/bin/perl
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
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

local @supported_input_formats = sort(keys( %RSAT::MarkovModel::supported_input_formats));
#local @supported_output_formats = sort(keys( %RSAT::MarkovModel::supported_output_format));
local @supported_output_formats = ("tab","transitions","tables","patser","oligo-analysis", "meme", "MotifSampler");

################################################################
### default values for filling the form
$default{output}="display";
$default{output_format} = "transitions";
$default{markov_order} = "2";
$default{bg_pseudo} = "0";
$default{bg_format}="oligo-analysis";
$default{bg_choose}="rsat";
$checked{$default{bg_choose}} = "CHECKED";
$default{bg_taxo}="organism";
$checked{$default{bg_taxo}} = "CHECKED";
$default{decimals}="3";
$default{organism} = "Saccharomyces cerevisiae";
$default{strands} = "single strand";
$default{noov} = "";

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
&RSA_header("convert-background-model", "form");
print "<CENTER>";
print "Interconversions between formats of background models supported by different programs.<P>\n";
#print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"convert-background-model.cgi");


################################################################
#### Background specification
print "<hr>";

my %bg_params = ("markov" => 1,
		 "title" => "RSAT pre-calculated background models",
		 "title_choose" => 1,
		 "noov" => 1,
		 "strands"=> 1,
		 "title2"=>"Custom background model",
		 #"taxon" => 1,
		 "sep_bg_pseudo" => 1
		);

&GetBackgroundModel(%bg_params);

print "<hr>";


### Output bg format
print "<BR>";
print "<B><A HREF='help.convert-background-model.html#item__2dto_output_format'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
			 -Values=>[@supported_output_formats],
			 -default=>$default{output_format});
print "<BR>\n";


print "<br/>";
print "<A HREF='help.convert-background-model.html#item__2ddecimals__23'><B>decimals</B></A>\n";
print $query->popup_menu(-name=>'decimals',
			 -Values=>[1..12],
			 -default=>$default{decimals});

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

print "<TD><B><A HREF='help.convert-background-model.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

