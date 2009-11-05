#!/usr/bin/perl
#### this cgi script fills the HTML form for the program seq-proba
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

### Read the CGI query
$query = new CGI;

local @supported_input_formats = ("ft","gft","gff","gff3","dnapat");
local @supported_output_formats = ("ft","fasta","gff","gff3","dnapat");

my $input_formats = join (",",@supported_input_formats);

################################################################
### default values for filling the form
$default{output}="display";
$default{sequence_format} = "fasta";
$default{sequence} = "";

## BG model
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{bg_method}="bgfile";
#$default{bg_method}="markov";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = "0";
$default{organism} = "Saccharomyces cerevisiae";

## Return fields
@return_fields = qw(id proba_b len seq);
$default{id}="checked";
$default{proba_b}="checked";
$default{len}="checked";
$default{seq}="";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 



################################################################
### header
&RSA_header("seq-proba", "form");
print "<CENTER>";
print "Inter-conversions between various sequence formats.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"seq-proba.cgi");


################################################################
#### sequence
print "<hr>";
&DisplaySequenceChoice();

################################################################
## Background model
print "<hr>";

my %bg_params =("markov" => 1,
	       );
&GetBackgroundModel(\%bg_params);

################################################################
#### Return fields
print "<hr>";
print "<p><b><a href='help.seq-proba.html#return'>Return fields</a></b>&nbsp;<br>\n";
my $i = 0;
foreach my $field (@return_fields) {
  print $query->checkbox(-name=>$field,
			 -checked=>$default{$field},
			 -label=>'');
  print "&nbsp;<A HREF='help.seq-proba.html#",$field,"'><B>", $field, "</B></A>\n";
  print "&nbsp\n";
}
print "<p>\n";

################################################################
### send results by email or display on the browser
print "<hr>";
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

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"seq-proba_form.cgi");
$demo=">seq1
CACGTG
>seq2
CCGCGG
>seq3
TATAAA
";
print "<TD><B>";
print $query->hidden(-name=>'sequence',-default=>$demo);
print $query->hidden(-name=>'sequence_format',-default=>'fasta');
print $query->hidden(-name=>'output_format',-default=>"wconsensus");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.seq-proba.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);
