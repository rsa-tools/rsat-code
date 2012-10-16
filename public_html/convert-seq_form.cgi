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
$default{output_format}="fasta";
$default{addrc}="";
$default{line_width}=60;
$default{short_action}="no treatment";
$default{short_size}=30;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}



################################################################
### header
&RSA_header("convert-seq", "form");
print "<CENTER>";
print "Inter-conversions between various sequence formats.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"convert-seq.cgi");


################################################################
#### sequence
print "<hr>";
&DisplaySequenceChoice();

print "<hr>";
print "<h4>Sequence processing</h4>";

## Short sequences
print "<B><A HREF='help.convert-seq.html'>Short sequences</A></B>&nbsp;";
print $query->popup_menu(-name=>'short_action',
			 -Values=>['no treatment',
				   'mask',
				   'skip'],
			 -default=>$default{short_action});
print "&nbsp;"x2, "<B><A HREF='help.convert-seq.html'>min size</A></b>\n";
print $query->textfield(-name=>'short_size',
			-default=>$default{short_size},
			-size=>3);


## Add reverse complement
print "<br/>";
print $query->checkbox(-name=>'addrc',
		       -checked=>$default{addrc},
		       -label=>'');
print "<B><A HREF='help.convert-seq.html'>Add reverse complement</A></b>\n";

### Output format
print "<hr>";

print "<B><A HREF='help.convert-seq.html'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
			 -Values=>['fasta',
				   'wconsensus',
				   'raw',
				   "tab",
				   'multi'],
			 -default=>$default{output_format});
print "&nbsp;"x2, "<B><A HREF='help.convert-seq.html'>Line width</A></b>\n";
print $query->textfield(-name=>'line_width',
			-default=>$default{line_width},
			-size=>3);

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

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"convert-seq_form.cgi");
$demo=">YBR020w	GAL1 upstream sequence, from -800 to -1, size 800
CAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGA
ACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAA
CTGAGCTGTCATTTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAA
AGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTT
ACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTT
TGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCC
GCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCG
GTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCT
ACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAAC
CTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCC
TTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAA
ATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTT
CTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAAC
GTCAAGGAGAAAAAACTATA
";
print "<TD><B>";
print $query->hidden(-name=>'sequence',-default=>$demo);
print $query->hidden(-name=>'sequence_format',-default=>'fasta');
print $query->hidden(-name=>'output_format',-default=>"wconsensus");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.convert-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);
