#!/usr/bin/perl
############################################################
#
# $Id: purge-sequence_form.cgi,v 1.1 2002/01/07 02:00:55 jvanheld Exp $
#
# Time-stamp: <2002-01-07 03:00:16 jvanheld>
#
############################################################
#### this cgi script fills the HTML form for the program purge-sequence
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_file} = "";
$default{match_len} = 40;
$default{mismatches} = 3;
$default{expected} = "auto";
$default{both_strands} = "on";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("purge-sequence");

print "<CENTER>\n";
print "Discards duplicated segments from a sequence set<BR>\n";
print "Program developed by <A HREF='mailto:kurtz\@TechFak.Uni-Bielefeld.DE (Stefan Kurtz)'>Stefan Kurtz</A><BR>\n";
print "Web interface by <A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>Jacques van Helden</A>).\n";
print "</CENTER>\n";


print $query->start_multipart_form(-action=>"purge-sequence.cgi");

#print "<FONT FACE='Helvetica'>";

#### input sequence
&DisplaySequenceChoice;

### add reverse complement strand
print "&nbsp;" x 5;
print $query->checkbox(-name=>'both_strands',
		       -checked=>$default{both_strands},
		       -label=>'');
print "<A HREF='help.purge-sequence.html#both_strands'><B>\n";
print "purge reverse complement strand\n";
print "</B></A>\n";

print "<BR>\n";

### match length
print "<B><A HREF='help.purge-sequence.html#match_len'>\n";
print "Minimal match length</A>\n";
print $query->textfield(-name=>'match_len',
		  -default=>$default{match_len},
		  -size=>5);
print "<BR>\n";

### mismatches
print "<B><A HREF='help.purge-sequence.html#mismatches'>\n";
print "Maximal number of mismatches</A>\n";
print $query->textfield(-name=>'mismatches',
		  -default=>$default{mismatches},
		  -size=>5);
print "<BR>\n";


### send results by e-mail or display on the browser
&SelectOutput;

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"purge-sequence_form.cgi");
$demo_sequence = ">YBR020w	GAL1 upstream sequence, from -800 to -1, size 800
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
>YBR019c	GAL10 upstream sequence, from -800 to -1, size 800
CGGTTTAGCATCATAAGCGCTTATAAATTTCTTAATTATGCTCGGGCACTTTTCGGCCAA
TGGTCTTGGTAATTCCTTTGCGCTAGAATTGAACTCAGGTACAATCACTTCTTCTGAATG
AGATTTAGTCATTATAGTTTTTTCTCCTTGACGTTAAAGTATAGAGGTATATTAACAATT
TTTTGTTGATACTTTTATGACATTTGAATAAGAAGTAATACAAACTGAAAATGTTGAAAG
TATTAGTTAAAGTGGTTATGCAGCTTTTCCATTTATATATCTGTTAATAGATCAAAAATC
ATCGCTTCGCTGATTAATTACCCCAGAAATAAGGCTAAAAAACTAATCGCATTATCATCC
TATGGTTGTTAATTTGATTCGTTAATTTGAAGGTTTGTGGGGCCAGGTTACTGCCAATTT
TTCCTCTTCATAACCATAAAAGCTAGTATTGTAGAATCTTTATTGTTCGGAGCAGTGCGG
CGCGAGGCACATCTGCGTTTCAGGAACGCGACCGGTGAAGACGAGGACGCACGGAGGAGA
GTCTTCCGTCGGAGGGCTGTCGCCCGCTCGGCGGCTTCTAATCCGTACTTCAATATAGCA
ATGAGCAGTTAAGCGTATTACTGAAAGTTCCAAAGAGAAGGTTTTTTTAGGCTAAGATAA
TGGGGCTCTTTACATTTCCACAACATATAAGTAAGATTAGATATGGATATGTATATGGTG
GTAATGCCATGTAATATGATTATTAAACTTCTTTGCGTCCATCCAAAAAAAAAGTAAGAA
TTTTTGAAAATTCAATATAA
";
print "<TD><B>";
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>"fasta");
print $query->hidden(-name=>'match_len',-default=>"40");
print $query->hidden(-name=>'mismatches',-default=>"3");
print $query->hidden(-name=>'both_strands',-default=>"on");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.purge-sequence.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_purge-sequence.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



