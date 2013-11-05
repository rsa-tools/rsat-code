#!/usr/bin/perl
############################################################
#
# $Id: purge-sequence_form.cgi,v 1.5 2009/11/03 10:06:12 jvanheld Exp $
#
# Time-stamp: <2003-10-01 12:17:22 jvanheld>
#
############################################################
#### this cgi script fills the HTML form for the program purge-sequence
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
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
$default{treatment} = "mask";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("purge-sequence", "form");

print "<CENTER>\n";
print "Discards duplicated segments from a sequence set<BR>\n";
print "Program developed by <a target=_blank href=http://www.techfak.uni-bielefeld.de/~kurtz/>Stefan Kurtz</a> (<A HREF='mailto:kurtz\@TechFak.Uni-Bielefeld.DE (Stefan Kurtz)'>kurtz\@TechFak.Uni-Bielefeld.DE</A>)<BR>\n";
print "Web interface by <a target=_blank href=http://www.bigre.ulb.ac.be/~jvanheld/>Jacques van Helden</a> (<A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>jvanheld\@bigre.ulb.ac.be</A>).\n";
print "</CENTER>\n";


print $query->start_multipart_form(-action=>"purge-sequence.cgi");

#print "<FONT FACE='Helvetica'>";

#### input sequence
&DisplaySequenceChoice();

### add reverse complement strand
print "&nbsp;" x 5;
print $query->checkbox(-name=>'both_strands',
		       -checked=>$default{both_strands},
		       -label=>'');
print "<A HREF='help.purge-sequence.html#both_strands'><B>\n";
print "purge reverse complement strand\n";
print "</B></A>\n";

print "<BR>\n";

### delete or mask repeats them
print "&nbsp;" x 5;
print "<B><A HREF='help.oligo-analysis.html#treatment'>Treatment for repeats</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'treatment',
			 -Values=>["delete","mask"],
			 -default=>$default{treatment});
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


### send results by email or display on the browser
&SelectOutput("server");

### action buttons
print "<UL><UL><TABLE class= 'formbutton'>\n";
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
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



