#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis_form.cgi,v 1.6 2001/12/24 01:45:19 jvanheld Exp $
#
# Time-stamp: <2001-12-24 02:45:11 jvanheld>
#
############################################################
#### this cgi script fills the HTML form for the program dyad-analysis
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
$default{organism} = "Saccharomyces cerevisiae";
#$default{title} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_file} = "";
$default{oligo_size} = 3;
$default{spacing_from} = 0;
$default{spacing_to} = 20;
$default{strand} = "both strands";
$default{noov} = '';
$default{purge} = 'checked';
$default{dyad_type} = "any dyad";
$default{exp_freq} = "dyad freq in non-coding sequences";
$default{occ_significance_threshold} = "0";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("dyad-analysis");

### head
print "<CENTER>";
print "Analysis of spaced dyads in a set of DNA sequences<P>\n";
print "</CENTER>";

print $query->start_multipart_form(-action=>"dyad-analysis.cgi");

#print "<FONT FACE='Helvetica'>";

### Title
#print "<B><A HREF='help.dyad-analysis.html#title'>Title</A></B>&nbsp;\n";
#print $query->textfield(-name=>'title',
#			-default=>$default{title},
#			-size=>50);
#
#print "<BR>\n";



&DisplaySequenceChoice;


#### purge sequences
print $query->checkbox(-name=>'purge',
  		       -checked=>$default{purge},
  		       -label=>'');
print "&nbsp;<A HREF='help.oligo-analysis.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";
print "<HR width=550 align=left>\n";

### oligo size
print "<B><A HREF='help.dyad-analysis.html#oligo_size'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_size',
			 -Values=>[1..4],
			 -default=>$default{oligo_size});

### spacing
print "<A HREF='help.dyad-analysis.html#spacing'><B>Spacing</B></A>&nbsp;\n";
print "&nbsp;", "from", "&nbsp;";
print $query->popup_menu(-name=>'spacing_from',
			 -Values=>[0..20],
			 -default=>$default{spacing_from});
print "&nbsp;", "to", "&nbsp;";
print $query->popup_menu(-name=>'spacing_to',
			 -Values=>[0..20],
			 -default=>$default{spacing_to});

print "<BR>\n";

### dyad type
print "<B><A HREF='help.dyad-analysis.html#dyad_type'>Dyad type</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'dyad_type',
			 -Values=>["inverted repeats",
				   "direct repeats",
				   "any repeat",
				   "any dyad"],
			 -default=>$default{dyad_type});

### strand ###
print "<BR>";
print "<B><A HREF='help.dyad-analysis.html#count_strands'>Count on</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

### prevent overlapping matches of the same pattern
print "&nbsp;" x 5;
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "<A HREF='help.dyad-analysis.html#noov'><B>\n";
print "prevent overlapping matches\n";
print "</B></A>\n";

print "<BR>\n";


print "<HR width=550 align=left>\n";


### expected frequency calculation
print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
		    $query->Tr($query->td("<A HREF='help.oligo-analysis.html#exp_freq'><B>Expected frequency calibration</B></A>&nbsp;<BR>")),
		    $query->Tr($query->td(["<INPUT TYPE='radio' NAME='exp_freq' VALUE='dyad freq in non-coding sequences' CHECKED>Dyad frequencies from all non-coding regions<BR>",
					   &OrganismPopUpString])),
		    $query->Tr($query->td([
					   "<INPUT TYPE='radio' NAME='exp_freq' VALUE='monad (word) freq in the input sequences.'>Monad (word) frequencies from the input sequences<BR>",
					   ])),
		    );

print "<HR width=550 align=left>\n";

### expected frequency calculation
#print "<A HREF='help.dyad-analysis.html#exp_freq'><B>Expected frequency</B></A>&nbsp;";
#print $query->popup_menu(-name=>'exp_freq',
#			 -values=>['dyad freq in non-coding sequences',
#				   'monad (word) freq in the input sequences'],
#			 -default=>$default{exp_freq});
#
#print "<BR>\n";


### significance threshold
print "<B><A HREF='help.dyad-analysis.html#threshold'>\n";
print "Threshold of significance</A> >= \n";
print $query->textfield(-name=>'occ_significance_threshold',
		  -default=>$default{occ_significance_threshold},
		  -size=>5);
print "<BR>\n";



### send results by e-mail or display on the browser
&SelectOutput;

print "<B>Warning !</B> dyad-analysis is time-consuming, especially if you select a wide spacing range. We recommend e-mail output.<BR>\n";

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"dyad-analysis_form.cgi");
$demo_sequence = ">GAL1	YBR020W upstream sequence, from -800 to -1
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
>GAL2	YLR081W upstream sequence, from -800 to -1
CATTAATTTTGCTTCCAAGACGACAGTAATATGTCTCCTACAATACCAGTTTCGCTGCAG
AAGGCACATCTATTACATTTACTGAGCATAACGGGCTGTACTAATCCAAGGAGGTTTACG
GACCAGGGGAACTTTCCAGATTCAGATCACAGCAATATAGGACTAGAAAATATCAGGTAG
CCGCACTCAACTTGTAACTGGCAACTACTTTGCATCAAACTCCAATTAAATGCGGTAGAA
TCTTTTCACAAAAGGTACTCAACGTCAATTCGGAAAGCTTCCTTCCGGAATGGCTTAAGT
AGGTTGCAATTTCTTTTTCTATTAGTAGCTAAAAATGGGTCACGTGATCTATATTCGAAA
GGGGCGGTTGCCTCAGGAAGGCACCGGCGGTCTTTCGTCCGTGCGGAGATATCTGCGCCG
TTCAGGGGTCCATGTGCCTTGGACGATATTAAGGCAGAAGGCAGTATCGGGGCGGATCAC
TCCGAACCGAGATTAGTTAAGCCCTTCCCATCTCAAGATGGGGAGCAAATGGCATTATAC
TCCTGCTAGAAAGTTAACTGTGCACATATTCTTAAATTATACAACATTCTGGAGAGCTAT
TGTTCAAAAAACAAACATTTCGCAGGCTAAAATGTGGAGATAGGATAAGTTTTGTAGACA
TATATAAACAATCAGTAATTGGATTGAAAATTTGGTGTTGTGAATTGCTCTTCATTATGC
ACCTTATTCAATTATCATCAAGAATAGTAATAGTTAAGTAAACACAAGATTAACATAATA
AAAAAAATAATTCTTTCATA
>GAL7	YBR018C upstream sequence, from -800 to -1
GAGAACTGGAAAGATTGTGTAACCTTGAAAAACGGTGAAACTTACGGGTCCAAGATTGTC
TACAGATTTTCCTGATTTGCCAGCTTACTATCCTTCTTGAAAATATGCACTCTATATCTT
TTAGTTCTTAATTGCAACACATAGATTTGCTGTATAACGAATTTTATGCTATTTTTTAAA
TTTGGAGTTCAGTGATAAAAGTGTCACAGCGAATTTCCTCACATGTAGGGACCGAATTGT
TTACAAGTTCTCTGTACCACCATGGAGACATCAAAAATTGAAAATCTATGGAAAGATATG
GACGGTAGCAACAAGAATATAGCACGAGCCGCGGAGTTCATTTCGTTACTTTTGATATCA
CTCACAACTATTGCGAAGCGCTTCAGTGAAAAAATCATAAGGAAAAGTTGTAAATATTAT
TGGTAGTATTCGTTTGGTAAAGTAGAGGGGGTAATTTTTCCCCTTTATTTTGTTCATACA
TTCTTAAATTGCTTTGCCTCTCCTTTTGGAAAGCTATACTTCGGAGCACTGTTGAGCGAA
GGCTCATTAGATATATTTTCTGTCATTTTCCTTAACCCAAAAATAAGGGAAAGGGTCCAA
AAAGCGCTCGGACAACTGTTGACCGTGATCCGAAGGACTGGCTATACAGTGTTCACAAAA
TAGCCAAGCTGAAAATAATGTGTAGCTATGTTCAGTTAGTTTGGCTAGCAAAGATATAAA
AGCAGGTCGGAAATATTTATGGGCATTATTATGCAGAGCATCAACATGATAAAAAAAAAC
AGTTGAATATTCCCTCAAAA
>GAL80	YML051W upstream sequence, from -800 to -1
TATCCTTTACGTTTTGACTTGGTGCTCGAAGATGCTTTCAGAGATGGTGCTTATCCTCAT
GTCTTTTGGGTTTGTCTTCAATACGGCAGCCGTTGTCTTGCAAACGGCCGCCTCTGCCAT
GGCAAAGAATGCTTTCCATGACGATCATCGTAGTGCCCAATTGGGTGCCTCTATGATGGG
TATGGCTTGGGCAAGTGTCTTTTTATGTATCGTGGAATTTATCCTGCTGGTCTTCTGGTC
TGTTAGGGCAAGGTTGGCCTCTACTTACTCCATCGACAATTCAAGATACAGAACCTCCTC
CAGATGGAATCCCTTCCATAGAGAGAAGGAGCAAGCAACTGACCCAATATTGACTGCCAC
TGGACCTGAAGACATGCAACAAAGTGCAAGCATAGTGGGGCCTTCTTCCAATGCTAATCC
GGTCACTGCCACTGCTGCTACGGAAAACCAACCTAAAGGTATTAACTTCTTCACTATAAG
AAAATCACACGAGCGCCCGGACGATGTCTCTGTTTAAATGGCGCAAGTTTTCCGCTTTGT
AATATATATTTATACCCCTTTCTTCTCTCCCCTGCAATATAATAGTTTAATTCTAATATT
AATAATATCCTATATTTTCTTCATTTACCGGCGCACTCTCGCCCGAACGACCTCAAAATG
TCTGCTACATTCATAATAACCAAAAGCTCATAACTTTTTTTTTTGAACCTGAATATATAT
ACATCACATATCACTGCTGGTCCTTGCCGACCAGCGTATACAATCTCGATAGTTGGTTTC
CCGTTCTTTCCACTCCCGTC
>MEL1	YBR184W upstream sequence, from -800 to -1
GCATACTCTACGTTATTTACAAAAATGTCGATATCCATCAAATTTTGTTTGGCGTACAGA
TTGTAGTTGTGGCTGCTACTGCAGGAAGTTTGACGTACAGATACGTCCATGATCCACTTG
CCAAAAGAAATCTCAAGGCTTCAATGGCGCTCGGCGCAATTTTGTTCTTATCTGGCTACA
TTTCGTGGCTACTTGATATACACTATTGTTCGTTCTGGGTGCACGTTAGAAGAAGTATTT
TGGCTTTACCACTTGGTGTACTGCTTGAACCACACGGATGGTGGCATATATTAACTGGTA
TGGGGATTTATTTCTACATTGTTTCTTTGGAACATTTAAGGGTCATTACGCTCAACGTCA
GCTGCAATTACCAGTTCATCTGGAGATGGAAAGTCTTCCCTGAACTGATATGGAAAGGGC
GCAAACCCTCAACAAGATATTCACTTGAACTATTTGGCCCATACGTAGAAGATCAATCAA
TTGAAGTTAAAAAGGAGAAGTAATAATTATAGCATAATATATATTCATAATGTATAGGCA
TATTTATTTTTTATTTTTTTTATTTCATGTTCTATTTAATGACGAATCACGAAGAAAATA
TATCTAAGAAAAGATCTTTTGAATCCTTGATTTGCGAATAGTTTAAATGACCCAGCTTAT
TGCTCTGGTGAAAAAAAACTTTGTGCGGTCTCAAAGCCGTCGGCGGCAAAATAACGTGAA
TTGATGAAAGTAAATAAACAAAACAAAATCTCTAATTGTTGTAACACAAATACTAAGAAA
TTTGTTAGCTAATTCGGGAC
>GCY1	YOR120W upstream sequence, from -800 to -1
GTCTTAGTATCTCATCTCATCTCAATTTCTATATTCCACTATAAAATTTTTCACTCTTTC
TGCGCGCGCCAATGTCCCCGCAACTACTCAATAGGTAACATGAGAATATTTCAGTTCGTA
AGAGAGAAGAGATGAAGTTATTTGGGCTCTTTGCTCGAGGTTACAGAAGGGCCGCATTAG
AGTGAATGAGCTGATGATATTTCGCCCAGTTCTACATTTTTTTTTTTTTGGAAGTATGAC
CTCTGTTAAATTTTTTTTTTTTTAAATTTCACTTTCTAAAGTCCCAGAAATCCGCTTGAA
TGTCTTACATATTGCAATGGATATGCTTGGGTGATCATACTTCCTGGCTTTAGATATTTG
AAACTTAACTCTTGTCAACAAACTTCCTATGGAGTGTATAAGAATTGTAAGTTATAACAC
CGGCGAACAATCGGGGCAGACTATTCCGGGGAAGAACAAGGAAGGGCGGTCTTTTCTCCC
TCATTGTCATAGCAAGGTCATTTCGCCTTCTCAGAAAGGGGTAGAATCAATCTAGCACGC
AGATTGCAAACACGGCTTAATAATATGCCTATCAGGCATTCACCCGTGTGACGAATCGCA
CACCGCTGCTCTCCTTAATTCCCTAGAGTAGAAACCGAGCTTTCAGGAAAAGACTACGGC
AGTAAAGAATTGCTTTACTGGGCGTATAAAACCGGGAGAATCAAGACATTCTAATGACTT
GATTCAGGATGAGAGCTTAATAGGTGCATCTTAGCAAGCTAAAATTTGGACAGCTCTCAT
TACTAAATTAAGATAGAAAA";
print "<TD><B>";
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>"fasta");
print $query->hidden(-name=>'spacing_from',-default=>"8");
print $query->hidden(-name=>'spacing_to',-default=>"12");
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
#print $query->hidden(-name=>'title',-default=>'upstream sequences from the yeast GAL genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.dyad-analysis.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.dyad-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_dyad-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);


