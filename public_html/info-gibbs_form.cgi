#!/usr/bin/perl
############################################################
#
# $Id: info-gibbs_form.cgi
#
# 
#
############################################################
#### this cgi script fills the HTML form for the program info-gibbs
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
$default{length} = 12;
$default{expected} = 2.0;
$default{motifs} = 1;
$default{nrun} = 5;
$default{iter} = 1000;
$default{two_strands} = "checked";
$default{bg_order} = 3;
$default{background} = "upstream-noorf";
$default{bg_level} = "organism";
$default{organism} = "Saccharomyces cerevisiae";
$default{taxon} = "Fungi";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("info-gibbs", "form");
print '<style><!-- textarea {height: 100px; width: 510px;}--></style>';
print "<center>\n";
print "<p>An enhanced gibbs sampler, based on a stochastic optimization of the information content of position-specific scoring matrices.</p>\n";
print "</center>\n";
print "<p><b>Reference:</b> Defrance M, van Helden J. (2009) <i>Info-gibbs</i>: a motif discovery algorithm that directly optimizes information content during sampling. Bioinformatics 25(20):2715-22.\n";
print "[<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/19689955'>Pubmed 19689955</a>]\n";
print "[<a target='_blank' href='http://bioinformatics.oxfordjournals.org/content/25/20/2715.long'>Open acces</a>]\n";
print "<hr>\n";

print $query->start_multipart_form(-action=>"info-gibbs.cgi");

#print "<FONT FACE='Helvetica'>";

#### input sequence
&DisplaySequenceChoice;

### add reverse complement strand
print $query->checkbox(-name=>'two_strands',
		       -checked=>$default{two_strands},
		       -label=>'');
print "<a href=\"help.info-gibbs.html#strand\">search both strands</a>\n";


print "<hr width=\"550\" align=\"left\"></hr>\n";
print "<b>Options</b><br />\n";

### matrix length
print "<a href=\"help.info-gibbs.html#length\">Matrix length</a>\n";
print $query->textfield(-name=>'length',
		  -default=>$default{length},
		  -size=>5);
print "<br />\n";

### expected number of sites per sequence
print "<a href=\"help.info-gibbs.html#expected\">Expected number of sites per sequence</a>\n";
print $query->textfield(-name=>'expected',
		  -default=>$default{expected},
		  -size=>5);
#print "&nbsp;&nbsp;";
print "<br />\n";

### motifs
print "<a href=\"help.info-gibbs.html#motifs\">Number of motifs to extract</a>\n";
print $query->popup_menu(-name=>'motifs',
			 -Values=>[1,2,3,4,5],
			 -default=>$default{motifs});
#print $query->textfield(-name=>'motifs',
#			-default=>$default{motifs},
#			-size=>5);
print "<br />\n";


### iterations

print "<a href=\"help.info-gibbs.html#iter\">Maximum number of iterations</a>\n";
print $query->textfield(-name=>'iter',
		  -default=>$default{iter},
		  -size=>5);
#print "&nbsp;&nbsp;";
print "<br />\n";

### nrun
print "<a href=\"help.info-gibbs.html#runs\">Number of runs</a>\n";
print $query->popup_menu(-name=>'nrun',
			 -Values=>[1..10],
			 -default=>$default{nrun});
#print $query->textfield(-name=>'nrun',
#		  -default=>$default{nrun},
#		  -size=>5);
print "<br />\n";


## Background model
print "<br />\n";
print "<a href=\"help.info-gibbs.html#background\">Background model</a>\n";
print "<br />\n";
print '<input type="radio" checked="checked" value="input" name="freq_estimate"/><b>Estimated from input sequences</b><br />';
&PrintGenomeSubsetBgOptions();
print "<ul>";
print "&nbsp;&nbsp;<b>Markov background order</b> \n";
print $query->popup_menu(-name=>'bg_order',
			 -Values=>[0,1,2,3,4,5],
			 -default=>$default{bg_order});
print "</ul>";
### send results by email or display on the browser
print "<hr width=\"550\" align=\"left\"></hr>\n";
&SelectOutput;

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"info-gibbs_form.cgi");
$demo_sequence = ">dinG
TAAACCGCATTATGTTGGTGGTTATTGCGAGCCGCTTTCCAGAAACAGAAAAACCATTACCCCTGAAAACCGAAAAATGC
CACAATATTGGCTGTTTATACAGTATTTCAGGTTTTCTC
>dinQ
AGTTGGTGGCTCTGGCTGGAGTGAGAGAGTTCTTATCTAACAGCTCAATACTATTAGCCAGTGGCTAACATTGCAAGAAT
TAGCGATTTCTTCAGCTGGCGGAGGGAGCACTTCGGTATATATCGCAGTGACGATTTATATACTTTCACTGGGTCATCGT
CATATTAAGCCTTGTCCGAATCCGGATGTAAAAAAACCGACTTTGCGTCGGTTTTTTACTTTCCAGCCCTGAGTTGGTGG
CTCTGACTGGAGTGAGATAGCCATCATCTAACGCATCAATGCTAATACCATACTGAAACTATTGCAAGGACGTGCTGGTT
TTATAACCTGCATGTACTGTATGATTATCCAGTTAGCTCTGAGGCATTTTCACTCTGGCAATGCGCATAAACGCTTTCAA
AGTCCTGGTCAGAAGTACGGGTGGTGCCGTTAACTGATGCTCTGGCCGGAGTGAGAGAGTTCTTATCTAACAATGAGACA
TGCGCCGTGACAGGCAGTGG
>ftsK
CACGGAACAGGTGCAAAATCGGCGTATTTTGATTACACTCCTGTTAATCCATACAGCAACAGTACTGGGGTAACCTGGTA
CTGTTGTCCGTTTTAGCATCGGGCAGGAAAAGCCTGTAACCTGGAGAGCCTTTC
>lexA
CCCTTCCAGAATTCGATAAATCTCTGGTTTATTGTGCAGTTTATGGTTCCAAAATCGCCTTTTGCTGTATATACTCACAG
CATAACTGTATATACACCCAGGGGGCGGA
>yehH
TCCCCCTTCAAGGCGCTGCATCGACAGCGCCTTTTCTTTATAAATTCCTAAAGTTGTTTTCTTGCGATTTTGTCTCTCTC
TAACCCGCATAAATACTGGTAGCATCTGCATTCAACTGGATAAAATTACAGGGATGCAGA
>polB
TGACTGTATAAAACCACAGCCAATCAAACGAAACCAGGCTATACTCAAGCCTGGTTTTTTGATGGATTTTCAGC
>recA
TACTGTATGAGCATACAGTATAATTGCTTCAACAGAACATATTGACTATCCGGTATTACCCGGCATGACAGGAGTAAAA
>recA
TACTGTATGAGCATACAGTATAATTGCTTCAACAGAACATATTGACTATCCGGTATTACCCGGCATGACAGGAGTAAAA
>recN
TTTTACGCCAGCCTCTTTACTGTATATAAAACCAGTTTATACTGTACACAATAACAGTAATGGTTTTTCATACAGGAAAA
CGACT
>rpsU
GACTTGTTTTACCTCGCTTTATTACCGCGCAGTGTAGGACCAATGCGGGTTGATGTAAAACTTTGTTCGCCCCTGGAGAA
AGCCTCGTGTATACTCCTCACCCTTATAAAAGTCCCTTTCAAAAAAGGCCGCGGTGCTTTACAAAGCAGCAGCAATTGCA
GTAAAATTCCGCACCATTTTGAAATAAGCTGGCGTTGATGCCAGCGGCAAACCGAATTAATCAAAGGTGAGAGGCAC
>ruvA
GTTGTCATTCCATTGAAATAGATACAAATATACGTAATAGCATAAAAAATTTATTTTTTCAGGAGGCGCAATGTATATTT
TTATTACGCATTTCTTCACTGAATATGTAATATTAAAATATTTGCTTCCAATATAACCTGTAGAATAAATTATACTGTGC
CATTTTTCAGTTCATCGAGACACCTCGCAAGTTTTCTTCATCCTTCGCTGGATATCTATCCAGCATTTTTTTATCATACA
GCATTATCTTTGATTCATTACGCAGGAGCGTCAT
>ssb
TCACCTTTCCCGGATTAAACGCTTTTTTGCCCGGTGGCATGGTGCTACCGGCGATCACAAACGGTTAATTATGACACAAA
TTGACCTGAATGAATATACAGTATTGGAATGCATTACCCGGAGTGTTGTGTAACAATGTCTGGCCAGGTTTGTTTCCCGG
AACCGAGGTCACAACATAGTAAAAGCGCTATTGGTAATGGTACAATCGCGCGTTTACACTTATTCAGAACGATTTTTTTC
AGGAGACACGAAC
>sulA
AAAATTCCTTTTAAAATCATAACATAAAAGAATGATTCACATTAACGGATCCGTTAACTACGAAAATAGGCAACTTATTC
TTAAGGGGCAAGATTAATTTATGTTTTCCCGTCACCAACGACAAAATTTGCGAGGCTCTTTCCGAAAATAGGGTTGATCT
TTGTTGTCACTGGATGTACTGTACATCCATACAGTAACTCACAGGGGCTGGATTGATT
>tisA
TGGGGACCTTGTTGGTTTTGTGTTTAACAATATTTATACAAGCACAGCTTTACAGGGGAGACAATGGAAAATTTTTCAGC
AAGGGAAAATTGAGGGGTTGATCACGTTTTGTACTGAATTGCAGATAACAAAAAACCCCGCCGGAGCGAGGTTTCGTCAG
TCGCCTGCGGCTGGTAACCGCAAAGCACACTGTATTATGTCAACACTGAAAGTATACGTGTTCCGCGCAGAACGCGCAAT
TTCGGCACGAATTTTGACGTATTTAGTGCATAGTTGAGTATCGATCACAGTTTGCGTTTTGTCCAAATATTACTGTTTAT
TTATACAGTAAACTTCTATAATATCACTGTACGCAATGTGTTATGCGGGGGCCGCATCGTTACCCGGCGCACTAAGTCCT
GGCTGAAACGGGTGGTGCCGTCAGCGCCTTAACCCCGCGTGAGCACACTGTGTT
>umuD
AATCATTCGCCTCTTTAAATATATAAATTGTAATGAAACTCCTGTTTTACAACTATTAATAAATTTTACTTCATCTAATT
CATAGTTAGCCGGGCGGGATGCGTCAATGTCTTTATTTCTATTAATATGATAAATATCAAACAATGTTTAATGTCATTAT
GGCGAATGCTTCTATTCTATTTTTTAGCCGGGTGATATTTTTCATTTCTGCTGGATGAGCGTCGTCGCCAGAAGGCCACG
TGAGCACAAGATAAGAGAACGAAAAATCAGCAGCCTATGCAGCGACAAATATTGATAGCCTGAATCAGTATTGATCTGCT
GGCAAGAACAGACTACTGTATATAAAAACAGTATAACTTCAGGCAGATTATT
>uvrA
GTTCGTGTCTCCTGAAAAAAATCGTTCTGAATAAGTGTAAACGCGCGATTGTACCATTACCAATAGCGCTTTTACTATGT
TGTGACCTCGGTTCCGGGAAACAAACCTGGCCAGACATTGTTACACAACACTCCGGGTAATGCATTCCAATACTGTATAT
TCATTCAGGTCAATTTGTGTCATAATTAACCGTTTGTGATCGCCGGTAGCACCATGCCACCGGGCAAAAAAGCGTTTAAT
CCGGGAAAGGTGA
>uvrB
TGGCAGTCACTGAACAGGCATCTCTTGCCATAAAACTGTCATCACTCATCTTGACAAATGTTAAAAAAGCCGTTGCTTTG
GGGATAACCCGGTAAGGCCGGAGTTTTATCTCGCCACAGAGTAAATTTTGCTCATGATTGACAGCGGAGTTTACGCTGTA
TCAGAAATATTATGGTGATGAACTGTTTTTTTATCCAGTATAATTTGTTGGCATAATTAAGTACGACGAGTAAAATTACA
TACCTGCCCGCCCAACTCCTTCAGGTAGCGACTC
>uvrD
TCAGCAAATCTGTATATATACCCAGCTTTTTGGCGGAGGGCGTTGCGCTTCTCCGCCCAACCTATTTTTACGCGGCGGTG
CCA
>uvrY
TTATTTTCCTCGTCATGTTGCAATGAAAATTTGCGGTGAAAAATGTTAAGGCGGCGGAGTATACCATAAGCTTTGCTAAA
AATAGCAGTGGTTGTTTTTTGAGCGTGATATCGGCAGTGCTATAAAATACGTTAATATAATGACCAATAAATATTTTTAT
CATGAATGTTTTTTGCCCGTATTGTTCGGGTTAATTAATGTTACATATTCAGCGGGCTGATTTTCATTTTTGCTGAATAA
AGTCAATTTTTGTCACATTTCATCGTAGGGCTTACTGTGAAACGATCCGGTAAGCCGTTGGTGACGGGCGTGACCATAAC
TGTGGACAATCGAATTGACAAAAACGAGAGAAAAATCGAATACCCACCATTTTTAACGTTTCAAAGTTGCAATAAAAACC
GCTAATATACGAATGACTAACTATCAGTAGCGTTATCCCTATTTCTGGAGATATTCCT
>ydjM
TGACCAGGGGCAGTGATCTCGCTGCCCCTGGTTCTTTATCTGAAATTGCATTCAACTGACGGATTAATCGTCAATTTAAG
AGAAAGAGTTACACCGTCACCACTTCCGTGCACTGTATAAAAATCCTATACTGTACGTATCGACAGTTTA
";
print "<TD><B>";
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>"fasta");
print $query->hidden(-name=>'length',-default=>"20");
print $query->hidden(-name=>'expected',-default=>"1.0");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.info-gibbs.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_info-gibbs.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:defrance\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);


