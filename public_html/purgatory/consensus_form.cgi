#!/usr/bin/perl
############################################################
#
# $Id: consensus_form.cgi,v 1.11 2011/04/10 13:49:57 jvanheld Exp $
#
# Time-stamp: <2003-07-11 15:07:43 jvanheld>
#
############################################################
#### this cgi script fills the HTML form for the program consensus
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
$default{length} = 10;
$default{cycles} = "auto";
$default{matrices_to_save} = 1;
$default{alphabet} = "a:t 0.3 c:g 0.2";
$default{strands} = "include as a single sequence";
$default{symmetrical} = '';
$default{one_per_seq} = '';
$default{prior_freq} = '';

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("consensus", "form");

print "<CENTER>\n";
print "Matrix-based motif discovery using CONSENSUS<BR>\n";
print "Extract shared motifs from a set of unaligned sequences<BR>\n";
print "Program developed by Jerry Hertz. \n";
print "Web interface by <A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>Jacques van Helden</A>.<br>";
print "The stand-alone version of <i>consensus</i> is available at <a target=_blank href='http://stormo.wustl.edu/src'>http://stormo.wustl.edu/src/</a><P>";
print "</CENTER><hr>";


print $query->start_multipart_form(-action=>"consensus.cgi");

#print "<FONT FACE='Helvetica'>";

#### input sequence
&DisplaySequenceChoice;

### seed with first sequence
print "&nbsp";
print $query->checkbox(-name=>'seed',
		       -checked=>$default{seed},
		       -label=>'');
print "<A class='iframe' HREF='help.consensus.html#seed'><B>\n";
print " seed with first sequence\n";
print "</B></A>\n";
print "<BR>\n";

### treatment of reverse complement strand
print "<B><A class='iframe' HREF=help.consensus.html#strands>Complement strand</A></B>";
print $query->popup_menu(-name=>'strands',
			 -Values=>["ignore","include as separate sequence","include as a single sequence"],
			 -default=>$default{strands});
print "<BR>\n";

### matrix length
print "<B><A class='iframe' HREF='help.consensus.html#length'>\n";
print "Matrix length</A>\n";
print $query->textfield(-name=>'length',-id=>'length',
		  -default=>$default{length},
		  -size=>5);
#print "<BR>\n";

### symmetrical pattern
print "&nbsp";
print $query->checkbox(-name=>'symmetrical',
		       -checked=>$default{symmetrical},
		       -label=>'');
print "<A class='iframe' HREF='help.consensus.html#symmetrical'><B>\n";
print " Assume that the pattern is symmetrical\n";
print "</B></A>\n";
print "<BR>\n";

### expected number of matches
print "<B><A class='iframe' HREF='help.consensus.html#cycles'>\n";
print "Expected number of matches</A>\n";
print $query->textfield(-name=>'cycles',-id=>'cycles',
		  -default=>$default{cycles},
		  -size=>5);
#print "<BR>\n";



### one match per sequence
print "&nbsp";
print $query->checkbox(-name=>'one_per_seq',
		       -checked=>$default{one_per_seq},
		       -label=>'');
print "<A class='iframe' HREF='help.consensus.html#one_per_seq'><B>\n";
print " at least one match within each sequence\n";
print "</B></A>\n";
print "<BR>\n";


### matrices to save
print "<B><A class='iframe' HREF='help.consensus.html#matrices_to_save'>\n";
print "Number of matrices to save</A>\n";
print $query->textfield(-name=>'matrices_to_save',
		  -default=>$default{matrices_to_save},
		  -size=>5);
print "<BR>\n";

### alphabet
print "<B><A class='iframe' HREF='help.consensus.html#alphabet'>\n";
print "Alphabet</A>\n";
print $query->textfield(-name=>'alphabet',-id=>'alphabet',
			-default=>$default{alphabet},
			-size=>50);
print "<BR>\n";

### prior frequencies
print $query->checkbox(-name=>'prior_freq',-id=>'prior_freq',
		       -checked=>$default{prior_freq},
		       -label=>'');
print "<A class='iframe' HREF='help.consensus.html#prior_freq'><B>\n";
print " use the designated prior frequencies\n";
print "</B></A>\n";
print "<BR>\n";

### send results by email or display on the browser
&SelectOutput;

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"consensus_form.cgi");
$demo_sequence_raw = ">PHO5	pho5 upstream sequence, from -800 to -1
TTTTACACATCGGACTGATAAGTTACTACTGCACATTGGCATTAGCTAGGAGGGCATCCA
AGTAATAATTGCGAGAAACGTGACCCAACTTTGTTGTAGGTCCGCTCCTTCTAATAATCG
CTTGTATCTCTACATATGTTCTATTTACTGACCGAAAGTAGCTCGCTACAATAATAATGT
TGACCTGATGTCAGTCCCCACGCTAATAGCGGCGTGTCGCACGCTCTCTTTACAGGACGC
CGGAGACCGGCATTACAAGGATCCGAAAGTTGTATTCAACAAGAATGCGCAAATATGTCA
ACGTATTTGGAAGTCATCTTATGTGCGCTGCTTTAATGTTTTCTCATGTAAGCGGACGTC
GTCTATAAACTTCAAACGAAGGTAAAAGGTTCATAGCGCTTTTTCTTTGTCTGCACAAAG
AAATATATATTAAATTAGCACGTTTTCGCATAGAACGCAACTGCACAATGCCAAAAAAAG
TAAAAGTGATTAAAAGAGTTAATTGAATAGGCAATCTCTAAATGAATCGATACAACCTTG
GCACTCACACGTGGGACTAGCACAGACTAAATTTATGATTCTGGTCCCTGTTTTCGAAGA
GATCGCACATGCCAAATTATCAAATTGGTCACCTTACTTGGCAAGGCATATACCCATTTG
GGATAAGGGTAAACATCTTTGAATTGTCGAAATGAAACGTATATAAGCGCTGATGTTTTG
CTAAGTCGAGGTTAGTATGGCTTCATCTCTCATGAGAATAAGAACAACAACAAATAGAGC
AAGCAAATTCGAGATTACCA
>PHO8	pho8 upstream sequence, from -800 to -1
TCTTCACCAAATTTCTTTTTTTTTTCCTACTAGAAGAAGGCGTAGCAGATAAGAAGGAAA
AATTATATTAAGCGTGCGGGTAAAGGCAAGGAAGAATCAAGTAAGACCTCAAGAATGGCA
CTATAAGTGTGGTATTATAATCTGTGTAATCCTAATTTGAGCTCTACACAATACCATTCG
ACGGTTAACAGCTACTGCATCACCGTCCAGTCATGTCGTACAACGGAATAGGGCTCAAGT
CGGCAAAAGGGTCATCTACGTCGGGCCACGTGCAGCGATCACTTGCTAGCAACAATAGGC
GCAGACCACAGGGTAGTCAACAGCAGCGGCAACAACGACAAAATGCGATCAAAAAGGCCA
GCCATGACAAGGCAAGCAGGCCTCTTGCTGTGCAGAAACAGATAGAGACTCATATGGAGA
AACGTGAGATTGAAGTACAAGTTAGCGAGCTACGGGACCGACTGGAGGAGGAAGAAACGC
TCTCGGAAGAGCAGATTGACAAGAAATGTGAAGCGTTGAGGGCAAAACTGACGAACGAGT
GGCAAGAACAGCAGCGGATGTCCTCTTTGTACACCCCTCGTAAGGCGCGTCTAACGGAAG
AGCAGCATCGACATGAATAGCAGCATTGACGATAGCGATAAGCTTCGCGCGTAGAGGAAA
AGTAAAGGGATTTTAGTATATAAAGAAAGAAGTGTATCTAAACGTTTATATTTTTTCGTG
CTCCACATTTTGCCAGCAAGTGGCTACATAAACATTTACATATCAGCATACGGGACATTA
TTTGAACGCGCATTAGCAGC
>PHO11	pho11 upstream sequence, from -800 to -1
GCAGCCTCTACCATGTTGCAAGTGCGAACCATACTGTGGCCACATAGATTACAAAAAAAG
TCCAGGATATCTTGCAAACCTAGCTTGTTTTGTAAACGACATTGAAAAAAGCGTATTAAG
GTGAAACAATCAAGATTATCTATGCCGATGAAAAATGAAAGGTATGATTTCTGCCACAAA
TATATAGTAGTTATTTTATACATCAAGATGAGAAAATAAAGGGATTTTTTCGTTCTTTTA
TCATTTTCTCTTTCTCACTTCCGACTACTTCTTATATCTACTTTCATCGTTTCATTCATC
GTGGGTGTCTAATAAAGTTTTAATGACAGAGATAACCTTGATAAGCTTTTTCTTATACGC
TGTGTCACGTATTTATTAAATTACCACGTTTTCGCATAACATTCTGTAGTTCATGTGTAC
TAAAAAAAAAAAAAAAAAAGAAATAGGAAGGAAAGAGTAAAAAGTTAATAGAAAACAGAA
CACATCCCTAAACGAAGCCGCACAATCTTGGCGTTCACACGTGGGTTTAAAAAGGCAAAT
TACACAGAATTTCAGACCCTGTTTACCGGAGAGATTCCATATTCCGCACGTCACATTGCC
AAATTGGTCATCTCACCAGATATGTTATACCCGTTTTGGAATGAGCATAAACAGCGTCGA
ATTGCCAAGTAAAACGTATATAAGCTCTTACATTTCGATAGATTCAAGCTCAGTTTCGCC
TTGGTTGTAAAGTAGGAAGAAGAAGAAGAAGAAGAGGAACAACAACAGCAAAGAGAGCAA
GAACATCATCAGAAATACCA
>PHO81	pho81 upstream sequence, from -800 to -1
AAACGAGCATGAGGGTTACAAAGAACTTCCGTTTCAAAAATGAATATAATCGTACGTTTA
CCTTGTGGCAGCACTAGCTAACGCTACGTGGAATGAACGTACCGTGCCCTATTATTCTTG
CTTGTGCTATCTCAAGAATTGCATTTTGTAATAACAACTGCATGGGAAAAATTATATAGA
TTTTCTACTATTATGTCCGCCTAAGTCAGTTAACCATCTTTATCACAAAATATACAATTA
ACCAACTACTTAATCAATTCGGTTATATTGCTTAGTATATACGTCTTTGGCACGCGATTG
AAACGCGCTAATTGCATCAGCCTATCTTTCTATGCAAGAATGCAAGAAAAATTGATGTGA
TGTGCCTTATCACAATTCATTACCTCCTATTTCCTCTGCAGCAACAAGTTTCCTTGATTA
TAAAGGTCTTTAGCGTGAGAGGTACAGGTGTTATGGCACGTGCGAATAAGGGCAGAAATT
AATCAAATTTATCAACTATTTGGCGATGGCTCGAGACAGGTATAGAACCACTACTAGGTG
ATATTGAGGCTTTTGTACAATTTATAGCAAGTTTTTGAGAGTCCCTTCAAGTTTGTTACA
TAATCTTCTTTGTGCAACGTACAAGAGCAAAGTAGAAAAATTTGGTTTTTATTTTTTTAA
GCAACATCAGCTGCACTAGTTGAGCTTTTGACAAGACATACTGCTCAAAAAATCTTCATA
ACATTATTTTTCGGTTCCACAGTGATTGAGCTTTTTGAGAGAATAACCCTTTGGAGGCAA
CATAGATAGATAAACGTGCA
>PHO84	PHO84 upstream sequence, from -800 to -1
AAAAAAAAAGATTCAATAAAAAAAGAAATGAGATCAAAAAAAAAAAAAATTAAAAAAAAA
AAGAAACTAATTTATCAGCCGCTCGTTTATCAACCGTTATTACCAAATTATGAATAAAAA
AACCATATTATTATGAAAAGACACAACCGGAAGGGGAGATCACAGACCTTGACCAAGAAA
ACATGCCAAGAAATGACAGCAATCAGTATTACGCACGTTGGTGCTGTTATAGGCGCCCTA
TACGTGCAGCATTTGCTCGTAAGGGCCCTTTCAACTCATCTAGCGGCTATGAAGAAAATG
TTGCCCGGCTGAAAAACACCCGTTCCTCTCACTGCCGCACCGCCCGATGCCAATTTAATA
GTTCCACGTGGACGTGTTATTTCCAGCACGTGGGGCGGAAATTAGCGACGGCAATTGATT
ATGGTTCGCCGCAGTCCATCGAAATCAGTGAGATCGGTGCAGTTATGCACCAAATGTCGT
GTGAAAGGCTTTCCTTATCCCTCTTCTCCCGTTTTGCCTGCTTATTAGCTAGATTAAAAA
CGTGCGTATTACTCATTAATTAACCGACCTCATCTATGAGCTAATTATTATTCCTTTTTG
GCAGCATGATGCAACCACATTGCACACCGGTAATGCCAACTTAGATCCACTTACTATTGT
GGCTCGTATACGTATATATATAAGCTCATCCTCATCTCTTGTATAAAGTAAAGTTCTAAG
TTCACTTCTAAATTTTATCTTTCCTCATCTCGTAGATCACCAGGGCACACAACAAACAAA
ACTCCACGAATACAATCCAA";

$demo_sequence = join("\\n", split(/\n/, $demo_sequence_raw));

print '<script>
function setDemo(demo_sequence){
    $("#reset").trigger("click");
    sequence.value = demo_sequence;
    sequence_format.value = "fasta";
    length.value = "10";
    alphabet.value = "a:t 0.325 c:g 0.175";
    cycles.value = "10";
    $("#prior_freq").prop("checked", true);
}
</script>';
print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_sequence'" .')">DEMO</button>';
print "</B></TD>\n";


#print "<TD><B><A HREF='demo.consensus.html'>DEMO</A></B></TD>\n";
print "<TD><B><A class='iframe' HREF='help.consensus.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A class='iframe' HREF='tutorials/tut_consensus.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);


