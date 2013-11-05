#!/usr/bin/perl
#### this cgi script fills the HTML form for the program patser
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc

# $default{matrix_format} = "consensus";
# $default{matrix} = ""; ### [-m <Name of matrix file---default name is "matrix">]
# $default{pseudo_counts} = 1; ### [-b <Correction added to the elements of the alignment matrix (default: 1)>]
# $default{strands} = "both"; ### [-c <Score the complementary strand>]


################################################################
#### STILL TO BE TREATED
# [-R <Set the range for approximating a weight matrix with integers (default: 10000)>]
# [-e <Small difference for considering 2 scores equal (default: 0.000001)>]
# [-li <Determine lower-threshold score from adjusted information content>]
# [-lp <Determine lower-threshold score from a maximum ln(p-value)>]


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

&ReadMatrixFromFile();

### print the form ###
&RSA_header("patser", "form");
#&ListParameters;

### head
print "<CENTER>";
print "Scan a DNA sequence with a posiion-specific scoring matrix (PSSM)<BR>\n";
print "Program developed by <A HREF='mailto:hertz\@colorado.edu (Jerry Hertz)'>Jerry Hertz</A>. \n";
print "Web interface by <A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>Jacques van Helden</A>.<br>";
print "The stand-alone version of <i>patser</i> is available at <a target=_blank href=http://ural.wustl.edu/software.html>http://ural.wustl.edu/software.html</a></p>";
print "</CENTER><hr>";

print $query->start_multipart_form(-action=>"patser.cgi");

################################################################
#### sequence
print "<BR>\n";
&DisplaySequenceChoice();


################################################################
#### patser options
&DisplayPatserOptions();

################################################################
#### origin for calculating position
print "<BR>";
print "<A HREF='help.dna-pattern.html#origin'><B>Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['start',
				   'end'],
			 -default=>$default{origin});

################################################################
#### flanking residues for the matching sequences
print "&nbsp;"x10;
print "<A HREF='help.all-upstream-search.html#flanking'><B> flanking</B></A>\n";
print $query->textfield(-name=>'flanking',
			-default=>$default{flanking},
			-size=>2);


################################################################
### send results by email or display on the browser
print "<BR>\n";
&SelectOutput();

################################################################
### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"patser_form.cgi");

$demo_sequence = ">PHO5   pho5 upstream sequence, from -800 to -1
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
>PHO8   pho8 upstream sequence, from -800 to -1
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
>PHO11  pho11 upstream sequence, from -800 to -1
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
>PHO81  pho81 upstream sequence, from -800 to -1
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
>PHO84  pho84 upstream sequence, from -800 to -1
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
ACTCCACGAATACAATCCAA
";
$demo_matrix = "A    |    1    3    2    0    8    0    0    0    0    0    1    2
C    |    2    2    3    8    0    8    0    0    0    2    0    2
G    |    1    2    3    0    0    0    8    0    5    4    5    2
T    |    4    1    0    0    0    0    0    8    3    2    2    2";

print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->hidden(-name=>'alphabet',-default=>"a:t 0.325 c:g 0.175");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.patser.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_patser.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





