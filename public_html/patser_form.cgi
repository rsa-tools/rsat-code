#!/usr/bin/perl
#### this cgi script fills the HTML form for the program patser
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
$default{matrix_format} = "consensus";
$default{matrix} = ""; ### [-m <Name of matrix file---default name is "matrix">]
$default{matrix_is_weight} = ""; ### [-w <Matrix is a weight matrix>]
$default{matrix_is_vertical} = ""; ### [-v <Vertical matrix---rows correspond to positions>]
$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc
$default{pseudo_counts} = 1; ### [-b <Correction added to the elements of the alignment matrix (default: 1)>]
$default{alphabet_file} = ""; ### [-a <Name of ascii alphabet file---default name is "alphabet">]
$default{alphabet} = "a:t 0.3 c:g 0.2"; ### [-A <Ascii alphabet information>]
$default{case} = "insensitive"; ### [-CS <Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)>]
$default{strands} = "both"; ### [-c <Score the complementary strand>]

$default{lthreshold_method} = "manual";
### [-li <Determine lower-threshold score from adjusted information content>]
### [-lp <Determine lower-threshold score from a maximum ln(p-value)>]
$default{lthreshold} = "0"; ### [-ls <Lower-threshold score, inclusive (formerly the -l option)>]
$default{uthreshold} = "none"; ### [-u <Upper-threshold score, exclusive>]


$default{return} = "all matching positions";
$default{top_scores} = "1"; ### [-t <Print only the top scores>]
$default{sort} = "checked"; ### [-ds <Print top scores in order of decreasing score (default: print in order of position)>]



### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### if a matrix file is specified in the query,
### read matrix from this file
if (($matrix_file = $query->param("matrix_file")) &&
     (-e $matrix_file)) {
  open MATRIX, $matrix_file;
  while (<MATRIX>) {
    $default{matrix} .= $_;
  }
  close MATRIX;
}

### print the form ###
&RSA_header("patser");
#&ListParameters;

### head
print "<CENTER>";
print "Scan a DNA sequence with a profile matrix<BR>\n";
print "Program developed by <A HREF='mailto:hertz\@colorado.edu (Jerry Hertz)'>Jerry Hertz</A><BR>";
print "Web interface by <A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>Jacques van Helden</A><P>";

print "</CENTER>";

print $query->start_multipart_form(-action=>"patser.cgi");


################################################################
#### Matrix specification
print "<A HREF='help.patser.html#matrix'><B>\n";
print "Matrix</B></A>\n";
print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
print "<B>Format</B>&nbsp;";

#### matrix format
print $query->popup_menu(-name=>'matrix_format',
			 -Values=>['consensus',
				   'gibbs',
				   'transfac'
				   ],
			 -default=>$matrix_format);

#### mwight matrix
print "&nbsp;"x6;
print $query->checkbox(-name=>'matrix_is_weight',
		       -label=>"contains weigths",
		       -checked=>$default{matrix_is_weight});

#### vertical matrix
print "&nbsp;"x6;
print $query->checkbox(-name=>'matrix_is_vertical',
		       -label=>"vertical",
		       -checked=>$default{matrix_is_vertical});


### text area to enter the matrix
print "<BR>\n";
print $query->textarea(-name=>'matrix',
		       -default=>$default{matrix},
		       -rows=>4,
		       -columns=>60);

#### sequence
print "<BR>\n";
&DisplaySequenceChoice;


### strands
print "<BR>\n";
print "<A HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['single',
				   'both'],
			 -default=>$default{strands});

### return
print "<BR>\n";
print "<A HREF='help.patser.html#return'><B>Return</B></A>&nbsp;\n";
print $query->radio_group(-name=>'return',
			  -Values=>[ 'all matching positions'],
			  -default=>$default{return});

print "&nbsp;"x6;
print $query->radio_group(-name=>'return',
			  -Values=>['top values for each sequence'],
			  -labels=>{'top values for each sequence'=>''},
			  -default=>$default{return});
print $query->textfield(-name=>top_scores,
			-default=>$default{top_scores},
			-size=>3
			);
print "top value(s) for each sequence";

### pseudo-counts and thresholds
print "<BR>\n";
print CGI::table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       CGI::Tr({-align=>left,-valign=>MIDDLE},
		       [
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       "<B><A HREF='help.patser.html#pseudo_counts'>Pseudo-counts</A>\n",
			       $query->textfield(-name=>'pseudo_counts',
						       -default=>$default{pseudo_counts},
						       -size=>2),
			       "&nbsp;&nbsp;<b>Thresholds</b>",
			       "<A HREF='help.patser.html#lthreshold'><B> lower</B></A>",
			       $query->textfield(-name=>'lthreshold',
						 -default=>$default{lthreshold},
						 -size=>2),
			       "<A HREF='help.patser.html#uthreshold'><B> upper</B></A>",
			       $query->textfield(-name=>'uthreshold',
						 -default=>$default{uthreshold},
						 -size=>2)
			       ]),
			])
		 );

### alphabet
print "<BR>\n";
print "<B><A HREF='help.patser.html#alphabet'>\n";
print "Alphabet</A></b>\n";
print $query->textfield(-name=>'alphabet',
			-default=>$default{alphabet},
			-size=>50);

#### case sensitivity
print "<BR>\n";
print "<B><A HREF='help.patser.html#case'>\n";
print "Case</A></b>\n";
print "&nbsp;"x6;
print $query->popup_menu(-name=>'case',
			 -Values=>[ "sensitive", "insensitive","insensitive, but mark lowercases"],
			 -default=>$default{case});


### send results by e-mail or display on the browser
print "<BR>\n";
&SelectOutput;

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

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


#print "<TD><B><A HREF='demo.patser.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.patser.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





