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

$default{lthreshold_method} = "adjusted information content (auto)";
$default{lthreshold} = "auto"; ### [-ls <Lower-threshold score, inclusive (formerly the -l option)>]
$default{uthreshold} = "none"; ### [-u <Upper-threshold score, exclusive>]


$default{return} = "all matches";
$default{top_scores} = "3"; ### [-t <Print only the top scores>]
$default{positions} = "checked"; ### convert the result into a score table
$default{table} = ""; ### convert the result into a score table
$default{sort} = "checked"; ### [-ds <Print top scores in order of decreasing score (default: print in order of position)>]

$default{unrecognized} = "discontinuities (with warning)"; ### [-d1 <Treat unrecognized characters as discontinuities, but print warning (the default)>]

$default{vertically_print} = "checked"; ### [-p <Vertically print the weight matrix>]
$default{min_calc_P} = 0; # [-M <Set the minimum score for calculating the p-value of scores (default: 0)>]

#### additional options
$default{origin} = "end";
$default{flanking} = "4";


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

#### weight matrix
print "&nbsp;"x6;
print $query->checkbox(-name=>'matrix_is_weight',
		       -label=>" contains weigths",
		       -checked=>$default{matrix_is_weight});

#### vertical matrix
print "&nbsp;"x6;
print $query->checkbox(-name=>'matrix_is_vertical',
		       -label=>" vertical",
		       -checked=>$default{matrix_is_vertical});


### text area to enter the matrix
print "<BR>\n";
print $query->textarea(-name=>'matrix',
		       -default=>$default{matrix},
		       -rows=>4,
		       -columns=>60);

################################################################
#### sequence
print "<BR>\n";
&DisplaySequenceChoice;


################################################################
### strands
print "<BR>\n";
print "<A HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['single',
				   'both'],
			 -default=>$default{strands});

################################################################
### return value
print "<BR>\n";
print "<table border=0>\n";
print "<tr><td>\n";
print "<A HREF='help.patser.html#return'><B>Return</B></A>&nbsp;\n";
print "</td><td>\n";
print $query->radio_group(-name=>'return',
			  -Values=>[ 'all matches'],
			  -default=>$default{return});

print "</td><td>\n";
print $query->radio_group(-name=>'return',
			  -Values=>['top values for each sequence'],
			  -labels=>{'top values for each sequence'=>''},
			  -default=>$default{return});
print $query->textfield(-name=>top_scores,
			-default=>$default{top_scores},
			-size=>3
			);
print "&nbsp;top value(s) for each sequence";

print "</td></tr><tr><td>\n";
print "&nbsp;";
print "</td><td>\n";
print $query->checkbox(-name=>'positions',
		       -label=>' matching positions',
		       -checked=>$default{positions});

print "</td><td>\n";
print $query->checkbox(-name=>'table',
		       -label=>' score table',
		       -checked=>$default{table});

print "</td></tr>\n";
print "</table>\n";


################################################################
### pseudo-counts
print "<BR>\n";
print "<B><A HREF='help.patser.html#pseudo_counts'>Pseudo-counts</A>\n";
print $query->textfield(-name=>'pseudo_counts',
			-default=>$default{pseudo_counts},
			-size=>2);

#### TEMPORARILY DISACTIVATED, BECAUSE INTERFERES WITH ADJUSTED INFO THRESHOLD
################################################################
### [-M <Set the minimum score for calculating the p-value of scores (default: 0)>]
# print "<BR>\n";
# print "<B><A HREF='help.patser.html#min_calc_P'>Minimum score for calculating the p-value</A>\n";
# print $query->textfield(-name=>'min_calc_P',
# 			-default=>$default{min_calc_P},
# 			-size=>5);

################################################################
#### Thresholds

#### lower threshold
print "<br>\n";
print "<A HREF='help.patser.html#lthreshold'><B>Lower threshold estimation</B></A>";
print $query->popup_menu(-name=>'lthreshold_method',
			 -Values=>['weight', 'maximum ln(p-value)', , 'adjusted information content (auto)', 'none'],
			 -default=>$default{lthreshold_method});

print $query->textfield(-name=>'lthreshold',
			-default=>$default{lthreshold},
			-size=>6);


#### upper threshold
print "<br>\n";
print "<A HREF='help.patser.html#uthreshold'><B>Upper threshold</B></A>", "&nbsp"x6;
print $query->textfield(-name=>'uthreshold',
			-default=>$default{uthreshold},
			-size=>6);

################################################################
### alphabet
print "<BR>\n";
print "<B><A HREF='help.patser.html#alphabet'>\n";
print "Alphabet</A></b>\n";
print $query->textfield(-name=>'alphabet',
			-default=>$default{alphabet},
			-size=>50);

################################################################
#### case sensitivity
print "<BR>\n";
print "<B><A HREF='help.patser.html#case'>\n";
print "Case</A></b>\n";
print "&nbsp;"x6;
print $query->popup_menu(-name=>'case',
			 -Values=>[ "sensitive", "insensitive","insensitive, but mark lowercases"],
			 -default=>$default{case});

################################################################
#### unrecognized characters
print "<BR>\n";
print "<B><A HREF='help.patser.html#unrecognized'>\n";
print "Treat unrecognized characters as</A></b>\n";
print "&nbsp;"x2;
print $query->popup_menu(-name=>'unrecognized',
			 -Values=>[ "errors", "discontinuities (with warning)","discontinuities (no warning)"],
			 -default=>$default{unrecognized});

################################################################
#### vertically print the matrix
print "<BR>\n";
print "<a href=help.patser.html#vertically_print>";
print $query->checkbox(-name=>'vertically_print',
		       -label=>' print the weight matrix',
		       -checked=>$default{vertically_print});
print "</a>";

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
### send results by e-mail or display on the browser
print "<BR>\n";
&SelectOutput;

################################################################
### action buttons
print "<UL><UL><TABLE>\n";
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
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





