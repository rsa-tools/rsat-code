#!/usr/bin/perl
#### this cgi script fills the HTML form for the program matrix-scan
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc
$default{bg_method} = "bginput";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = 0;
$default{window_size} = 200;
$default{organism} = "Saccharomyces cerevisiae";
$default{matrix_format} = "tab";
$default{pseudo_counts} = 1;

## Return fields
$default{return_sites} = "CHECKED";
$default{return_limits} = "CHECKED";
$default{return_rank} = "CHECKED";
$default{return_normw} = "";
$default{return_matrix} = "CHECKED";
$default{return_freq_matrix} = "";
$default{return_weight_matrix} = "";
$default{return_bg_model} = "";

## Threshold values
$default{lth_score} = "6";
$default{uth_score} = "none";
$default{lth_rank} = "none";
$default{uth_rank} = "none";
$default{lth_proba_M} = "none";
$default{uth_proba_M} = "none";
$default{lth_proba_B} = "none";
$default{uth_proba_B} = "none";
$default{lth_normw} = "none";
$default{uth_normw} = "none";

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
&RSA_header("matrix-scan");
#&ListParameters;

### head
print "<CENTER>";
print "Scan a DNA sequence with a profile matrix<BR>\n";
print "Program developed by <A HREF='mailto:jturatsi\@scmbb.ulb.ac.be (Jean Valery Turatsinze)'>Jean Val&eacute;ry Turatsinze</A><BR>";
print "Web interface by <A HREF='mailto:jvanheld\@scmbb.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A><P>";

print "</CENTER>";

print $query->start_multipart_form(-action=>"matrix-scan.cgi");

################################################################
#### sequence
print "<hr>";
&DisplaySequenceChoice();


################################################################
#### Matrix specification
print "<hr>";
print "<A HREF='help.patser.html#matrix'><B>\n";
print "Matrix</B></A>\n";
print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
print "<B>Format</B>&nbsp;";

#### matrix format
print $query->popup_menu(-name=>'matrix_format',
			 -Values=>[
				   'tab',
				   'meme',
				   'consensus',
				   'gibbs',
				   'transfac'
				  ],
			 -default=>$default{matrix_format});

#### text area to enter the matrix
print "<BR>\n";
print $query->textarea(-name=>'matrix',
		       -default=>$default{matrix},
		       -rows=>4,
		       -columns=>60);

################################################################
## Background model
print "<hr>";
print "<p><B>Background model</B><br>";

#### Markov order
print ("<b><a href=help.matrix-scan.html#markov_order>Markov Chain order</a></b> &nbsp;");
print $query->popup_menu(-name=>'markov_order',
			 -Values=>[0..8],
			 -default=>$default{markov_order});
print "<p>";


## Method to estimate the bg model

print ("<b><a href=help.matrix-scan.html#bg_method>Estimation method</a></b>");

#### Input sequences
print ("<br><INPUT TYPE='radio' NAME='bg_method' VALUE='bginput' $checked{'bginput'}>", 
       "<b>Estimate from input sequences</b>");

# ## Sliding window
# print ("<br><INPUT TYPE='radio' NAME='bg_method' VALUE='window' $checked{'window'}>",
#        "<b>Sliding window</b> &nbsp;");
# print $query->textfield(-name=>'window_size',
# 			-default=>$default{window_size},
# 			-size=>5);

#### Pre-defined background frequencies
print ("<br><INPUT TYPE='radio' NAME='bg_method' VALUE='bgfile' $checked{bgfile}>");
print ("<b>Genome subset</b> &nbsp; ", 
       $query->popup_menu(-name=>'background',
			  -Values=>["upstream","upstream-noorf","intergenic"],
			  -default=>$default{background}));
print  &OrganismPopUpString();
	
print "<hr>";

################################################################
#### pseudo-counts
print "<BR>\n";
print "<B><A HREF='help.patser.html#pseudo_counts'>Pseudo-counts</A>\n";
print $query->textfield(-name=>'pseudo_counts',
			-default=>$default{pseudo_counts},
			-size=>2);

################################################################
#### strands
print "<BR>\n";
print "<A HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['single',
				   'both'],
			 -default=>$default{strands});

################################################################
#### origin for calculating position
print "&nbsp;"x4,  "<A HREF='help.dna-pattern.html#origin'><B>Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['start',
				   'end'],
			 -default=>$default{origin});


################################################################
## Fields to return + thresholds
&ReturnTable();

################################################################
### send results by email or display on the browser
print "<hr>";
print "<BR>\n";
&SelectOutput();

################################################################
### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"matrix-scan_form.cgi");

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
$demo_matrix = "A | 1  3  2  0  8  0  0  0  0  0  1  2
C | 2  2  3  8  0  8  0  0  0  2  0  2
G | 1  2  3  0  0  0  8  0  5  4  5  2
T | 4  1  0  0  0  0  0  8  3  2  2  2";

print "<TD><B>";
print $query->hidden(-name=>'organism',-default=>'Saccharomyces_cerevisiae');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->hidden(-name=>'alphabet',-default=>"a:t 0.325 c:g 0.175");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.matrix-scan.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_matrix-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);


################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {
  print "<p><b>Return</b>\n";
  
  ### Sites
  @return_fields = qw(sites rank normw limits matrix freq_matrix weight_matrix bg_model);
  foreach my $field (@return_fields) {
    print $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }

#    print $query->checkbox(-name=>'return_limits',
#			 -checked=>$default{return_limits},
#			 -label=>' Limits ');
#  print $query->checkbox(-name=>'return_rank',
#			 -checked=>$default{return_rank},
#			 -label=>' Rank ');
#  print $query->checkbox(-name=>'return_normw',
#			 -checked=>$default{return_normw},
#			 -label=>' Normalized weights ');
#  print $query->checkbox(-name=>'return_model',
#			 -checked=>$default{return_model},
#			 -label=>' Background model ');


  print "<p><b>Thresholds</b>\n";
  print "<blockquote>\n";

  print $query->table({-border=>1,-cellpadding=>3,-cellspacing=>0},
		      $query->Tr({-align=>left,-valign=>TOP},
				 [
				  $query->th([" <A HREF='help.matrix-scan.html#return_fields'>Fields</A> ",
					      " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					      " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A> "]),

				  ### Threshold on score
				  $query->td(['score',
					      $query->textfield(-name=>'lth_score',
								-default=>$default{lth_score},
								-size=>5),
					      $query->textfield(-name=>'uth_score',
								-default=>$default{uth_score},
								-size=>5)
					     ]),

				  ### Threshold on proba_M
				  $query->td(['proba_M',
					      $query->textfield(-name=>'lth_proba_M',
								-default=>$default{lth_proba_M},
								-size=>5),
					      $query->textfield(-name=>'uth_proba_M',
								-default=>$default{uth_proba_M},
								-size=>5)
					     ]),

				  ### Threshold on proba_B
				  $query->td(['proba_B',
					      $query->textfield(-name=>'lth_proba_B',
								-default=>$default{lth_proba_B},
								-size=>5),
					      $query->textfield(-name=>'uth_proba_B',
								-default=>$default{uth_proba_B},
								-size=>5)
					     ]),

				  ### Threshold on normw
				  $query->td(['normw',
					      $query->textfield(-name=>'lth_normw',
								-default=>$default{lth_normw},
								-size=>5),
					      $query->textfield(-name=>'uth_normw',
								-default=>$default{uth_normw},
								-size=>5)
					     ]),

				  ### Threshold on rank
				  $query->td(['rank',
					      $query->textfield(-name=>'lth_rank',
								-default=>$default{lth_rank},
								-size=>5),
					      $query->textfield(-name=>'uth_rank',
								-default=>$default{uth_rank},
								-size=>5)
					     ]),

				 ]
				)
		     );
  print "</blockquote>\n";

}





