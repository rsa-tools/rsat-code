#!/usr/bin/perl
#### this cgi script fills the HTML form for the program oligo-analysis
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
$default{title} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_freq_file} = "";
$default{sequence_type} = "dna";
$default{oligo_length} = 6;
$default{background} = "upstream-noorf";
$default{markov_order} = 2;
$default{pseudo_weight} = "0.05";
$default{strand} = "both strands";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{purge} = 'checked';
#$default{purge} = '';

$default{freq_estimate} = "Background model";

$default{zscore} = '';
$default{lth_zscore} = 'none';
$default{uth_zscore} = 'none';

$default{rank} = 'checked';
$default{lth_rank} = "none";
$default{uth_rank} = "none";

$default{ratio} = '';
$default{lth_ratio} = "none";
$default{uth_ratio} = "none";

$default{occ} = 'checked';
$default{lth_occ} = "1";
$default{uth_occ} = "none";

$default{mseq} = '';
$default{uth_mseq} = "none";
$default{lth_mseq} = "none";

$default{proba} = 'checked';
$default{lth_occ_pro} = "none";
$default{uth_occ_pro} = "none";

$default{lth_occ_sig} = "0";
$default{uth_occ_sig} = "none";

$default{freq} = '';
$default{lth_observed_freq} = "none";
$default{uth_observed_freq} = "none";

$default{return}="fields";

### print the form ###
&RSA_header("oligo-analysis");
print "<CENTER>";
print "Analysis of oligonucleotide representation in a set of DNA sequences<P>\n";
print "</CENTER>";
print "<HR>";
print "<blockquote>";

#&ListParameters;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

print $query->start_multipart_form(-action=>"oligo-analysis.cgi");

### Title
#print "<B><A HREF='help.oligo-analysis.html#title'>Title</A></B>&nbsp;\n";
#print $query->textfield(-name=>'title',
#			-default=>$default{title},
#			-size=>50);
#print "<BR>\n";
#print "<HR width=550 align=left>\n";

print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td([&SequenceChoice()])
			]),
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td(["<B><A HREF='help.oligo-analysis.html#sequence_type'>Sequence type</A></B>".
			       $query->popup_menu(-name=>'sequence_type',
						  -Values=>["dna","protein","other"],
						  -default=>$default{sequence_type})
			       ])
			])

		 );

#### purge sequences
print $query->checkbox(-name=>'purge',
		       -checked=>$default{purge},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-analysis.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";

print "<HR width=550 align=left>\n";


### oligo size
print "<B><A HREF='help.oligo-analysis.html#oligo_length'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_length',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_length});

### prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-analysis.html#noov'><B>prevent overlapping matches</B></A>";
print "<BR>\n";

### strand ###
print "<B><A HREF='help.oligo-analysis.html#count_strands'>Count on</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

#### group patterns by pairs of reverse complement
print $query->checkbox(-name=>'grouprc',
		       -checked=>$default{grouprc},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-analysis.html#grouprc'><B>return reverse complements together in the output</B></A>";
print "<BR>";


print "<HR width=550 align=left>\n";



################################################################
#### estimation of expected frequencies
print "<A HREF='help.oligo-analysis.html#exp_freq'><B>Expected frequency calibration</B></A>&nbsp;<p>";


#### pre-defined background frequencies
print ( "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='background' CHECKED>", 
	"Predefined background frequencies");
print "<ul>";
print ( "<a href='help.oligo-analysis.html#background'>Background model</a> &nbsp;&nbsp;&nbsp;&nbsp;", 
	$query->popup_menu(-name=>'background',
			   -Values=>["upstream","upstream-noorf","intergenic"],
			   -default=>$default{background}));
	
print "<br>", &OrganismPopUpString;
print "</ul>";


print "<p>";

#### Markov chain model
print ("<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Markov Chain (higher order dependencies)'>", 
       "Markov Chain (higher order dependencies) ");

print "order &nbsp;";
print $query->textfield(-name=>'markov_order',
			-default=>$default{markov_order},
			-size=>5);
print "<p>";

#### Lexicon partitioning
print "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Lexicon partitioning'>Lexicon partitioning<p>";

#### Bernouilli model
print "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Residue frequencies from input sequence'>Residue frequencies from input sequence<p>";

#### equiprobable residues
print "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Equiprobable residues'>Equiprobable residues (<A HREF='help.oligo-analysis.html#equiprobable'>usually NOT recommended</a>)<p>";

#### custom expected frequency file
print "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='file_upload'><a href='help.oligo-analysis.html#upload_freq_file'>Upload your own expected frequency file</a><BR>";

print $query->filefield(-name=>'upload_freq_file',
			-default=>'starting value',
			-size=>30,
			-maxlength=>200);
print "<p>";

################################################################
#### pseudo-weights
print ("<b><a href=help.oligo-analysis.html#pseudo>Pseudo-weight</a></b> &nbsp;",
       $query->textfield(-name=>'pseudo_weight',
			 -default=>$default{pseudo_weight},
			 -size=>5));

print "<HR width=550 align=left>\n";




#print "<A HREF='help.oligo-analysis.html#exp_freq'><B>Expected frequency</B></A>&nbsp;";
#print $query->radio_group(-name=>'freq_estimate',
#			  -Values=>['Equiprobable residues',
#				    'Residue frequencies from input sequence',
#				    'Markov Chain (higher order dependencies)',
#				    'Lexicon partitioning',
#				    'Oligo frequencies from all intergenic regions'],
#			  -default=>$default{freq_estimate});
#print "<BR>";



#### table with all the statistics and thresholds
print "<h3>Return</h3>\n";

print ("<INPUT TYPE='radio' NAME='return' VALUE='fields' checked>", 
       "One row per pattern");


print "<BLOCKQUOTE>\n";
print $query->table({-border=>1,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			 [
			  $query->th([" <A HREF='help.oligo-analysis.html#return_fields'>Fields</A> ",
				   " <A HREF='help.oligo-analysis.html#thresholds'>Lower<BR>Threshold</A> ",
				   " <A HREF='help.oligo-analysis.html#thresholds'>Upper<BR>Threshold</A> "]),

			  ### occurrences
			  $query->td([$query->checkbox(-name=>'occ',
						    -checked=>$default{occ},
						    -label=>' Occurrences '),
				   $query->textfield(-name=>'lth_occ',
						     -default=>$default{lth_occ},
						     -size=>5),
				   $query->textfield(-name=>'uth_occ',
						     -default=>$default{uth_occ},
						     -size=>5)
				   ]),

			  ### binomial proba
			  $query->td([$query->checkbox(-name=>'proba',
						    -checked=>$default{proba},
						    -label=>' Binomial proba '),
				   $query->textfield(-name=>'lth_occ_pro',
						     -default=>$default{lth_occ_pro},
						     -size=>5),
				   $query->textfield(-name=>'uth_occ_pro',
						     -default=>$default{uth_occ_pro},
						     -size=>5)]),

			  ### significance index
			  $query->td([$query->checkbox(-name=>'proba',
						    -checked=>$default{proba},
						    -label=>' Significance '),
				   $query->textfield(-name=>'lth_occ_sig',
						     -default=>$default{lth_occ_sig},
						     -size=>5),
				   $query->textfield(-name=>'uth_occ_sig',
						     -default=>$default{uth_occ_sig},
						     -size=>5)
				   ]),

			  ### Z-scores
			  $query->td([$query->checkbox(-name=>'zscore',
						    -checked=>$default{zscore},
						    -label=>' Z-scores '),
				   $query->textfield(-name=>'lth_zscore',
						     -default=>$default{lth_zscore},
						     -size=>5),
				   $query->textfield(-name=>'uth_zscore',
						     -default=>$default{uth_zscore},
						     -size=>5)
				      ]),

			  ### frequencies
			  $query->td([$query->checkbox(-name=>'freq',
						    -checked=>$default{freq},
						    -label=>' Frequencies '),
				   $query->textfield(-name=>'lth_observed_freq',
						     -default=>$default{lth_observed_freq},
						     -size=>5),
				   $query->textfield(-name=>'uth_observed_freq',
						     -default=>$default{uth_observed_freq},
						     -size=>5)
				      ]),

			  ### matching sequences
			  $query->td([$query->checkbox(-name=>'mseq',
						    -checked=>$default{mseq},
						    -label=>' Matching sequences '),
				   $query->textfield(-name=>'lth_mseq',
						     -default=>$default{lth_mseq},
						     -size=>5),
				   $query->textfield(-name=>'uth_mseq',
						     -default=>$default{uth_mseq},
						     -size=>5)
				      ]),

			  ### ratio
			  $query->td([$query->checkbox(-name=>'ratio',
						    -checked=>$default{ratio},
						    -label=>' Obs/exp ratio '),
				   $query->textfield(-name=>'lth_ratio',
						     -default=>$default{lth_ratio},
						     -size=>5),
				   $query->textfield(-name=>'uth_ratio',
						     -default=>$default{uth_ratio},
						     -size=>5)
				      ]),

			  ### rank
			  $query->td([$query->checkbox(-name=>'rank',
						       -checked=>$default{rank},
						       -label=>' Rank '),
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
print "</BLOCKQUOTE>\n";
print ("<INPUT TYPE='radio' NAME='return' VALUE='table'>", 
       "One row per gene (occurrence counts only, email output recommended)", "<P>\n");

print ("<INPUT TYPE='radio' NAME='return' VALUE='distrib'>", 
       "Pattern count distribubtions, one row per pattern (occurrence counts only, email output recommended)", "<P>\n");


print "<HR width=550 align=left>\n";

### send results by email or display on the browser
&SelectOutput();


### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"oligo-analysis_form.cgi");
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
print "<TD><B>";
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'title',-default=>'upstream sequences from the yeast PHO genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.oligo-analysis.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.oligo-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_oligo-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


