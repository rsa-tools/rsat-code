#!/usr/bin/perl
#### this cgi script fills the HTML form for the program oligo-analysis
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$output_context = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{organism} = "Saccharomyces cerevisiae";
$default{title} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{exp_freq_file} = "";
$default{sequence_type} = "dna";
$default{oligo_size} = 6;
$default{markov_order} = 2;
$default{strand} = "both strands";
$default{noov} = '';
$default{grouprc} = 'checked';
$default{freq_estimate} = "Oligo frequencies from all non-coding regions";
$default{occ} = 'checked';
$default{proba} = 'checked';
$default{zscore} = '';
$default{freq} = '';
$default{mseq} = '';
$default{ratio} = '';
$default{occurrence_threshold} = "1";
$default{ms_threshold} = "none";
$default{proba_occ_threshold} = "none";
$default{occ_significance_threshold} = "0";

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

print "<HR width=550 align=left>\n";

### sequence type
#print "<B><A HREF='help.oligo-analysis.html#sequence_type'>Sequence type</A>&nbsp;</B>\n";
#print $query->popup_menu(-name=>'sequence_type',
#			 -Values=>["dna","protein","other"],
#			 -default=>$default{sequence_type});
#print "<BR>\n";

### oligo size
print "<B><A HREF='help.oligo-analysis.html#oligo_size'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_size',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_size});

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

print $query->checkbox(-name=>'grouprc',
		       -checked=>$default{grouprc},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-analysis.html#grouprc'><B>return reverse complements together in the output</B></A>";
print "<BR>";



print "<HR width=550 align=left>\n";

### expected frequency calculation
print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
		    $query->Tr($query->td("<A HREF='help.oligo-analysis.html#exp_freq'><B>Expected frequency calibration</B></A>&nbsp;<BR>")),
		    $query->Tr($query->td(["<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Oligo frequencies from all non-coding regions' CHECKED>Oligo frequencies from all non-coding regions<BR>",
					   &OrganismPopUpString])),
		    $query->Tr($query->td([
					   "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Markov Chain (higher order dependencies)'>Markov Chain (higher order dependencies)<BR>",
					   "Markov order &nbsp;".$query->textfield(-name=>'markov_order',
										   -default=>$default{markov_order},
										   -size=>5),
					   ])),
		    $query->Tr($query->td("<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Lexicon partitioning'>Lexicon partitioning<BR>")),
		    $query->Tr($query->td("<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Residue frequencies from input sequence'>Residue frequencies from input sequence<BR>")),
		    $query->Tr($query->td("<INPUT TYPE='radio' NAME='freq_estimate' VALUE='Equiprobable residues'>Equiprobable residues<BR>"))
		    );

print "<HR width=550 align=left>\n";






#print "<A HREF='help.oligo-analysis.html#exp_freq'><B>Expected frequency</B></A>&nbsp;";
#print $query->radio_group(-name=>'freq_estimate',
#			  -Values=>['Equiprobable residues',
#				    'Residue frequencies from input sequence',
#				    'Markov Chain (higher order dependencies)',
#				    'Lexicon partitioning',
#				    'Oligo frequencies from all non-coding regions'],
#			  -default=>$default{freq_estimate});
#print "<BR>";



#### table with all the statistics and thresholds
print "<BLOCKQUOTE>\n";
print $query->table({-border=>1,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			 [
			  $query->th([" <A HREF='help.oligo-analysis.html#return'>Return</A> ",
				   " <A HREF='help.oligo-analysis.html#thresholds'>Lower<BR>Threshold</A> ",
				   " <A HREF='help.oligo-analysis.html#thresholds'>Upper<BR>Threshold</A> "]),

			  ### occurrences
			  $query->td([$query->checkbox(-name=>'occ',
						    -checked=>$default{occ},
						    -label=>' Occurrences '),
				   $query->textfield(-name=>'occurrence_threshold',
						     -default=>$default{occurrence_threshold},
						     -size=>5),
				   '']),

			  ### binomial proba
			  $query->td([$query->checkbox(-name=>'proba',
						    -checked=>$default{proba},
						    -label=>' Binomial proba '),
				   '',
				   $query->textfield(-name=>'proba_occ_threshold',
						     -default=>$default{proba_occ_threshold},
						     -size=>5)]),

			  ### significance index
			  $query->td([$query->checkbox(-name=>'proba',
						    -checked=>$default{proba},
						    -label=>' Significance '),
				   $query->textfield(-name=>'occ_significance_threshold',
						     -default=>$default{occ_significance_threshold},
						     -size=>5),
				   '']),

			  ### Z-scores
			  $query->td([$query->checkbox(-name=>'zscore',
						    -checked=>$default{zscore},
						    -label=>' Z-scores '),
				   '',
				   '']),

			  ### frequencies
			  $query->td([$query->checkbox(-name=>'freq',
						    -checked=>$default{freq},
						    -label=>' Frequencies '),
				   '',
				   '']),

			  ### matching sequences
			  $query->td([$query->checkbox(-name=>'mseq',
						    -checked=>$default{mseq},
						    -label=>' Matching sequences '),
				   $query->textfield(-name=>'ms_threshold',
						     -default=>$default{ms_threshold},
						     -size=>5),
				   '']),

			  ### ratio
			  $query->td([$query->checkbox(-name=>'ratio',
						    -checked=>$default{ratio},
						    -label=>' Obs/exp ratio '),
				   '',
				   '']),

			  ### rank
			  $query->td([$query->checkbox(-name=>'rank',
						    -checked=>$default{rank},
						    -label=>' Rank '),
				   '',
				   ''])


			 ]
			)
		);
print "</BLOCKQUOTE>\n";


print "<HR width=550 align=left>\n";

### send results by e-mail or display on the browser
&SelectOutput;


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
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


