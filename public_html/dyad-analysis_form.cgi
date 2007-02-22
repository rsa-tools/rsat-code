#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis_form.cgi,v 1.13 2007/02/22 09:40:51 jvanheld Exp $
#
# Time-stamp: <2003-07-11 15:08:24 jvanheld>
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
$default{freq_estimate} = "background";
#$default{title} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_file} = "";
$default{oligo_size} = 3;
$default{spacing_from} = 0;
$default{spacing_to} = 20;
$default{strand} = "both strands";
$default{noov} = 'checked';
$default{purge} = 'checked';
$default{dyad_type} = "any dyad";
$default{exp_freq} = "background";
$default{upload_freq_file} = "";
$default{background} = "upstream-noorf";
#$default{lth_occ_sig} = "0";

## Return values and thresholds
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

$default{proba} = 'checked';
$default{lth_occ_P} = "none";
$default{uth_occ_P} = "none";

$default{eval} = 'checked';
$default{lth_occ_E} = "none";
$default{uth_occ_E} = "none";

$default{lth_occ_sig} = "0";
$default{uth_occ_sig} = "none";

$default{freq} = '';
$default{lth_observed_freq} = "none";
$default{uth_observed_freq} = "none";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 
$checked{$default{freq_estimate}} = "CHECKED";

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



&DisplaySequenceChoice();


#### purge sequences
print $query->checkbox(-name=>'purge',
  		       -checked=>$default{purge},
  		       -label=>'');
print "&nbsp;<A HREF='help.dyad-analysis.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";
print "<HR width=550 align=left>\n";

### oligo size
print "<B><A HREF='help.dyad-analysis.html#oligo_size'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_size',
			 -Values=>[3..3],
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
print "<A HREF='help.dyad-analysis.html#exp_freq'><B>Expected frequency calibration</B></A>&nbsp;<p>";
print ( "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='background' $checked{background}>", 
	"Background model");
print "<ul>";
print ( "<a href='help.dyad-analysis.html#background'>Sequence type</a> &nbsp;&nbsp;&nbsp;&nbsp;", 
	$query->popup_menu(-name=>'background',
			   -Values=>["upstream","upstream-noorf","intergenic"],
			   -default=>$default{background}));

print "<br>", &OrganismPopUpString();
print "</ul>";

print ( "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='monads' $checked{monads}>", 
	"Monad frequencies in the input sequence");

print "<BR>\n";

#### custom expected frequency file
print "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='file_upload' $checked{file_upload}><a href='help.oligo-analysis.html#upload_freq_file'>Upload your own expected frequency file</a><BR>";

print $query->filefield(-name=>'upload_freq_file',
			-default=>'starting value',
			-size=>30,
			-maxlength=>200);
print "<p>";


# print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
# 		    $query->Tr($query->td("<A HREF='help.dyad-analysis.html#exp_freq'><B>Expected frequency calibration</B></A>&nbsp;<BR>")),
# 		    $query->Tr($query->td(["<INPUT TYPE='radio' NAME='exp_freq' VALUE='dyad freq from intergenic sequences' CHECKED>Dyad frequencies from all intergenic regions<BR>",
# 					   &OrganismPopUpString])),
# 		    $query->Tr($query->td([
# 					   "<INPUT TYPE='radio' NAME='exp_freq' VALUE='monad (word) freq in the input sequences.'>Monad (word) frequencies from the input sequences<BR>",
# 					   ])),
# 		    );

 print "<HR width=550 align=left>\n";

# ### significance threshold
# print "<B><A HREF='help.dyad-analysis.html#threshold'>\n";
# print "Threshold of significance</A> >= \n";
# print $query->textfield(-name=>'lth_occ_sig',
# 		  -default=>$default{lth_occ_sig},
# 		  -size=>5);
# print "<BR>\n";

&ReturnTable();



### send results by email or display on the browser
print "<HR width=550 align=left>\n";

&SelectOutput();

print "<font color=red><B>Warning !</B> dyad-analysis is time-consuming, especially if you select a wide spacing range. If you don't obtain any result after 5 minutes, we recommend email output.</font><BR>\n";

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
#print $query->hidden(-name=>'spacing_from',-default=>"8");
#print $query->hidden(-name=>'spacing_to',-default=>"12");
print $query->hidden(-name=>'background',-default=>"upstream");
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
#print $query->hidden(-name=>'title',-default=>'upstream sequences from the yeast GAL genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.dyad-analysis.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.dyad-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_dyad-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {

    print "<h3>Return</h3>\n";

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
						$query->textfield(-name=>'lth_occ_P',
								  -default=>$default{lth_occ_P},
								  -size=>5),
						$query->textfield(-name=>'uth_occ_P',
								  -default=>$default{uth_occ_P},
								  -size=>5)]),
				    ### binomial E-value
				    $query->td([$query->checkbox(-name=>'eval',
								 -checked=>$default{eval},
								 -label=>' Binomial E-value '),
						$query->textfield(-name=>'lth_occ_E',
								  -default=>$default{lth_occ_E},
								  -size=>5),
						$query->textfield(-name=>'uth_occ_E',
								  -default=>$default{uth_occ_E},
								  -size=>5),
					       ]),

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

    print "</blockquote>";
}
