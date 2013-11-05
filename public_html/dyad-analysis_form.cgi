#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis_form.cgi,v 1.36 2012/07/08 20:22:17 jvanheld Exp $
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
##require "RSA.cgi.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{organism} = "Saccharomyces cerevisiae";
$default{taxon} = "Saccharomycetales";
$default{bg_method} = "monads";
#$default{bg_method} = "background";
$default{background} = "upstream-noorf";
$default{bg_level} = "organism";
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
#$default{lth_occ_sig} = "0";
$default{to_matrix} = '1';
$default{side} = 'over-represented';

## Return values and thresholds
$default{zscore} = '';
$default{lth_zscore} = 'none';
$default{uth_zscore} = 'none';

$default{rank} = 'checked';
$default{lth_rank} = "none";
$default{uth_rank} = "50";

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

&MatrixFromPatterns_defaults();

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("dyad-analysis", "form");

### head
print "<center>";
print "Analysis of spaced dyads in a set of DNA sequences\n<p/>";
print "Reference: <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/10734201'>van Helden, J., Rios, A. F. and Collado-Vides, J. (2000). Nucleic Acids Res 28, 1808-18.</a><p>\n";
print "</center>";
print "<u>Warning</u> !! For <b>vertebrate</b> genomes, analyses of complete promoters from <b>co-expressed gene groups</b> return <b>many false positive</b> (i.e. if you submit a random set of genes, you always get plenty of highly 'significant' motifs). This is likely to come from the heterogeneity of human sequences (mixtures of GC-rich and GC-poor promoters).
<br/>
However, analyses of <b>ChIP-seq peaks</b> return <b>very good</b> results. See the program <i><a href='peak-motifs_form.cgi'>peak-motifs</a></i>.
<p/>";

print "<hr>";
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

################################################################
## Dyad counting options
print "<b>Dyad counting mode</b><br>\n";

### Monad length
print "<B><A HREF='help.dyad-analysis.html#oligo_size'>Monad length</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_size',
			 -Values=>[3..3],
			 -default=>$default{oligo_size});

### spacing
print "<A HREF='help.dyad-analysis.html#spacing'><B>Spacing</B></A>&nbsp;\n";
print "&nbsp;", "from", "&nbsp;";
print $query->popup_menu(-name=>'spacing_from',
			 -Values=>[0..22],
			 -default=>$default{spacing_from});
print "&nbsp;", "to", "&nbsp;";
print $query->popup_menu(-name=>'spacing_to',
			 -Values=>[0..22],
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

&PrintDyadBackgroundOptions();

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

#print "<font color=red><B>Warning !</B> dyad-analysis is time-consuming, especially if you select a wide spacing range. If you don't obtain any result after 5 minutes, we recommend email output.</font><BR>\n";

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
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
print $query->hidden(-name=>'background',-default=>"upstream-noorf");
print $query->hidden(-name=>'bg_level',-default=>"organism");
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'gibbs_msps',-default=>"2");
print $query->hidden(-name=>'gibbs_flanks',-default=>"2");
#print $query->hidden(-name=>'title',-default=>'upstream sequences from the yeast GAL genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.dyad-analysis.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.dyad-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_dyad-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {

    print "<h4>Return</h4>\n";

    print "<BLOCKQUOTE>\n";
    print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
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
								  -size=>5),
						$query->popup_menu(-name=>'side',
								   -Values=>['over-represented','under-represented','both'],
								   -default=>$default{side})
					       ])
,

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

#### Convert patterns to matrix
&MatrixFromPatterns_print_form();
# print $query->checkbox(-name=>'to_matrix',
# 		       -checked=>$default{to_matrix},
# 		       -label=>'');
# print "&nbsp;Convert assembled patterns to Position-Specific Scoring Matrices (<font color=red>Can be time-consuming for large sequence files</font>).";
# print "<BR>";


    print "</blockquote>";
}
