#!/usr/bin/perl
#### this cgi script fills the HTML form for the program ORM
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
require "RSA2.cgi.lib";
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
#$default{sequence_type} = "dna";
$default{oligo_length} = 6;
$default{background} = "upstream-noorf";
$default{markov_order} = 2;
$default{strand} = "both strands";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{window_group} = 'checked';
$default{purge} = 'checked';
#$default{side} = 'over-represented';
$default{align} = 'right';
$default{freq_estimate} = "background";
$default{bg_level} = "organism";


$default{rank} = 'checked';
$default{lth_rank} = "none";
$default{uth_rank} = "20";

$default{lth_w_rank} = "none";
$default{uth_w_rank} = "1";


$default{occ} = 'checked';
$default{lth_occ} = "2";
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


$default{window_width} = '50';
$default{bg_window_width} = '500';

#$default{return}="fields";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 
$checked{$default{freq_estimate}} = "CHECKED";

### print the form ###
&RSA_header("ORM", "form");
print '<style><!-- textarea {height: 100px; width: 550px;}--></style>';
print "<CENTER>";
print "Analysis of oligonucleotide representation in a set of DNA sequences<P>\n";
print "</CENTER>";
print "<HR>";
print "<blockquote>";

#&ListParameters();

print $query->start_multipart_form(-action=>"ORM.cgi");

### Title
#print "<B><A HREF='help.ORM.html#title'>Title</A></B>&nbsp;\n";
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

			])

		 );

#### purge sequences
print $query->checkbox(-name=>'purge',
		       -checked=>$default{purge},
		       -label=>'');
print "&nbsp;<A HREF='help.ORM.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";

print "<HR width=550 align=left>\n";


### oligo size
print "<B><A HREF='help.ORM.html#oligo_length'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_length',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_length});

### prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "&nbsp;<A HREF='help.ORM.html#noov'><B>prevent overlapping matches</B></A>";
print "<BR>\n";

### strand ###
print "<B><A HREF='help.ORM.html#count_strands'>Count on</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

### align ###
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.ORM.html#align'>Align</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'align',
			 -Values=>['right',
 				  'left'],
 			 -default=>$default{align});


#### group patterns by pairs of reverse complement
#print $query->checkbox(-name=>'grouprc',
#		       -checked=>$default{grouprc},
#		       -label=>'');
#print "&nbsp;<A HREF='help.ORM.html#grouprc'><B>return reverse complements together in the output</B></A>";
#print "<BR>";
print '<br />';
print '<b><a href="help.ORM.html#window_width">Window width</a></b>';
print $query->textfield(-name=>'window_width',
			-default=>$default{window_width},
			-size=>5);
print '<b><a href="help.ORM.html#window_width">group windows</a></b>';
print $query->checkbox(-name=>'window_group',
		       -checked=>$default{window_group},
		       -label=>'');


print '&nbsp;&nbsp;&nbsp;&nbsp;<b><a href="help.ORM.html#bgwindow_width">Background Window width</a></b>';
print $query->textfield(-name=>'bg_window_width',
			-default=>$default{bg_window_width},
			-size=>5);




print "<HR width=550 align=left>\n";

################################################################
## Background model
&PrintOligoBackgroundOptions();

################################################################


#### estimation of expected frequencies
#print "<A HREF='help.ORM.html#exp_freq'><B>Expected frequency calibration</B></A>&nbsp;<br />";
#### pre-defined background frequencies
#print ( "<INPUT TYPE='radio' NAME='freq_estimate' VALUE='background' $checked{background}>", 
#	"Predefined background frequencies");
#print "<ul>";
#print ( "<a href='help.ORM.html#background'>Background model</a> &nbsp;&nbsp;&nbsp;&nbsp;", 
#	$query->popup_menu(-name=>'background',
#			   -Values=>["upstream","upstream-noorf","intergenic"],
#			   -default=>$default{background}));
#	
#print "<br>", &OrganismPopUpString();
#print "</ul>";


#### Markov chain model
#print ("<p><INPUT TYPE='radio' NAME='freq_estimate' VALUE='Markov Chain (higher order dependencies)' $checked{'Markov Chain (higher order dependencies)'}>", 
#       "Markov Chain (higher order dependencies)");

#print "order &nbsp;";
#print $query->textfield(-name=>'markov_order',
#			-default=>$default{markov_order},
#			-size=>5);

#print '</p>';

#### Bernouilli model
#print "<p><INPUT TYPE='radio' NAME='freq_estimate' VALUE='Residue frequencies from input sequence' $checked{'Residue frequencies from input sequence'}>Residue frequencies from input sequence</p>";

#### custom expected frequency file
#print "<p><INPUT TYPE='radio' NAME='freq_estimate' VALUE='file_upload' $checked{'file_upload'}><a href='help.ORM.html#upload_freq_file'>Upload your own expected frequency file</a><br />";
#print $query->filefield(-name=>'upload_freq_file',
#			-default=>'starting value',
#			-size=>30,
#			-maxlength=>200);
#print '</p>';

################################################################
print "<HR width=550 align=left>\n";

#print "<A HREF='help.ORM.html#exp_freq'><B>Expected frequency</B></A>&nbsp;";
#print $query->radio_group(-name=>'freq_estimate',
#			  -Values=>['Equiprobable residues',
#				    'Residue frequencies from input sequence',
#				    'Markov Chain (higher order dependencies)',
#				    'Lexicon partitioning',
#				    'Oligo frequencies from all intergenic regions'],
#			  -default=>$default{freq_estimate});
#print "<BR>";

&ReturnTable();

print "<HR width=550 align=left>\n";

### send results by email or display on the browser
&SelectOutput();


### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"ORM_form.cgi");
$demo_sequence = ">MET8	YBR213W; upstream from -463 to -1; size: 463; location: NC_001134.7 649900 650362 D; upstream neighbour: YBR212W (distance: 463)
TTACAAAAGACAAAAAAAGAAAATTTTAATCTTGTCCGCAGTTTTATCTGCGTCTCTACG
TTCTTACGTTTCTTCTATTAATGCCATTTCAGTTACAACCTAGTCAATTGTCGATCCATA
ATTCTAATCAAATTTGTTTTTCCTCTATACTACCTATCTATTTTTATCTATCTAAGTACA
TTTATTTACTCAAACAGTTCCGTTTCAAAGTGTTTTATATTAACTATATATGCGAAAAGC
TGGCGTCATAATTTCACGTGTTATAATAGCCATGCTGACGGAAAAAAAATGTGAAAATCG
CTACAAAGTCCGATGACTACGGGCAGTAGCATGTAAATGATGGACACACACACACACATA
TATATATATATACATTTACTTCAATAAAAGGCTGTGCCAGACATTTTTGCCATACATTGT
TCATGAAGTGTGCAAAATAAGAGAGTGTATAATAGGATAAAAA
>MET32	YDR253C; upstream from -547 to -1; size: 547; location: NC_001136.8 964562 965108 R; upstream neighbour: YDR254W (distance: 547)
TAATTGCTACTCAAATATACTAGTCAAAGATAGTATCCACCAAAATCTTTCCCCGCTAAA
ATAACGCCAGATGCTTTCTATGCTTCTAATCTTTTACCATTTACCTTTGTTTATTTCAAT
ATAAACTTTAATTTACAGTCCCTATCTATTGCCCGACTGGACTAACATGCACGTGACATT
TTGTGATGGTTTTTCGTCCCTTACTTAGTACGCTTAGTACGCCACAGTTTATATTTTCTT
GACAATAATAAAGAACCTGATTGTGGGTTAGAACTTGCTATACTTTTAGTTTAAAATAAG
CAGGAAATAATCTTGAGTTCTGTATCATTATTATAAATAAAACTATATTTGTTCTCTTTG
TCGCCCTCGGAACTTTCCTCATTACATTGACGAGGTATATATAGATATAGTAGATATACA
TATCTATCCATGGTATATATGTATGCATCTGGATAATTGAATAGGGTTTCATGTCATATG
CCAAGAATTTGTTAATAATATAGTGGAAAAAAGTCAAGAGGTATTATAAATTTCAAAAAA
GTACCAA
>MET18	YIL128W; upstream from -568 to -1; size: 568; location: NC_001141.1 113238 113805 D; upstream neighbour: YIL129C (distance: 568)
TTTGATGTATAACAAAACTAAAAAGGGTTATTAAAATGGGAACACAACAAACAACCAGAA
TTTTCACACTTTAACCAGTCACGTCCTATTATGAAGACCTAAATCCACATTTGCTTTCTC
TCTTCATTTGCCTAATCCTTTATCCCAATTTCTACAGTTCTATATGTATTTTCCTGTGTG
GCTGTCGTTTCGTGGTTAGTGATACAACCATAACGATTCAACCAACTCCCAATGTATGTG
ATGTTGATACCGCTAATTTGGAAGGGATGGTATACTCTAGGTGACCTCAATGAGTCAAAG
AGAGCTAGGACATACTTCGAGATAGGTAATACCACTTTGCAGCTTCTTTTTAGGCCTTCA
TGAGTGAGTAGCCAAGAAAAAGTTAAAAGCGGGTAATAGGTATGAATTTTTCAAATACTG
AAATTTGGTTTAGTTATTTAAGTGAATTGTAGATTATGTACATTTTACGTGCAATGAAGG
AGTCACCTCTATGATCATCTAGTTATTAGCTGTTAGTTTTCATTGAACTTGTTTTAACTG
GGAAAAAGCGGAACAATTGGGCCTTACA
>MET30	YIL046W; upstream from -177 to -1; size: 177; location: NC_001141.1 268473 268649 D; upstream neighbour: YIL046W-A (distance: 177)
CACGTGATCGGGAAGCCACAGTTTGCGCGGAGATATTTTATTTTTTTTCATCAGCGTAAG
AAGAAAGCAACCTTGCAGTCTGTATCGTAAGAGAAGACTGCAGTTAAAGAAGTTTAGAGA
AGAGGCTTGAGTATCGGTAAAGGGGTGTGTGTTTGGTGATTTATAAAGGAGAAGGGC
>MET28	YIR017C; upstream from -489 to -1; size: 489; location: NC_001141.1 384117 384605 R; upstream neighbour: YIR018W (distance: 489)
GACTGTGATAATATGCTAGTTACACTGTTTATGTTGTGTGAACTTGTTGTAATATGGTTA
ACTTCACTTTCAGTGATTGATATGATAGCGACATCACTGCCGTGCAAAAAGACCATTCCA
TTACTGCACCTTTTTGTCCTTTTCCGTGGAATAAAAGTTCACTCGTCAGTTCCATGCATT
CTGGAAAAAAATGATCTGAAAGATGCCACAGTTGTGGGGCCCGCCCGGCCCAATAGGTAA
ACTAAAATACAATAGAAGGGGTACTGAGTGCACGTGACTTATTTTTTTTTTTTGGTTTTA
GGTTTCGCTTTTTTCACCTTTTTCTACTTTCTAACACCACAGTTTTGGGCGGGAAGCGGA
AACGCCATAGTTGTAGGTCACTGGCGTGAGTCAAGGCCGGGCAGCCAATGACTAAGAACA
CGAGGTAACTTGAATTTAACTATTTATAACCAGTGGTAGTTACGAAGACAAATTGTTTTG
TTCGTCAAT
>MET6	YER091C; upstream from -687 to -1; size: 687; location: NC_001137.2 342164 342850 R; upstream neighbour: YER092W (distance: 687)
TTTTTTCTTGTTTTATAATCAGTCAAGTATTGGTTTCCCACAGCCATTCAACTCAGGTTC
ATCATCTTTTTCGCTTCCAAAAATGCAGTTGATTTCACACAATTTTTCATGAACCAGGGT
CCCGCACTCCGGGTAAAGGACCATCACGCCACATCACGTGCACATTACTAGTAAAAGCCA
CAGGAAATATTTCACGTGACTTACAAACAGAGTCGTACGTCAGGACCGGAGTCAGGTGAA
AAAATGTGGGCCGGTAAAGGGAAAAAACCAGAAACGGGACTACTATCGAACTCGTTTAGT
CGCGAACGTGCAAAAGGCCAATATTTTTCGCTAGAGTCATCGCAGTCATGGCAGCTCTTT
CGCTCTATCTCCCGGTCGCAAAACTGTGGTAGTCATAGCTCGTTCTGCTCAATTGAGAAC
TGTGAATGTGAATATGGAACAAATGCGATAGATGCACTAATTTAAGGGAAGCTAGCTAGT
TTTCCCAACTGCGAAAGAAAAAAAGGAAAGAAAAAAAAATTCTATATAAGTGATAGATAT
TTCCATCTTTACTAGCATTAGTTTCTCTTTTACGTATTCAATATTTTTGTTAAACTCTTC
CTTTATCATAAAAAAGCAAGCATCTAAGAGCATTGACAACACTCTAAGAAACAAAATACC
AATATAATTTCAAAGTACATATCAAAA
>MET10	YFR030W; upstream from -338 to -1; size: 338; location: NC_001138.4 212962 213299 D; upstream neighbour: YFR029W (distance: 338)
TGCATCTAAATATATACGTATGTTTAAGGTTCTGGTATACAGGTATTAAAAGAAAACACT
ATCAACATTCCCAATAAGATATACCACACCACGTGAGCTTATAGAAGCACGTGACCACAA
TTCACCCCACAGGTGTGGCTTTTTTGGTGCCGTAGAAAAGACTCATTCATGAATCGTCGG
AAACCCATAGTCATCTTCGAGCAAAAGGTATATATAAGCAACAGAGGGCAGTAGTTCTCG
AGACCACCATCTTTTGATTGGAAATAGTTTCGTTTAGATGGGGTGCACATAGTTTTTTTC
AACTGCTTTTCCTCGAGGTCACCCAAATATACAACGAG
>MET13	YGL125W; upstream from -380 to -1; size: 380; location: NC_001139.7 272146 272525 D; upstream neighbour: YGL126W (distance: 380)
CTCAGGAAAAGTTGGCGATAGACCACGAGCGACTGAAAAAATAACAGCGACTTTTCTCCC
GGTAGCGGGCCGTCGTTTAGTCATTCTATCCCTCGGATTATAGACTGTGAATATTGCATA
TGCAACTTTGACTCAAATTTTTCCAAAATTTGATATATATATATATATATATATGTTTGT
ATGTATATATATATATACGTATATATATCATATATACGAAAAGTAGAAAAAAAAAGGTGA
TATTTCGCTCGTGGAAAAGCTAATGCCACAGCTTGTGTTTCGTGTAGTTTGCCTTGCTCC
CCTTGATTGAAATAGTCTCCCTAAACTAAAGTTATCAGCAAACAGAACCACCACAGTTAC
TACTACAACCACATCGCAAT
>MET3	YJR010W; upstream from -800 to -1; size: 800; location: NC_001142.6 455354 456153 D; upstream neighbour: YJR009C (distance: 1557)
AAGAGTACAATTTATAAATTAATGAAAACACAGAAGTATTTAGATCGGCTCAAATGTTTT
TGGACATTAAAAGATCTTGAAACTGAGTAAGATGCTCAGAATACCCGTCAAGATAAGAGT
ATAATGTAGAGTAATATACCAAGTATTCAGCATATTCTCCTCTTCTTTTGTATAAATCAC
GGAAGGGATGATTTATAAGAAAAATGAATACTATTACACTTCATTTACCACCCTCTGATC
TAGATTTTCCAACGATATGTACGTAGTGGTATAAGGTGAGGGGGTCCACAGATATAACAT
CGTTTAATTTAGTACTAACAGAGACTTTTGTCACAACTACATATAAGTGTACAAATATAG
TACAGATATGACACACTTGTAGCGCCAACGCGCATCCTACGGATTGCTGACAGAAAAAAA
GGTCACGTGACCAGAAAAGTCACGTGTAATTTTGTAACTCACCGCATTCTAGCGGTCCCT
GTCGTGCACACTGCACTCAACACCATAAACCTTAGCAACCTCCAAAGGAAATCACCGTAT
AACAAAGCCACAGTTTTACAACTTAGTCTCTTATGAAGTTACTTACCAATGAGAAATAGA
GGCTCTTTCTCGAGAAATATGAATATGGATATATATATATATATATATATATATATATAT
ATATATGTAAACTTGGTTCTTTTTTAGCTTGTGATCTCTAGCTTGGGTCTCTCTCTGTCG
TAACAGTTGTGATATCGTTTCTTAACAATTGAAAAGGAACTAAGAAAGTATAATAATAAC
AAGAATAAAGTATAATTAAC
>MET14	YKL001C; upstream from -800 to -1; size: 800; location: NC_001143.7 439029 439828 R; upstream neighbour: YKR001C (distance: 1222)
TATTTTTTTAATTACATAATCATAAAAATAAATGTTCATGATTTCCGAACGTATAAAATA
AGAATGTTACGAGAATTTGTTTTCTTGGTAATTAAAATAATCAAATACACATAGAAAGGA
GAGTAAACTGCTTCCTCTGTATAAATCAAAGCAAAATTGTAAATAGCGTTGACAAGTGAT
TACAGAAGTTAGGTGAGGTTAATTACCAATTTCTTTTTTTAAAATTGGTGAAATAAGATT
ACGTTTAAAGGAGCATTAACAGGTTTACTCATAACAATCATTTTCAAATTTCCCTATGCA
TGTTTAGAGCAAGCGCCTTTGTGAGCCCTCCCGGTTACGACGCCTTGGCAATGTAGCAGA
TAACTCTGCACTTCTAGAATCATTCCACTACGACATTTGGCTCATCACCAGCTCGCGAGA
AATGTAAATAAGCCAACAACCAAGAATGCGTAACATTAAAGAATACAGTTGCTTTCATTT
CGGCGTGATGGTACGGCACCCACGGTACCTTACATTATTCTCGAAAAATAGCTGCACGCT
TTTCCAGGAATAAAAGACCGTGCCACTAATTTCACGTGATCAATATATTTACAAGCCACC
TCAAAAAATGTGGCAATGGAGAAGAGGATGAACGACTCAATATGACTTCAACTTCATGAA
TTTGTCAAAATATCTATATAAGATGCAAAATTTCTATACAACATCAGTTGCGTATCCGTT
AATGTCGTTCATTTTCTCTCTTTGTTCGAACTTGACATCAAGAAAAGTTGGAATTATTTC
TCCAAGCACACTGTACACCA
>MET1	YKR069W; upstream from -702 to -1; size: 702; location: NC_001143.7 570552 571253 D; upstream neighbour: YKR068C (distance: 702)
TTTTGACCCAGTTTTGGTTATCAATGAACACTTGAAGCTTTACTCTGCATTCCCATCTCT
ATAGCTATGGGTAATCACAGCTACGATCACTTACTCTGTTATTATTATATTAAGTTCAAT
GTTGGCCAAACCGGGTAACATGTAACACTTTCAGGTTGGCCTTACCTTTGGCTTGGAGTT
TCGCAAGTTTTCAAATTTTTGGCTCCTGCTGTCAAGGTGCATAGAATAGCGCTTATTTAT
CTATTTATATCCAAGATGTACAATCCTCGTTCTCTGAGTCCAACATATTTGCTCGCAACT
GTAGAAATCACAACTACAGCAACAGTAAAGATATCATTTTCTATTTTCGTTATTGGTTTC
TCGACCTTTTTATATACGATACGTCAAACTTGAATCATTTTATACGTTTTTCTCTTTCTA
GAAATGCCATTATGCACGTGACATTACAAATTGTGGTGAAAAAAGGCTCTCATAATAAAC
TGTGAACGGACTCATAATGAAATTTGCTTCACTATGTGAATCATCGCTAATAAACTCGCT
ACAAAAGTCGAGTATGCTTAAGTCAAAAAAATGATATATATATATAATTTACTTATGTGT
TTCTGCAAAGTTGTAGGCTTCATTTAGAATTGCTCAGATATTCCATCCCAATTAAAAAAA
GCACGGATAGAGTGATAAATAAACTAAGAAAATTTCAAAAGA
>MET17	YLR303W; upstream from -800 to -1; size: 800; location: NC_001144.4 731744 732543 D; upstream neighbour: YLR301W (distance: 982)
TATACTAGAAGTTCTCCTCGAGGATTTAGGAATCCATAAAAGGGAATCTGCAATTCTACA
CAATTCTATAAATATTATTATCATCGTTTTATATGTTAATATTCATTGATCCTATTACAT
TATCAATCCTTGCGTTTCAGCTTCCACTAATTTAGATGACTATTTCTCATCATTTGCGTC
ATCTTCTAACACCGTATATGATAATATACTAGTAACGTAAATACTAGTTAGTAGATGATA
GTTGATTTTTATTCCAACACTAAGAAATAATTTCGCCATTTCTTGAATGTATTTAAAGAT
ATTTAATGCTATAATAGACATTTAAATCCAATTCTTCCAACATACAATGGGAGTTTGGCC
GAGTGGTTTAAGGCGTCAGATTTAGGTGGATTTAACCTCTAAAATCTCTGATATCTTCGG
ATGCAAGGGTTCGAATCCCTTAGCTCTCATTATTTTTTGCTTTTTCTCTTGAGGTCACAT
GATCGCAAAATGGCAAATGGCACGTGAAGCTGTCGATATTGGGGAACTGTGGTGGTTGGC
AAATGACTAATTAAGTTAGTCAAGGCGCCATCCTCATGAAAACTGTGTAACATAATAACC
GAAGTGTCGAAAAGGTGGCACCTTGTCCAATTGAACACGCTCGATGAAAAAAATAAGATA
TATATAAGGTTAAGTAAAGCGTCTGTTAGAAAGGAAGTTTTTCCTTTTTCTTGCTCTCTT
GTCTTTTCATCTACTATTTCCTTCGTGTAATACAGGGTCGTCAGATACATAGATACAATT
CTATTACCCCCATCCATACA
>MET2	YNL277W; upstream from -481 to -1; size: 481; location: NC_001146.5 116868 117348 D; upstream neighbour: YNL277W-A (distance: 481)
GCAGTATAAATTGTACTTCAAAGCACTAGTCATGAAAAACGCTTACATTAGTTCAGTTTG
TCAAGGTTATGCTATTACTTGTACTTATTTCTTGCTATTGTTAGTGGCTCCCCACATTGA
CGTATTTTCACGTGATGCGCCTCACTGCGGAAGGCGCCACACATTGCCTGCAAAAAATTG
TGGATGCACTCATTTGATAGTAAACTAAGTCATGTTAATCGTTTGGATTTGGCACACACC
CACAAATATACACATTACATATATATATATATTCAAAATACAGCTGCGTCCAATAGATGA
GCTTCCGCTTCGTTGTACAACCTACCTGCTATCTTGTTCACGGATATTTCTTGCTTTTAA
TAAACAAAAGTAACTCTAGAACAGTCAAGTCTTCGATAATTTTTTTAGTCACAGGGTCCG
TCTAAAGTTTCTCTTTATTTGGAATAATAGAAAAGAAAGAAAAAAACGTAGTATAAAAGG
A
>MET4	YNL103W; upstream from -800 to -1; size: 800; location: NC_001146.5 426937 427736 D; upstream neighbour: YNL104C (distance: 980)
AGAGGCTGACCCAAGAGGAGAAAACATCGAACCTACGCGACTCCATAGCGCACATCTCCC
ATGCGCCCGTGCTTATATATATATAAATATACATACACACATACATGCACGCACATACAT
GCACGCGCATGACAGTGACTGGCCGTTACTGTACAATTTTTTCAGCCAAGTATGACACAC
ATTCAACTCAGCTTTTCTGAGGCCTTCTTTCTTTTCCTGCGCGTCGGTAGAGCGATGACT
AACCTACTACTGTCTCAGAGCCGGTCCCGCTCCGGTAGCAATCCTGGGGCTGGTCATAAC
AGCCGAGTGGAAGTGTCAAAGCGGAGAACAGAAGCATAAGCTCAATCGCTGGACATACGG
ATGCTTATATACGTCTTATTGTCGTTGAAAAATATCGAATTTTTTACTTCATTTATCGAG
GCTTCTTCGAGCACTTTTCCGCTATGGCTTTTTCCCCGTTTCCTTTTAATCACGTGCGCG
GGTAGCACCCGGCACACAGCTGGTGTCTCGTCGCACATGCTATTGTGTGTCATCGGGCCA
CACAAGCATATTGCTTGAATTTTCTTTCATCGTTCAACTTAAATCCACCCAATCTAGATG
TAGCCGTAGCATGTAATAACGTATATCCTTGTTTACATGCATCTGTGCCAGGTGAAACGG
TCTGTTTGAACGCACATCATTTCATATATTAGTCAACTCCTGAAGGTCTTCTTGCCTGTC
CGTCAACTGTTTAGACAGACTCTCGTCAATAAAGCGCACTTCTGATAAGCACTTTTATTC
CTTTTTTTCCACTGTGAACG
>MET22	YOL064C; upstream from -215 to -1; size: 215; location: NC_001147.5 207177 207391 R; upstream neighbour: YOL063C (distance: 215)
CATATTTTGACATTTACATAGCCATCTATATAATAATCCTTCCTTCATTGAAATGCGCGA
ATGACTCAGACGAGCAATATCACGTGTTGCGATTTTACTTTCAGTTTGCAAAGAAGAAAC
CATGCCACATATAAAGAACGTTTGCACTTCTCTTTAATTTATTAGTTAGTAAGTAAGAAG
TTTAAAGACAACTCAGAAGACATCAGCACTTTACT
>MET7	YOR241W; upstream from -250 to -1; size: 250; location: NC_001147.5 786746 786995 D; upstream neighbour: YOR239W (distance: 250)
GAAGTTCTGAGACAAGTACCACCTCCTCTCTCATCATAAAACAAGTAAAAGTTTTCTCGT
CGCGCATATTATTTTGGTGATTGATTGTTTTTTCCTCCGATATCATCACTTATTACCTGT
AATTTTATCTTTTTCTACCCCATAGAATTCGTCTTATAAGTCTATACCCTCAAAACTATC
TATCATTTTAATATTATCTGTCGCTTTAATTGTCTTATTTCTGAAGCTCACTGAAGAACA
TTGCTTTATT
>MET31	YPL038W; upstream from -161 to -1; size: 161; location: NC_001148.3 480371 480531 D; upstream neighbour: YPL038W-A (distance: 161)
ATTTTAGTCTAAAAATTTTGCTAGCCCATCAATTTTTTTTTTGTTCTAATGCAAAATATA
ACATGGGTAAGAAAAAGAAAAAGCCGTTCCTCAGTACGTAAAGAGATTTGATCATTAACA
AGTTGGGCTCAATATACACAGTCGATAGTCTATATGTGCAT
>MET12	YPL023C; upstream from -384 to -1; size: 384; location: NC_001148.3 506311 506694 R; upstream neighbour: YPL022W (distance: 384)
CTGGAAAGATATTTTCAACAGGATAGTGCAATATTATTTTTACACATTTAGCAAATGCTC
TACCAACGTCCTGAACCCTCCAGTACACCTGTGCTCTCTCCTCTCCTATGCGCCACGCAG
ACAAAAGTTTACTCGTCCCGACTTTTTTTTTTCATTAGACGCGATATTGACTGTGGCTAT
AGCTTACTCCAGGAGATAAGCGAGTAAAGCTTTTTCTAACTTCAATGATGAAGAAAAGTC
GCAAAATAAAGGCAAACAGAGAACACTTCAGGTTGTTGGTTACATTGGAAGAGGGACTTA
AGCTCTCACATCATCTATTTTGTTTCAAGTTCGTACATTTTTTGAAGCGTGTTGGACGGG
ACAGGTTGATTACATTTTTTAAAC
>MET16	YPR167C; upstream from -443 to -1; size: 443; location: NC_001148.3 877629 878071 R; upstream neighbour: YPR168W (distance: 443)
CTTATCGGTTTATTTTTCTATATATTTGCCTCTTTCTCAAACAGGAGTTAGTAGTTAAAA
GTACGAAGTTCTTGTTCTTTAATGCGCGCTGACAAAAGAATTGGATAAAAGAGAATGGTG
GGGGGACAAGAAGGAAATTTGTCCTAGTTTAACATGAATGGCATCTTGTTACCGGGTGGA
CATCACCTATTGATTCTAAATATCTTTACGGTTTATCATACTGTTCTTTATTCCGTCGTT
ATTCTTTTTATTTTTATCATCATTTCACGTGGCTAGTAAAAGAAAAGCCACAACATGACT
CAGCAAATCTCGACAAAGTAAAAGCTCATAGAGATAGTATTATATTGATATAAAAAAAGT
ATACTGTACTGTTTGTAACCTTTTCAATGCTTTAAGATCAAAACTAAGGCCAGCAAAGGT
ATCAACCCATAGCAACTCATAAA
";
print "<TD><B>";
print $query->hidden(-name=>'lth_occ',-default=>'2');
print $query->hidden(-name=>'window_width',-default=>'800');
print $query->hidden(-name=>'bg_level',-default=>"organism");
print $query->hidden(-name=>'bg_window_width',-default=>'800');
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'title',-default=>'upstream sequences from the yeast MET genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.ORM.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.ORM.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_ORM.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {

print "<b>Thresholds</b><br />\n";

#print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>CENTER},
			 [
			  $query->th([" <A HREF='help.ORM.html#return_fields'>Fields</A> ",
				   " <A HREF='help.ORM.html#thresholds'>Lower<BR>Threshold</A> ",
				   " <A HREF='help.ORM.html#thresholds'>Upper<BR>Threshold</A> "]),

			  ### occurrences
			  $query->td(['Occurrences',
				      $query->textfield(-name=>'lth_occ',
							-default=>$default{lth_occ},
							-size=>5),
				      $query->textfield(-name=>'uth_occ',
							-default=>$default{uth_occ},
							-size=>5)
				   ]),

			  ### binomial proba
#			  $query->td(["Probability",
#				      $query->textfield(-name=>'lth_occ_P',
#							-default=>$default{lth_occ_P},
#							-size=>5),
#				      $query->textfield(-name=>'uth_occ_P',
#							-default=>$default{uth_occ_P},
#							-size=>5),
#				      $query->popup_menu(-name=>'side',
#							 -Values=>['over-represented','under-represented','both'],
#							 -default=>$default{side})
#				     ]),

			  ### binomial E-value
			  $query->td(["E-value",
#				      $query->checkbox(-name=>'proba',
#						       -checked=>$default{proba},
#						       -label=>' Binomial E-value '),
				      $query->textfield(-name=>'lth_occ_E',
							-default=>$default{lth_occ_E},
							-size=>5),
				      $query->textfield(-name=>'uth_occ_E',
							-default=>$default{uth_occ_E},
							-size=>5),
				     ]),

			  ### significance index
			  $query->td(["Significance",
#				      $query->checkbox(-name=>'proba',
#						    -checked=>$default{proba},
#						    -label=>' Significance '),
				   $query->textfield(-name=>'lth_occ_sig',
						     -default=>$default{lth_occ_sig},
						     -size=>5),
				   $query->textfield(-name=>'uth_occ_sig',
						     -default=>$default{uth_occ_sig},
						     -size=>5)
				   ]),


			  ### frequencies
#			  $query->td(["Frequencies",
#				   $query->textfield(-name=>'lth_observed_freq',
#						     -default=>$default{lth_observed_freq},
#						     -size=>5),
#				   $query->textfield(-name=>'uth_observed_freq',
#						     -default=>$default{uth_observed_freq},
#						     -size=>5)
#				      ]),


			  ### rank
			  $query->td(["Rank",
				      $query->textfield(-name=>'lth_rank',
							-default=>$default{lth_rank},
							-size=>5),
				      $query->textfield(-name=>'uth_rank',
							-default=>$default{uth_rank},
							-size=>5)
				      ]),
			  ### w_rank
			  $query->td(["Window rank",
				      $query->textfield(-name=>'lth_w_rank',
							-default=>$default{lth_w_rank},
							-size=>5),
				      $query->textfield(-name=>'uth_w_rank',
							-default=>$default{uth_w_rank},
							-size=>5)
				      ]),



			  
			  ]
			       )
		    );
#print "</BLOCKQUOTE>\n";

}
