#!/usr/bin/perl
#### this cgi script fills the HTML form for the program oligo-analysis
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
$default{organism} = "Saccharomyces cerevisiae";
$default{taxon} = "Saccharomycetales";
$default{title} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_freq_file} = "";
$default{sequence_type} = "dna";
#$default{oligo_length} = 6;
$default{oligo_length1}="";
$default{oligo_length2}="";
$default{oligo_length3}="";
$default{oligo_length4}="";
$default{oligo_length5}="";
$default{oligo_length6}="checked";
$default{oligo_length7}="checked";
$default{oligo_length8}="checked";
$default{oligo_length9}="";

$default{bg_method} = "Markov model (higher order dependencies)";
##$default{bg_method} = "background";
$default{background} = "upstream-noorf";
$default{bg_level} = "organism";
$default{markov_order} = 2;
$default{pseudo_freq} = "0.01";
$default{strand} = "both strands";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{purge} = 'checked';
$default{side} = 'over-represented';

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
$default{lth_occ} = "none";
$default{uth_occ} = "none";

$default{mseq} = '';
$default{uth_mseq} = "none";
$default{lth_mseq} = "none";

$default{proba} = 'checked';
$default{lth_occ_P} = "none";
$default{uth_occ_P} = "none";

#$default{eval} = 'checked';
$default{lth_occ_E} = "none";
$default{uth_occ_E} = "none";

$default{occ_FWER} = '';
$default{lth_occ_FWER} = "none";
$default{uth_occ_FWER} = "none";

$default{lth_occ_sig} = "0";
$default{uth_occ_sig} = "none";

$default{freq} = '';
$default{lth_observed_freq} = "none";
$default{uth_observed_freq} = "none";

$default{return}="fields";

&MatrixFromPatterns_defaults();

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("oligo-analysis", "form");
print "<center>";
print "Analysis of oligomer occurrences in nucleotidic of peptidic sequences<P>\n";
print "Reference: <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/9719638'>van Helden, J., Andr&eacute;, B. and Collado-Vides, J. (1998). . J Mol Biol 281, 827-42.</a><p>";
print "</center>";
print "<u>Warning</u> !! For <b>vertebrate</b> genomes, analyses of complete promoters from <b>co-expressed gene groups</b> return <b>many false positive</b> (i.e. if you submit a random set of genes, you always get plenty of highly 'significant' motifs). This is likely to come from the heterogeneity of human sequences (mixtures of GC-rich and GC-poor promoters).
<br/>
However, analyses of <b>ChIP-seq peaks</b> return <b>very good</b> results. See the program <i>peak-motifs</i>.
<p/>";
print "<hr>";


&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->start_multipart_form(-action=>"oligo-analysis.cgi");

print &SequenceChoice();
print "<b><a href='help.oligo-analysis.html#sequence_type'>Sequence type</a></b>";
print "&nbsp;"x3,  $query->popup_menu(-name=>'sequence_type',
				      -Values=>["dna","protein","other"],
				      -default=>$default{sequence_type});

#### purge sequences
print "<br>", $query->checkbox(-name=>'purge',
		       -checked=>$default{purge},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-analysis.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";


################################################################
## Oligomer counting options
print "<hr>\n";
print "<b>Oligomer counting mode</b><br>\n";

# ### oligo size
# print "<B><A HREF='help.oligo-analysis.html#oligo_length'>Oligomer length</A>&nbsp;</B>\n";
# print $query->popup_menu(-name=>'oligo_length',
# 			 -Values=>[1,2,3,4,5,6,7,8],
# 			 -default=>$default{oligo_length});



## Oligo sizes
print "<p><b><a href='help.oligo-analysis.html#oligo_length'>Oligomer lengths</a>&nbsp;</b>\n";
@oligo_lengths = 1..9;
for my $len (@oligo_lengths) {
    print "&nbsp;"x2;
    print $query->checkbox(-name=>"oligo_length".$len,
			   -checked=>$default{"oligo_length".$len},
			   -label=>$len);
}
print "<br><i>Note: motifs can be larger than oligo sizes (oligos are used as seed for building matrices)</i>";

### prevent overlapping matches of the same pattern
print "\n<br>", $query->checkbox(-name=>'noov',
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


print "<hr>\n";

################################################################
## Background model
&PrintOligoBackgroundOptions();

################################################################
#### pseudo frequencies for the BG model
print ("<b><a href=help.oligo-analysis.html#pseudo>Pseudo-frequency</a></b> &nbsp;",
       $query->textfield(-name=>'pseudo_freq',
			 -default=>$default{pseudo_freq},
			 -size=>5));

print "<hr>\n";

&OligoReturnTable();

print "<hr width=550 align=left>\n";

### send results by email or display on the browser
&SelectOutput();


### action buttons
print "<ul><ul><table class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"oligo-analysis_form.cgi");
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
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'bg_method', -default=>'background');
print $query->hidden(-name=>'background',-default=>"upstream-noorf");
print $query->hidden(-name=>'bg_level',-default=>"organism");
print $query->hidden(-name=>'gibbs_msps',-default=>"2");
print $query->hidden(-name=>'gibbs_flanks',-default=>"2");
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'title',-default=>'upstream sequences from the yeast PHO genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.oligo-analysis.html'>DEMO</A></B></TD>\n";
print "<td><b><a href='help.oligo-analysis.html'>MANUAL</A></B></TD>\n";
print "<td><b><a href='tutorials/tut_oligo-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<td><b><a href='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</a></b></td>\n";
print "</tr></table>
</ul></ul>\n";

print "</font>\n";
print "<hr>";

print $query->end_html;

exit(0);


################################################################
## Table with all the supported statistics and thresholds
sub OligoReturnTable {

print "<p><b>Result</b></p>\n";

print ("<INPUT TYPE='radio' NAME='return' VALUE='fields' checked>", "One row per pattern");

print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			 [
			  $query->th([" <A HREF='help.oligo-analysis.html#return_fields'>Return fields</A> ",
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
				     ]),

			  ### binomial E-value
			  $query->td(["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E-value",
				      $query->textfield(-name=>'lth_occ_E',
							-default=>$default{lth_occ_E},
							-size=>5),
				      $query->textfield(-name=>'uth_occ_E',
							-default=>$default{uth_occ_E},
							-size=>5),
				     ]),

			  ### significance index
			  $query->td(["&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Significance",
				   $query->textfield(-name=>'lth_occ_sig',
						     -default=>$default{lth_occ_sig},
						     -size=>5),
				   $query->textfield(-name=>'uth_occ_sig',
						     -default=>$default{uth_occ_sig},
						     -size=>5)
				   ]),

			  ### binomial Family-wise error rate (FWER)
			  $query->td([#$query->checkbox(-name=>'FWER',
#						       -checked=>$default{occ_FWER},
#						       -label=>' FWER '),
				      "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FWER",
				      $query->textfield(-name=>'lth_occ_FWER',
							-default=>$default{lth_occ_FWER},
							-size=>5),
				      $query->textfield(-name=>'uth_occ_FWER',
							-default=>$default{uth_occ_FWER},
							-size=>5),
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

#### Convert patterns to matrix
&MatrixFromPatterns_print_form();

print "</BLOCKQUOTE>\n";
print ("<INPUT TYPE='radio' NAME='return' VALUE='table'>", 
       "<b>Occurrence table</b>: one row per sequence, one column  per oligo (occurrence counts only, email output recommended)", "<P>\n");

print ("<INPUT TYPE='radio' NAME='return' VALUE='distrib'>", 
       "Pattern count distribubtions, one row per pattern (occurrence counts only, email output recommended)", "<P>\n");

}
