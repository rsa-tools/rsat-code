#!/usr/bin/perl
#### this cgi script fills the HTML form for the program matrix-scan
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

$default{demo_descr1} = "";
$default{demo_descr2} = "";
$default{demo_descr3} = "";

$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc

$default{origin}="end";
$default{offset}="0";

$default{bg_method}="bginput";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = "0";
$default{organism} = "Saccharomyces cerevisiae";
$default{matrix_format} = "tab";
$default{pseudo_counts} = 1;
$default{consensus_as_name} = "";
$default{pseudo_distribution} = "pseudo_prior";
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{decimals} = "1";
$default{n_score} = "score";
$default{crer_ids} = "";

## Return fields
## matches
$default{return_sites} = "CHECKED";
$default{return_pval} = "CHECKED";
$default{return_site_limits} = "CHECKED";
$default{return_rank} = "CHECKED";
$default{return_normw} = "";
$default{return_bg_residues} = "";
$default{return_matrix} = "";
$default{return_freq_matrix} = "";
$default{return_weight_matrix} = "";
$default{return_bg_model} = "";

$default{return_distrib} = "CHECKED";
$default{return_occ_proba} = "CHECKED";
$default{sort_distrib} ="occ_sig";

$default{return_crer} = "CHECKED";
$default{return_crer_limits} = "CHECKED";
$default{return_crer_sites} = "CHECKED";

$default{analysis_type} = "analysis_sites";
$checked{$default{analysis_type}} = "CHECKED";


## Threshold values for site detection
$default{lth_score} = "1";
$default{uth_score} = "none";
$default{lth_rank} = "none";
$default{uth_rank} = "none";
$default{lth_proba_M} = "none";
$default{uth_proba_M} = "none";
$default{lth_proba_B} = "none";
$default{uth_proba_B} = "none";
$default{lth_normw} = "none";
$default{uth_normw} = "none";
$default{lth_sig} = "none";
$default{uth_sig} = "none";
$default{lth_pval} = "none";
$default{uth_pval} = "1e-4";

## Threshold values for CRER detection
$default{lth_site_pval} = "none";
$default{uth_site_pval} = "1e-3";
$default{lth_crer_size} = "30";
$default{uth_crer_size} = "500";
$default{lth_crer_sites} = "none";
$default{uth_crer_sites} = "none";
$default{lth_crer_pval} = "none";
$default{uth_crer_pval} = "none";
$default{lth_crer_sig} = "2";
$default{uth_crer_sig} = "none";

## Threshold values for occurrence statistics
$default{lth_occ_score} = "0";
$default{uth_occ_score} = "none";
$default{lth_inv_cum} = "none";
$default{uth_inv_cum} = "none";
$default{lth_exp_occ} = "none";
$default{uth_exp_occ} = "none";
$default{lth_occ_pval} = "none";
$default{uth_occ_pval} = "none";
$default{lth_occ_eval} = "none";
$default{uth_occ_eval} = "none";
$default{lth_occ_sig} = "0";
$default{uth_occ_sig} = "none";
$default{lth_occ_sig_rank} = "none";
$default{uth_occ_sig_rank} = "3";




################################################################
#### STILL TO BE TREATED
# [-R <Set the range for approximating a weight matrix with integers (default: 10000)>]
# [-e <Small difference for considering 2 scores equal (default: 0.000001)>]
# [-li <Determine lower-threshold score from adjusted information content>]
# [-lp <Determine lower-threshold score from a maximum ln(p-value)>]


### print the form ###
&RSA_header("matrix-scan");
&ListParameters() if ($ENV{rsat_echo} >=2);

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
  if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
  if ($key eq "bg_method"){
  	$checked{$query->param($key)} = "CHECKED";
  }
  if ($key eq "analysis_type"){
  	$checked{$query->param($key)} = "CHECKED";
  }
}

&ReadMatrixFromFile();

### head
print "<center>";
print "Scan a DNA sequence with a profile matrix<br>\n";
print "</p>";
print "</CENTER>";
print "<b>Citation</b>: <a href='mailto:jturatsi\@bigre.ulb.ac.be (Jean Valery Turatsinze)'>Jean Val&eacute;ry Turatsinze</A>, <A HREF='mailto:morgane\@bigre.ulb.ac.be (Morgane Thomas-Chollier)'>Morgane Thomas-Chollier</A>, <a href='mailto:defrance@bigre.ulb.ac.be'>Matthieu Defrance</a> and <A HREF='mailto:jvanheld\@bigre.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</a> (2008). Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nat Protoc, 3, 1578-1588. <a href='http://www.ncbi.nlm.nih.gov/pubmed/18802439'>Pubmed 18802439</a>";

## demo description
print $default{demo_descr1};
print $default{demo_descr2};
print $default{demo_descr3};

print $query->start_multipart_form(-action=>"matrix-scan.cgi");

################################################################
#### sequence
print "<hr>";
&DisplaySequenceChoice();


################################################################
#### Matrix specification
print "<hr>";
&GetMatrix("consensus"=>1);

################################################################
## Background model
print "<hr>";

my %bg_params =("markov" => 1,
		"bg_input" => 1,
		"bg_window" => 1,
		"markov_message" => 1
	       );
&GetBackgroundModel(%bg_params);

print "<hr>";



################################################################
#### strands

print "<p><B>Scanning options</B><br>\n";
print "<BR>\n";
print "<A HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['single',
				   'both'],
			 -default=>$default{strands});

################################################################
#### origin for calculating positions
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#origin'><B>Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['start',
				   'center',
				   'end',
				   'genomic'],
			 -default=>$default{origin});

################################################################
#### Offset for calculating positions
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#offset'><B>Offset</B></A>\n";
print $query->textfield(-name=>'offset',
			-default=>$default{offset},
			-size=>8);

################################################################
#### decimals
print "&nbsp;"x2,  "<A HREF='help.matrix-scan.html'><B>score decimals</B></A>\n";
print $query->popup_menu(-name=>'decimals',
			 -Values=>['0',
				   '1','2'],
			 -default=>$default{decimals});

################################################################
#### decimals
print "&nbsp;"x2,  "<A HREF='help.matrix-scan.html'><B>handling of N characters</B></A>\n";
print $query->popup_menu(-name=>'n_score',
			 -Values=>['score',
				   'skip'],
			 -default=>$default{n_score});


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


## Load demo sequences
$demo_sequence_file = $ENV{RSAT}."/public_html/demo_files/matrix_scan_demo_sequences.fasta";
$demo_sequence = `cat $demo_sequence_file`;

## Load demo matrices
$demo_matrix_file = $ENV{RSAT}."/public_html/demo_files/matrix_scan_demo_matrices.tf";
$demo_matrix = `cat $demo_matrix_file`;



$descr = "<H4>Comment on the demonstration example : </H4>";

$descr .= "<blockquote class ='demo'>";

$descr .= "In this demonstration, we will analyse the promoter of
<i>Drosophila melanogaster</i> even-skipped gene (<i>eve</i>). We will scan the 5500
bp sequence upstream the transcription start site with matrices
representing the binding specificity of 12 transcription factors known
to regulate <i>eve</i>. These matrices were built from binding sites
annotated in the <a target=_blank
href='http://www.oreganno.org'>ORegAnno</a> database by Jean-Valery
Turatsinze.<p/>";

## demo 1
print "<TD><B>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr."The program will return individual matches, i.e. sequence segments scoring above the predefined threshold. In this example, threshold is set on the P-value.
</blockquote>");
print $query->hidden(-name=>'bg_method',-default=>'bgfile');
print $query->hidden(-name=>'uth_pval',-default=>'1e-4');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'markov_order',-default=>'0');
print $query->hidden(-name=>'organism',-default=>'Drosophila_melanogaster');
print $query->hidden(-name=>'analysis_type',-default=>'analysis_sites');
print $query->hidden(-name=>'return_rank',-default=>'');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'consensus_as_name',-default=>'');
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'origin',-default=>'genomic');
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->submit(-label=>"DEMO 1 (sites)");
print "</B></TD>";
print $query->end_form;

## demo2
print $query->start_multipart_form(-action=>"matrix-scan_form.cgi");
print "<TD><B>";

print $query->hidden(-name=>'demo_descr2',-default=>$descr."The program will return CRERs: regions of a few hundreds residues that have a higher density of matches than expected by chance.
</blockquote>");
print $query->hidden(-name=>'bg_method',-default=>'bgfile');
print $query->hidden(-name=>'uth_site_pval',-default=>'1e-4');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'markov_order',-default=>'0');
print $query->hidden(-name=>'organism',-default=>'Drosophila_melanogaster');
print $query->hidden(-name=>'analysis_type',-default=>'analysis_crer');
print $query->hidden(-name=>'return_rank',-default=>'');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'consensus_as_name',-default=>'');
print $query->hidden(-name=>'crer_ids',-default=>'');
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'origin',-default=>'genomic');
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->submit(-label=>"DEMO 2 (CRERs)");
print "</B></TD>\n";
print $query->end_form;

## demo3: detect enrichment of hits for PSSMs
print $query->start_multipart_form(-action=>"matrix-scan_form.cgi");
print "<TD><B>";
print $query->hidden(-name=>'demo_descr3',-default=>$descr."The program will return matrices for which the total number of hits in the input sequences is higher than expected by chance.</blockquote>");
print $query->hidden(-name=>'bg_method',-default=>'bgfile');
print $query->hidden(-name=>'uth_site_pval',-default=>'1e-3');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'markov_order',-default=>'0');
print $query->hidden(-name=>'organism',-default=>'Drosophila_melanogaster');
print $query->hidden(-name=>'analysis_type',-default=>'analysis_occ');
print $query->hidden(-name=>'return_rank',-default=>'');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'consensus_as_name',-default=>'');
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'origin',-default=>'genomic');
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->hidden(-name=>'uth_occ_sig_rank',-default=>1);
print $query->hidden(-name=>'lth_occ_score',-default=>5);
print $query->submit(-label=>"DEMO 3 (enrichment)");
print "</B></TD>\n";
print $query->end_form;

#$demo_sequence = ">MET8	YBR213W; upstream from -463 to -1; size: 463; location: NC_001134.7 649900 650362 D; upstream neighbour: YBR212W (distance: 463)
#TTACAAAAGACAAAAAAAGAAAATTTTAATCTTGTCCGCAGTTTTATCTGCGTCTCTACG
#TTCTTACGTTTCTTCTATTAATGCCATTTCAGTTACAACCTAGTCAATTGTCGATCCATA
#ATTCTAATCAAATTTGTTTTTCCTCTATACTACCTATCTATTTTTATCTATCTAAGTACA
#TTTATTTACTCAAACAGTTCCGTTTCAAAGTGTTTTATATTAACTATATATGCGAAAAGC
#TGGCGTCATAATTTCACGTGTTATAATAGCCATGCTGACGGAAAAAAAATGTGAAAATCG
#CTACAAAGTCCGATGACTACGGGCAGTAGCATGTAAATGATGGACACACACACACACATA
#TATATATATATACATTTACTTCAATAAAAGGCTGTGCCAGACATTTTTGCCATACATTGT
#TCATGAAGTGTGCAAAATAAGAGAGTGTATAATAGGATAAAAA
#>MET32	YDR253C; upstream from -547 to -1; size: 547; location: NC_001136.8 964562 965108 R; upstream neighbour: YDR254W (distance: 547)
#TAATTGCTACTCAAATATACTAGTCAAAGATAGTATCCACCAAAATCTTTCCCCGCTAAA
#ATAACGCCAGATGCTTTCTATGCTTCTAATCTTTTACCATTTACCTTTGTTTATTTCAAT
#ATAAACTTTAATTTACAGTCCCTATCTATTGCCCGACTGGACTAACATGCACGTGACATT
#TTGTGATGGTTTTTCGTCCCTTACTTAGTACGCTTAGTACGCCACAGTTTATATTTTCTT
#GACAATAATAAAGAACCTGATTGTGGGTTAGAACTTGCTATACTTTTAGTTTAAAATAAG
#CAGGAAATAATCTTGAGTTCTGTATCATTATTATAAATAAAACTATATTTGTTCTCTTTG
#TCGCCCTCGGAACTTTCCTCATTACATTGACGAGGTATATATAGATATAGTAGATATACA
#TATCTATCCATGGTATATATGTATGCATCTGGATAATTGAATAGGGTTTCATGTCATATG
#CCAAGAATTTGTTAATAATATAGTGGAAAAAAGTCAAGAGGTATTATAAATTTCAAAAAA
#GTACCAA
#>MET18	YIL128W; upstream from -568 to -1; size: 568; location: NC_001141.1 113238 113805 D; upstream neighbour: YIL129C (distance: 568)
#TTTGATGTATAACAAAACTAAAAAGGGTTATTAAAATGGGAACACAACAAACAACCAGAA
#TTTTCACACTTTAACCAGTCACGTCCTATTATGAAGACCTAAATCCACATTTGCTTTCTC
#TCTTCATTTGCCTAATCCTTTATCCCAATTTCTACAGTTCTATATGTATTTTCCTGTGTG
#GCTGTCGTTTCGTGGTTAGTGATACAACCATAACGATTCAACCAACTCCCAATGTATGTG
#ATGTTGATACCGCTAATTTGGAAGGGATGGTATACTCTAGGTGACCTCAATGAGTCAAAG
#AGAGCTAGGACATACTTCGAGATAGGTAATACCACTTTGCAGCTTCTTTTTAGGCCTTCA
#TGAGTGAGTAGCCAAGAAAAAGTTAAAAGCGGGTAATAGGTATGAATTTTTCAAATACTG
#AAATTTGGTTTAGTTATTTAAGTGAATTGTAGATTATGTACATTTTACGTGCAATGAAGG
#AGTCACCTCTATGATCATCTAGTTATTAGCTGTTAGTTTTCATTGAACTTGTTTTAACTG
#GGAAAAAGCGGAACAATTGGGCCTTACA
#>MET30	YIL046W; upstream from -177 to -1; size: 177; location: NC_001141.1 268473 268649 D; upstream neighbour: YIL046W-A (distance: 177)
#CACGTGATCGGGAAGCCACAGTTTGCGCGGAGATATTTTATTTTTTTTCATCAGCGTAAG
#AAGAAAGCAACCTTGCAGTCTGTATCGTAAGAGAAGACTGCAGTTAAAGAAGTTTAGAGA
#AGAGGCTTGAGTATCGGTAAAGGGGTGTGTGTTTGGTGATTTATAAAGGAGAAGGGC
#>MET28	YIR017C; upstream from -489 to -1; size: 489; location: NC_001141.1 384117 384605 R; upstream neighbour: YIR018W (distance: 489)
#GACTGTGATAATATGCTAGTTACACTGTTTATGTTGTGTGAACTTGTTGTAATATGGTTA
#ACTTCACTTTCAGTGATTGATATGATAGCGACATCACTGCCGTGCAAAAAGACCATTCCA
#TTACTGCACCTTTTTGTCCTTTTCCGTGGAATAAAAGTTCACTCGTCAGTTCCATGCATT
#CTGGAAAAAAATGATCTGAAAGATGCCACAGTTGTGGGGCCCGCCCGGCCCAATAGGTAA
#ACTAAAATACAATAGAAGGGGTACTGAGTGCACGTGACTTATTTTTTTTTTTTGGTTTTA
#GGTTTCGCTTTTTTCACCTTTTTCTACTTTCTAACACCACAGTTTTGGGCGGGAAGCGGA
#AACGCCATAGTTGTAGGTCACTGGCGTGAGTCAAGGCCGGGCAGCCAATGACTAAGAACA
#CGAGGTAACTTGAATTTAACTATTTATAACCAGTGGTAGTTACGAAGACAAATTGTTTTG
#TTCGTCAAT
#>MET6	YER091C; upstream from -687 to -1; size: 687; location: NC_001137.2 342164 342850 R; upstream neighbour: YER092W (distance: 687)
#TTTTTTCTTGTTTTATAATCAGTCAAGTATTGGTTTCCCACAGCCATTCAACTCAGGTTC
#ATCATCTTTTTCGCTTCCAAAAATGCAGTTGATTTCACACAATTTTTCATGAACCAGGGT
#CCCGCACTCCGGGTAAAGGACCATCACGCCACATCACGTGCACATTACTAGTAAAAGCCA
#CAGGAAATATTTCACGTGACTTACAAACAGAGTCGTACGTCAGGACCGGAGTCAGGTGAA
#AAAATGTGGGCCGGTAAAGGGAAAAAACCAGAAACGGGACTACTATCGAACTCGTTTAGT
#CGCGAACGTGCAAAAGGCCAATATTTTTCGCTAGAGTCATCGCAGTCATGGCAGCTCTTT
#CGCTCTATCTCCCGGTCGCAAAACTGTGGTAGTCATAGCTCGTTCTGCTCAATTGAGAAC
#TGTGAATGTGAATATGGAACAAATGCGATAGATGCACTAATTTAAGGGAAGCTAGCTAGT
#TTTCCCAACTGCGAAAGAAAAAAAGGAAAGAAAAAAAAATTCTATATAAGTGATAGATAT
#TTCCATCTTTACTAGCATTAGTTTCTCTTTTACGTATTCAATATTTTTGTTAAACTCTTC
#CTTTATCATAAAAAAGCAAGCATCTAAGAGCATTGACAACACTCTAAGAAACAAAATACC
#AATATAATTTCAAAGTACATATCAAAA
#>MET10	YFR030W; upstream from -338 to -1; size: 338; location: NC_001138.4 212962 213299 D; upstream neighbour: YFR029W (distance: 338)
#TGCATCTAAATATATACGTATGTTTAAGGTTCTGGTATACAGGTATTAAAAGAAAACACT
#ATCAACATTCCCAATAAGATATACCACACCACGTGAGCTTATAGAAGCACGTGACCACAA
#TTCACCCCACAGGTGTGGCTTTTTTGGTGCCGTAGAAAAGACTCATTCATGAATCGTCGG
#AAACCCATAGTCATCTTCGAGCAAAAGGTATATATAAGCAACAGAGGGCAGTAGTTCTCG
#AGACCACCATCTTTTGATTGGAAATAGTTTCGTTTAGATGGGGTGCACATAGTTTTTTTC
#AACTGCTTTTCCTCGAGGTCACCCAAATATACAACGAG
#>MET13	YGL125W; upstream from -380 to -1; size: 380; location: NC_001139.7 272146 272525 D; upstream neighbour: YGL126W (distance: 380)
#CTCAGGAAAAGTTGGCGATAGACCACGAGCGACTGAAAAAATAACAGCGACTTTTCTCCC
#GGTAGCGGGCCGTCGTTTAGTCATTCTATCCCTCGGATTATAGACTGTGAATATTGCATA
#TGCAACTTTGACTCAAATTTTTCCAAAATTTGATATATATATATATATATATATGTTTGT
#ATGTATATATATATATACGTATATATATCATATATACGAAAAGTAGAAAAAAAAAGGTGA
#TATTTCGCTCGTGGAAAAGCTAATGCCACAGCTTGTGTTTCGTGTAGTTTGCCTTGCTCC
#CCTTGATTGAAATAGTCTCCCTAAACTAAAGTTATCAGCAAACAGAACCACCACAGTTAC
#TACTACAACCACATCGCAAT
#>MET3	YJR010W; upstream from -800 to -1; size: 800; location: NC_001142.6 455354 456153 D; upstream neighbour: YJR009C (distance: 1557)
#AAGAGTACAATTTATAAATTAATGAAAACACAGAAGTATTTAGATCGGCTCAAATGTTTT
#TGGACATTAAAAGATCTTGAAACTGAGTAAGATGCTCAGAATACCCGTCAAGATAAGAGT
#ATAATGTAGAGTAATATACCAAGTATTCAGCATATTCTCCTCTTCTTTTGTATAAATCAC
#GGAAGGGATGATTTATAAGAAAAATGAATACTATTACACTTCATTTACCACCCTCTGATC
#TAGATTTTCCAACGATATGTACGTAGTGGTATAAGGTGAGGGGGTCCACAGATATAACAT
#CGTTTAATTTAGTACTAACAGAGACTTTTGTCACAACTACATATAAGTGTACAAATATAG
#TACAGATATGACACACTTGTAGCGCCAACGCGCATCCTACGGATTGCTGACAGAAAAAAA
#GGTCACGTGACCAGAAAAGTCACGTGTAATTTTGTAACTCACCGCATTCTAGCGGTCCCT
#GTCGTGCACACTGCACTCAACACCATAAACCTTAGCAACCTCCAAAGGAAATCACCGTAT
#AACAAAGCCACAGTTTTACAACTTAGTCTCTTATGAAGTTACTTACCAATGAGAAATAGA
#GGCTCTTTCTCGAGAAATATGAATATGGATATATATATATATATATATATATATATATAT
#ATATATGTAAACTTGGTTCTTTTTTAGCTTGTGATCTCTAGCTTGGGTCTCTCTCTGTCG
#TAACAGTTGTGATATCGTTTCTTAACAATTGAAAAGGAACTAAGAAAGTATAATAATAAC
#AAGAATAAAGTATAATTAAC
#>MET14	YKL001C; upstream from -800 to -1; size: 800; location: NC_001143.7 439029 439828 R; upstream neighbour: YKR001C (distance: 1222)
#TATTTTTTTAATTACATAATCATAAAAATAAATGTTCATGATTTCCGAACGTATAAAATA
#AGAATGTTACGAGAATTTGTTTTCTTGGTAATTAAAATAATCAAATACACATAGAAAGGA
#GAGTAAACTGCTTCCTCTGTATAAATCAAAGCAAAATTGTAAATAGCGTTGACAAGTGAT
#TACAGAAGTTAGGTGAGGTTAATTACCAATTTCTTTTTTTAAAATTGGTGAAATAAGATT
#ACGTTTAAAGGAGCATTAACAGGTTTACTCATAACAATCATTTTCAAATTTCCCTATGCA
#TGTTTAGAGCAAGCGCCTTTGTGAGCCCTCCCGGTTACGACGCCTTGGCAATGTAGCAGA
#TAACTCTGCACTTCTAGAATCATTCCACTACGACATTTGGCTCATCACCAGCTCGCGAGA
#AATGTAAATAAGCCAACAACCAAGAATGCGTAACATTAAAGAATACAGTTGCTTTCATTT
#CGGCGTGATGGTACGGCACCCACGGTACCTTACATTATTCTCGAAAAATAGCTGCACGCT
#TTTCCAGGAATAAAAGACCGTGCCACTAATTTCACGTGATCAATATATTTACAAGCCACC
#TCAAAAAATGTGGCAATGGAGAAGAGGATGAACGACTCAATATGACTTCAACTTCATGAA
#TTTGTCAAAATATCTATATAAGATGCAAAATTTCTATACAACATCAGTTGCGTATCCGTT
#AATGTCGTTCATTTTCTCTCTTTGTTCGAACTTGACATCAAGAAAAGTTGGAATTATTTC
#TCCAAGCACACTGTACACCA
#>MET1	YKR069W; upstream from -702 to -1; size: 702; location: NC_001143.7 570552 571253 D; upstream neighbour: YKR068C (distance: 702)
#TTTTGACCCAGTTTTGGTTATCAATGAACACTTGAAGCTTTACTCTGCATTCCCATCTCT
#ATAGCTATGGGTAATCACAGCTACGATCACTTACTCTGTTATTATTATATTAAGTTCAAT
#GTTGGCCAAACCGGGTAACATGTAACACTTTCAGGTTGGCCTTACCTTTGGCTTGGAGTT
#TCGCAAGTTTTCAAATTTTTGGCTCCTGCTGTCAAGGTGCATAGAATAGCGCTTATTTAT
#CTATTTATATCCAAGATGTACAATCCTCGTTCTCTGAGTCCAACATATTTGCTCGCAACT
#GTAGAAATCACAACTACAGCAACAGTAAAGATATCATTTTCTATTTTCGTTATTGGTTTC
#TCGACCTTTTTATATACGATACGTCAAACTTGAATCATTTTATACGTTTTTCTCTTTCTA
#GAAATGCCATTATGCACGTGACATTACAAATTGTGGTGAAAAAAGGCTCTCATAATAAAC
#TGTGAACGGACTCATAATGAAATTTGCTTCACTATGTGAATCATCGCTAATAAACTCGCT
#ACAAAAGTCGAGTATGCTTAAGTCAAAAAAATGATATATATATATAATTTACTTATGTGT
#TTCTGCAAAGTTGTAGGCTTCATTTAGAATTGCTCAGATATTCCATCCCAATTAAAAAAA
#GCACGGATAGAGTGATAAATAAACTAAGAAAATTTCAAAAGA
#>MET17	YLR303W; upstream from -800 to -1; size: 800; location: NC_001144.4 731744 732543 D; upstream neighbour: YLR301W (distance: 982)
#TATACTAGAAGTTCTCCTCGAGGATTTAGGAATCCATAAAAGGGAATCTGCAATTCTACA
#CAATTCTATAAATATTATTATCATCGTTTTATATGTTAATATTCATTGATCCTATTACAT
#TATCAATCCTTGCGTTTCAGCTTCCACTAATTTAGATGACTATTTCTCATCATTTGCGTC
#ATCTTCTAACACCGTATATGATAATATACTAGTAACGTAAATACTAGTTAGTAGATGATA
#GTTGATTTTTATTCCAACACTAAGAAATAATTTCGCCATTTCTTGAATGTATTTAAAGAT
#ATTTAATGCTATAATAGACATTTAAATCCAATTCTTCCAACATACAATGGGAGTTTGGCC
#GAGTGGTTTAAGGCGTCAGATTTAGGTGGATTTAACCTCTAAAATCTCTGATATCTTCGG
#ATGCAAGGGTTCGAATCCCTTAGCTCTCATTATTTTTTGCTTTTTCTCTTGAGGTCACAT
#GATCGCAAAATGGCAAATGGCACGTGAAGCTGTCGATATTGGGGAACTGTGGTGGTTGGC
#AAATGACTAATTAAGTTAGTCAAGGCGCCATCCTCATGAAAACTGTGTAACATAATAACC
#GAAGTGTCGAAAAGGTGGCACCTTGTCCAATTGAACACGCTCGATGAAAAAAATAAGATA
#TATATAAGGTTAAGTAAAGCGTCTGTTAGAAAGGAAGTTTTTCCTTTTTCTTGCTCTCTT
#GTCTTTTCATCTACTATTTCCTTCGTGTAATACAGGGTCGTCAGATACATAGATACAATT
#CTATTACCCCCATCCATACA
#>MET2	YNL277W; upstream from -481 to -1; size: 481; location: NC_001146.5 116868 117348 D; upstream neighbour: YNL277W-A (distance: 481)
#GCAGTATAAATTGTACTTCAAAGCACTAGTCATGAAAAACGCTTACATTAGTTCAGTTTG
#TCAAGGTTATGCTATTACTTGTACTTATTTCTTGCTATTGTTAGTGGCTCCCCACATTGA
#CGTATTTTCACGTGATGCGCCTCACTGCGGAAGGCGCCACACATTGCCTGCAAAAAATTG
#TGGATGCACTCATTTGATAGTAAACTAAGTCATGTTAATCGTTTGGATTTGGCACACACC
#CACAAATATACACATTACATATATATATATATTCAAAATACAGCTGCGTCCAATAGATGA
#GCTTCCGCTTCGTTGTACAACCTACCTGCTATCTTGTTCACGGATATTTCTTGCTTTTAA
#TAAACAAAAGTAACTCTAGAACAGTCAAGTCTTCGATAATTTTTTTAGTCACAGGGTCCG
#TCTAAAGTTTCTCTTTATTTGGAATAATAGAAAAGAAAGAAAAAAACGTAGTATAAAAGG
#A
#>MET4	YNL103W; upstream from -800 to -1; size: 800; location: NC_001146.5 426937 427736 D; upstream neighbour: YNL104C (distance: 980)
#AGAGGCTGACCCAAGAGGAGAAAACATCGAACCTACGCGACTCCATAGCGCACATCTCCC
#ATGCGCCCGTGCTTATATATATATAAATATACATACACACATACATGCACGCACATACAT
#GCACGCGCATGACAGTGACTGGCCGTTACTGTACAATTTTTTCAGCCAAGTATGACACAC
#ATTCAACTCAGCTTTTCTGAGGCCTTCTTTCTTTTCCTGCGCGTCGGTAGAGCGATGACT
#AACCTACTACTGTCTCAGAGCCGGTCCCGCTCCGGTAGCAATCCTGGGGCTGGTCATAAC
#AGCCGAGTGGAAGTGTCAAAGCGGAGAACAGAAGCATAAGCTCAATCGCTGGACATACGG
#ATGCTTATATACGTCTTATTGTCGTTGAAAAATATCGAATTTTTTACTTCATTTATCGAG
#GCTTCTTCGAGCACTTTTCCGCTATGGCTTTTTCCCCGTTTCCTTTTAATCACGTGCGCG
#GGTAGCACCCGGCACACAGCTGGTGTCTCGTCGCACATGCTATTGTGTGTCATCGGGCCA
#CACAAGCATATTGCTTGAATTTTCTTTCATCGTTCAACTTAAATCCACCCAATCTAGATG
#TAGCCGTAGCATGTAATAACGTATATCCTTGTTTACATGCATCTGTGCCAGGTGAAACGG
#TCTGTTTGAACGCACATCATTTCATATATTAGTCAACTCCTGAAGGTCTTCTTGCCTGTC
#CGTCAACTGTTTAGACAGACTCTCGTCAATAAAGCGCACTTCTGATAAGCACTTTTATTC
#CTTTTTTTCCACTGTGAACG
#>MET22	YOL064C; upstream from -215 to -1; size: 215; location: NC_001147.5 207177 207391 R; upstream neighbour: YOL063C (distance: 215)
#CATATTTTGACATTTACATAGCCATCTATATAATAATCCTTCCTTCATTGAAATGCGCGA
#ATGACTCAGACGAGCAATATCACGTGTTGCGATTTTACTTTCAGTTTGCAAAGAAGAAAC
#CATGCCACATATAAAGAACGTTTGCACTTCTCTTTAATTTATTAGTTAGTAAGTAAGAAG
#TTTAAAGACAACTCAGAAGACATCAGCACTTTACT
#>MET7	YOR241W; upstream from -250 to -1; size: 250; location: NC_001147.5 786746 786995 D; upstream neighbour: YOR239W (distance: 250)
#GAAGTTCTGAGACAAGTACCACCTCCTCTCTCATCATAAAACAAGTAAAAGTTTTCTCGT
#CGCGCATATTATTTTGGTGATTGATTGTTTTTTCCTCCGATATCATCACTTATTACCTGT
#AATTTTATCTTTTTCTACCCCATAGAATTCGTCTTATAAGTCTATACCCTCAAAACTATC
#TATCATTTTAATATTATCTGTCGCTTTAATTGTCTTATTTCTGAAGCTCACTGAAGAACA
#TTGCTTTATT
#>MET31	YPL038W; upstream from -161 to -1; size: 161; location: NC_001148.3 480371 480531 D; upstream neighbour: YPL038W-A (distance: 161)
#ATTTTAGTCTAAAAATTTTGCTAGCCCATCAATTTTTTTTTTGTTCTAATGCAAAATATA
#ACATGGGTAAGAAAAAGAAAAAGCCGTTCCTCAGTACGTAAAGAGATTTGATCATTAACA
#AGTTGGGCTCAATATACACAGTCGATAGTCTATATGTGCAT
#>MET12	YPL023C; upstream from -384 to -1; size: 384; location: NC_001148.3 506311 506694 R; upstream neighbour: YPL022W (distance: 384)
#CTGGAAAGATATTTTCAACAGGATAGTGCAATATTATTTTTACACATTTAGCAAATGCTC
#TACCAACGTCCTGAACCCTCCAGTACACCTGTGCTCTCTCCTCTCCTATGCGCCACGCAG
#ACAAAAGTTTACTCGTCCCGACTTTTTTTTTTCATTAGACGCGATATTGACTGTGGCTAT
#AGCTTACTCCAGGAGATAAGCGAGTAAAGCTTTTTCTAACTTCAATGATGAAGAAAAGTC
#GCAAAATAAAGGCAAACAGAGAACACTTCAGGTTGTTGGTTACATTGGAAGAGGGACTTA
#AGCTCTCACATCATCTATTTTGTTTCAAGTTCGTACATTTTTTGAAGCGTGTTGGACGGG
#ACAGGTTGATTACATTTTTTAAAC
#>MET16	YPR167C; upstream from -443 to -1; size: 443; location: NC_001148.3 877629 878071 R; upstream neighbour: YPR168W (distance: 443)
#CTTATCGGTTTATTTTTCTATATATTTGCCTCTTTCTCAAACAGGAGTTAGTAGTTAAAA
#GTACGAAGTTCTTGTTCTTTAATGCGCGCTGACAAAAGAATTGGATAAAAGAGAATGGTG
#GGGGGACAAGAAGGAAATTTGTCCTAGTTTAACATGAATGGCATCTTGTTACCGGGTGGA
#CATCACCTATTGATTCTAAATATCTTTACGGTTTATCATACTGTTCTTTATTCCGTCGTT
#ATTCTTTTTATTTTTATCATCATTTCACGTGGCTAGTAAAAGAAAAGCCACAACATGACT
#CAGCAAATCTCGACAAAGTAAAAGCTCATAGAGATAGTATTATATTGATATAAAAAAAGT
#ATACTGTACTGTTTGTAACCTTTTCAATGCTTTAAGATCAAAACTAAGGCCAGCAAAGGT
#ATCAACCCATAGCAACTCATAAA
#";
#
#$demo_matrix = "
#; MET4 matrix, from Gonze et al. (2005). Bioinformatics 21, 3490-500.
#A |   7   9   0   0  16   0   1   0   0  11   6   9   6   1   8
#C |   5   1   4  16   0  15   0   0   0   3   5   5   0   2   0
#G |   4   4   1   0   0   0  15   0  16   0   3   0   0   2   0
#T |   0   2  11   0   0   1   0  16   0   2   2   2  10  11   8
#//
#; MET31 matrix, from Gonze et al. (2005). Bioinformatics 21, 3490-500.
#A |   3   6   9   6  14  18  16  18   2   0   0   0   1   3   8
#C |   8   3   3   2   3   0   1   0  13   2   0   1   0   3   6
#G |   4   3   4   8   0   0   1   0   2   0  17   1  17  11   1
#T |   3   6   2   2   1   0   0   0   1  16   1  16   0   1   3
#";
#
#print "<TD><B>";
#print $query->hidden(-name=>'organism',-default=>'Saccharomyces_cerevisiae');
#print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
#print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
#print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
#print $query->hidden(-name=>'alphabet',-default=>"a:t 0.325 c:g 0.175");
#print $query->submit(-label=>"DEMO");
#print "</B></TD>\n";
#print $query->end_form;


print "<TD><B><A HREF='help.matrix-scan.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_matrix-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

&ListParameters() if ($ENV{rsat_echo} >= 2);
&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->end_html;

exit(0);


################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {
  print "<p><b>Return</b> (Select one return type) </p>\n";


  #############################################
  ## Return fields
  #
  my $boxes_matches = "";
  @return_fields_matches = qw(sites pval rank );
  foreach my $field (@return_fields_matches) {
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,
				       -CHECKED=>$default{'return_'.$field},
				       -label=>' '.$field.' ');
  }
  $boxes_matches .= "<BR/>";
  @return_fields_matches = qw( site_limits normw);
  foreach my $field (@return_fields_matches) {
  		my $display_field = $field;
  		$display_field =~ s/site_//;
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,
				       -CHECKED=>$default{'return_'.$field},
				       -label=>' '.$display_field.' ');
  }
  $boxes_matches .= "<BR/>";
  @return_fields_matches = qw(weight_limits bg_residues);
  foreach my $field (@return_fields_matches) {
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,
				       -CHECKED=>$default{'return_'.$field},
				       -label=>' '.$field.' ');
    $boxes_matches .= "<BR/>";
  }


  ### Return fields
  my $boxes_occ = "";
  @return_fields_occ = qw(distrib);
  foreach my $field (@return_fields_occ) {
    $boxes_occ .= $query->checkbox(-name=>'return_'.$field,
				   -CHECKED=>$default{'return_'.$field},
				   -label=>' '.$field.' ');
  }
  $boxes_occ .= "<BR/>";
  @return_fields_occ = qw(occ_proba);
  foreach my $field (@return_fields_occ) {
    $boxes_occ .= $query->checkbox(-name=>'return_'.$field,
				   -CHECKED=>$default{'return_'.$field},
				   -label=>' '.$field.' ');
  }

  $boxes_occ .= "&nbsp;&nbsp; <b> sort by </b> &nbsp;&nbsp;".$query->popup_menu(-name=>'sort_distrib',
										-Values=>['scores',
											  'occ_sig'],
										-default=>$default{sort_distrib});

  ### Return fields
  my $boxes_crer = "";
  @return_fields_crer = qw(crer normw);
  foreach my $field (@return_fields_crer) {
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,
				    -CHECKED=>$default{'return_'.$field},
				    -label=>' '.$field.' ');
  }
  $boxes_crer .= "<BR/>";
  @return_fields_crer = qw(crer_limits);
  foreach my $field (@return_fields_crer) {
  	  	my $display_field = $field;
  		$display_field =~ s/crer_//;
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,
				    -CHECKED=>$default{'return_'.$field},
				    -label=>' '.$display_field.' ');
  }
  @return_fields_crer = qw(crer_sites);
  foreach my $field (@return_fields_crer) {
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,
				    -CHECKED=>$default{'return_'.$field},
				    -label=>' '."sites".' ');
  }
  $boxes_crer .= "<BR/>";
  $boxes_crer .= $query->checkbox(-name=>'crer_ids',
				  -CHECKED=>$default{crer_ids},
				  -label=>' '."crer-specific identifiers".' ');

  ### Return fields
  my $boxes_add = "";
  @return_fields_add = qw(matrix freq_matrix weight_matrix bg_model);
  foreach my $field (@return_fields_add) {
    $boxes_add.= $query->checkbox(-name=>'return_'.$field,
				  -CHECKED=>$default{'return_'.$field},
				  -label=>' '.$field.' ');
  }

  #############################################
  ## Thresholds
  #
  my $thresh_matches =
    $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
		  $query->Tr({-align=>center,-valign=>MIDDLE},
			     [
			      $query->th([" <A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A>"]),

			      ### Threshold on score
			      $query->td(['Weight<br>score',
					  $query->textfield(-name=>'lth_score',
							    -default=>$default{lth_score},
							    -size=>5),
					  $query->textfield(-name=>'uth_score',
							    -default=>$default{uth_score},
							    -size=>5)
					 ]),

			      ### Threshold on P-value of the score
			      $query->td(['P-value',
					  $query->textfield(-name=>'lth_pval',
							    -default=>$default{lth_pval},
							    -size=>5),
					  $query->textfield(-name=>'uth_pval',
							    -default=>$default{uth_pval},
							    -size=>5)
					 ]),
					  
				### Threshold on Sig of the score
			      $query->td(['Sig',
					  $query->textfield(-name=>'lth_sig',
							    -default=>$default{lth_sig},
							    -size=>5),
					  $query->textfield(-name=>'uth_sig',
							    -default=>$default{uth_sig},
							    -size=>5)
					 ]),

			      ### Threshold on proba_M
			      $query->td(['P(S|M)<br>proba_M',
					  $query->textfield(-name=>'lth_proba_M',
							    -default=>$default{lth_proba_M},
							    -size=>5),
					  $query->textfield(-name=>'uth_proba_M',
							    -default=>$default{uth_proba_M},
							    -size=>5)
					 ]),

			      ### Threshold on proba_B
			      $query->td(['P(S|B)<br>proba_B',
					  $query->textfield(-name=>'lth_proba_B',
							    -default=>$default{lth_proba_B},
							    -size=>5),
					  $query->textfield(-name=>'uth_proba_B',
							    -default=>$default{uth_proba_B},
							    -size=>5)
					 ]),

			      ### Threshold on normw
			      $query->td(['Normalized<br>weight',
					  $query->textfield(-name=>'lth_normw',
							    -default=>$default{lth_normw},
							    -size=>5),
					  $query->textfield(-name=>'uth_normw',
							    -default=>$default{uth_normw},
							    -size=>5)
					 ]),

			      ### Threshold on rank
			      $query->td(['Rank',
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

  ## Occurrences
  my $thresh_occ =
    $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
		  $query->Tr({-align=>center,-valign=>MIDDLE},
			     [
			      $query->th(["<A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A> "]),

			      $query->th({-colspan=>3,-align=>left},["Occurrences"
								    ]),
			      ### Threshold on score
			      $query->td(['Weight<br>score',
					  $query->textfield(-name=>'lth_occ_score',
							    -default=>$default{lth_occ_score},
							    -size=>5),
					  $query->textfield(-name=>'uth_occ_score',
							    -default=>$default{uth_occ_score},
							    -size=>5)
					 ]),

				### Threshold on occ_inv_cum
			      $query->td(['Occurrences<br>above the score',
					  $query->textfield(-name=>'lth_inv_cum',
							    -default=>$default{lth_inv_cum},
							    -size=>5),
					  $query->textfield(-name=>'uth_inv_cum',
							    -default=>$default{uth_inv_cum},
							    -size=>5)
					 ]),

			      $query->th({-colspan=>3,-align=>left},["Enrichment"
								    ]),

				### Threshold on exp_occ
			      $query->td(['Expected<br>occurrences',
					  $query->textfield(-name=>'lth_exp_occ',
							    -default=>$default{lth_exp_occ},
							    -size=>5),
					  $query->textfield(-name=>'uth_exp_occ',
							    -default=>$default{uth_exp_occ},
							    -size=>5)
					 ]),

			      ### Threshold on P-value of the score
			      $query->td(['Occ P-value',
					  $query->textfield(-name=>'lth_occ_pval',
							    -default=>$default{lth_occ_pval},
							    -size=>5),
					  $query->textfield(-name=>'uth_occ_pval',
							    -default=>$default{uth_occ_pval},
							    -size=>5)
					 ]),

			      ### Threshold on P-value of the score
			      $query->td(['Occ E-value',
					  $query->textfield(-name=>'lth_occ_eval',
							    -default=>$default{lth_occ_eval},
							    -size=>5),
					  $query->textfield(-name=>'uth_occ_eval',
							    -default=>$default{uth_occ_eval},
							    -size=>5)
					 ]),
					  
				### Threshold on Sig of the score
			      $query->td(['Occ sig',
					  $query->textfield(-name=>'lth_occ_sig',
							    -default=>$default{lth_occ_sig},
							    -size=>5),
					  $query->textfield(-name=>'uth_occ_sig',
							    -default=>$default{uth_occ_sig},
							    -size=>5)
					 ]),

			      ### Threshold on rank
			      $query->td(['Rank',
					  $query->textfield(-name=>'lth_occ_sig_rank',
							    -default=>$default{lth_occ_sig_rank},
							    -size=>5),
					  $query->textfield(-name=>'uth_occ_sig_rank',
							    -default=>$default{uth_occ_sig_rank},
							    -size=>5)
					 ]),
			     ]
			    )
		 );
  ## CRERs
  my $thresh_crer =
    $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
		  $query->Tr({-align=>center,-valign=>MIDDLE},
			     [
			      $query->th([" <A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A> "]),


			      ### Threshold on score
			      $query->td(['CRER size<b>*</b>',
					  $query->textfield(-name=>'lth_crer_size',
							    -default=>$default{lth_crer_size},
							    -size=>5),
					  $query->textfield(-name=>'uth_crer_size',
							    -default=>$default{uth_crer_size},
							    -size=>5)
					 ]),

				### Threshold on P-value of the score
			      $query->td(['site P-value<b>*</b>',
					  $query->textfield(-name=>'lth_site_pval',
							    -default=>$default{lth_site_pval},
							    -size=>5),
					  $query->textfield(-name=>'uth_site_pval',
							    -default=>$default{uth_site_pval},
							    -size=>5)
					 ]),

				### Threshold on crer_sites
			      $query->td(['CRER sites',
					  $query->textfield(-name=>'lth_crer_sites',
							    -default=>$default{lth_crer_sites},
							    -size=>5),
					  $query->textfield(-name=>'uth_crer_sites',
							    -default=>$default{uth_crer_sites},
							    -size=>5)
					 ]),

				### Threshold on crer_pval
			      $query->td(['CRER pval',
					  $query->textfield(-name=>'lth_crer_pval',
							    -default=>$default{lth_crer_pval},
							    -size=>5),
					  $query->textfield(-name=>'uth_crer_pval',
							    -default=>$default{uth_crer_pval},
							    -size=>5)
					 ]),
			      ### Threshold on crer_pval
			      $query->td(['CRER sig',
					  $query->textfield(-name=>'lth_crer_sig',
							    -default=>$default{lth_crer_sig},
							    -size=>5),
					  $query->textfield(-name=>'uth_crer_sig',
							    -default=>$default{uth_crer_sig},
							    -size=>5)
					 ]),
			      $query->Tr({-align=>middle,-valign=>TOP},
					 [
					  $query->td({-colspan=>4},[ "<b>*</b> =mandatory field"],
						    )
					 ]),
			     ]
			    )
		 );


  #############################################
  ## Table

  print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>3},
		      "<tr><td/>",
		      "<th bgcolor='#CCCCCC'>
  			<INPUT TYPE='radio' NAME='analysis_type' VALUE='analysis_sites' $checked{'analysis_sites'}><BR/>
  			<A HREF='help.matrix-scan.html#return_fields'>Individual matches</A></th>",
		      "<th bgcolor='#D6EEFA'>
  			<INPUT TYPE='radio' NAME='analysis_type' VALUE='analysis_crer' $checked{analysis_crer}><BR/>
			<A HREF='help.matrix-scan.html#return_fields'>CRERs <BR/> (Cis-Regulatory element <BR>Enriched Regions)</A> </th>",
		      "<th bgcolor='#F6E6CA'>
  			<INPUT TYPE='radio' NAME='analysis_type' VALUE='analysis_occ' $checked{analysis_occ}><BR/>
  			<A HREF='help.matrix-scan.html#return_fields'>Enrichment of hits<br>in the whole input sequence set</A></th> ",
		      "</tr>",
		      "<tr align='left' valign='top'><td><b>Fields to <BR/> return</b></td>",
		      "<td bgcolor='#CCCCCC'>$boxes_matches</td>",
		      "<td bgcolor='#D6EEFA'> $boxes_crer </td> ",
		      "<td bgcolor='#F6E6CA'>$boxes_occ</td>",
		      "</tr>",

		      $query->Tr({-align=>middle,-valign=>TOP},
				 [
				  $query->td({-colspan=>4},[ "<b>Other fields to return</b>  $boxes_add"],
					    )
				 ]),
		      "<tr align='left' valign='top'><td><b>Thresholds</b></td>",
		      "<td bgcolor='#CCCCCC'>$thresh_matches</td>",
		      "<td bgcolor='#D6EEFA'>$thresh_crer</td> ",
		      "<td bgcolor='#F6E6CA'>$thresh_occ</td>",
		      "</tr>",
		     );

}





