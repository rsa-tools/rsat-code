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

$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc

$default{bg_method}="bginput";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = 0;
$default{window_size} = 200;
$default{organism} = "Saccharomyces cerevisiae";
$default{matrix_format} = "tab";
$default{pseudo_counts} = 1;
$default{consensus_as_name} = "CHECKED";
$default{pseudo_distribution} = "bginput";
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{decimals} = "1";


## Return fields
## matches
$default{return_sites} = "CHECKED";
$default{return_pval} = "CHECKED";
$default{return_limits} = "CHECKED";
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
$default{return_crer_sites} = "";

$default{analysis_type} = "analysis_sites";
$checked{$default{analysis_type}} = "CHECKED";


## Threshold values
$default{lth_score} = "0";
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
$default{uth_pval} = "1e-5";

$default{lth_inv_cum} = "none";
$default{uth_inv_cum} = "none";
$default{lth_exp_occ} = "none";
$default{uth_exp_occ} = "none";
$default{lth_occ_pval} = "none";
$default{uth_occ_pval} = "none";
$default{lth_occ_eval} = "none";
$default{uth_occ_eval} = "none";
$default{lth_occ_sig} = "none";
$default{uth_occ_sig} = "none";
$default{lth_occ_sig_rank} = "none";
$default{uth_occ_sig_rank} = "3";

$default{lth_crer_size} = "none";
$default{uth_crer_size} = "200";
$default{lth_crer_sites} = "none";
$default{uth_crer_sites} = "none";
$default{lth_crer_pval} = "none";
$default{uth_crer_pval} = "none";
$default{lth_crer_sig} = "3";
$default{uth_crer_sig} = "none";




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
  if ($key eq "bg_method"){
  	$checked{$query->param($key)} = "CHECKED";
  }
  if ($key eq "analysis_type"){
  	$checked{$query->param($key)} = "CHECKED";
  }
}

&ReadMatrixFromFile();

### print the form ###
&RSA_header("matrix-scan");
#&ListParameters;

### head
print "<CENTER>";
print "Scan a DNA sequence with a profile matrix<BR>\n";
print "Program developed by <A HREF='mailto:jturatsi\@scmbb.ulb.ac.be (Jean Valery Turatsinze)'>Jean Val&eacute;ry Turatsinze</A>, <A HREF='mailto:morgane\@scmbb.ulb.ac.be (Morgane Thomas-Chollier)'>Morgane Thomas-Chollier</A> and <A HREF='mailto:jvanheld\@scmbb.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A><P>";

print "</CENTER>";

print $query->start_multipart_form(-action=>"matrix-scan.cgi");

################################################################
#### sequence
print "<hr>";
&DisplaySequenceChoice();


################################################################
#### Matrix specification
print "<hr>";
&GetMatrix();

################################################################
## Background model
print "<hr>";
&GetBackgroundModel("1","1");
			
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

$demo_sequence =">eve-RA	eve; upstream from -5500 to -1; size: 5500; location: chromosome:BDGP4.3:2R:1:20766785:1 5485728 5491227 D; upstream neighbour: CG12134-PB (distance: 6673)
CGTTCTCCAGCCGAAGATTTTTTCGAGCAACCAAAATATTATGGTGTGCCCCGCTGTTCT
CGCACAGTCAGCGCGAATTTGCTGCGGTGAGTCGATGCTGTTTCGCAGGACCTTCTTCCA
TTTTCGTCTCCCTCTTGCTCAGCCTGTCCCTGTTCCTCTGCAGTTCCCTATCTCCTGATG
CCTGTGCTCCTTTGGCGGCACTGTGTCCTGTCGTCGTTGTTTTCCTGTGATTTGACATGT
CTGTTAGCAGGATGCCTGACCCTGAGGCCGAGCCCTGGTCTCAGTGTCCACTGTTCCACT
TTGATGTGATTCGTCAGTGCGGTGGACTACTGCTACTGCTCTCTTGCTGGACTGCGTCTT
GAGTCCTGTTCGGCTGCCCCCTCCCGTGACCTCTGACCCTGCACTCTGCGGCTTTCCAGC
GGCGTTTGTTGGCGAATCTGACCCCGAGCTCCTGCTGCTCCTTCGCTCCTTCGCTCCTTC
TCCGCATCTCCGCTCTTTGGACTTCGTACGAATCAAAATTGGTCACAGCACCGAGTGAAT
TGCCCCGGAGACCGCAATGCGCTGTATTTATAGTAAACGTGTCCGATTGATTTGGCCACC
CGTGGCGGCTCTGTCACAGATGCCTCAATTTGCATCTATCGAATGGTTTACATGGCTCTA
AAAAGGTACCTCGATGGGTTGGTCACAATGTGGTGGCCTCTCAACATTGCAAGGCTCTTA
CTTGTGAATTATTAAGTTATTAACTGCTGCGATGTAAGTCATGGCAGTTTCTGTTTTCTT
TATAGGATATATATAGGAAGGATTAAAGGAGGCATGTACAATAATATGAGTATGATTTAG
CTCAAATTCCAAATATGATAAAAGTACAAAGCATACGATAATATAATCAAATTACGCTGA
CAATCACGATAATGTTCTTGTAGTAGTATTTGTGTAATATTTATGTTTTTTTAAGATAAG
AAACGGTAATAAAATCCACGTAAGTGTAAAAAATGGATGCCCTAATCTATGCCATGATGT
GTTCTACTTTCGAGATTTCGCCTCTGCCCTCATTGATGGTTTCCCGGGGCTACTTGGCCC
AAAAATCCCGGCCGTCCAAAAAGACGATCCTTAAAAAAGAAACCGCTAATCATTGGGCCG
CACAAAGAGCGGACAATCGCTCACCTAATTATTTGGCCCGATTGTGAGGAGCGGACAGTC
GGCTCGTGGACGCTTTTTGTGGCCTCTTTTTGTTTCGACAAAAAGCGAGCCAATTTTTTT
TCTTTCTGGGCCACTTTGTTGCTCTTTTTATGAGTTTTTTCCATTGTCAGTTTTTCCGGG
CCTGTCTCGCAGCCCTCGATTCCCGCGATGCCTGCCCTACAAACCTCCTAATTACGGCAG
TTAGTCGTTGTCCGGGACAGGAGAGTATGCGGAAGGACATGCGTGAGTTTATTGCCCGCT
CGAATTTCCACTAAAAATTGGGCCGAAAAAAAAACAACTAGGTAGGACTAGGAACTGCAA
ACTAGCAAAGCGGACGCGCCTTTTTATTGGTGCACCTTCGGCGGAACCGCAGGATAACAG
CAGTAAAAGCGACGACGAGGACACAAGGATCCTCGAAATCGAGAGCGACCTCGCTGCATT
AGAAAACTAGATCAGTTTTTTGTTTTGGCCGACCGATTTTTGTGCCCGGTGCTCTCTTTA
CGGTTTATGGCCGCGTTCCCATTTCCCAGCTTCTTTGTTCCGGGCTCAGAAATCTGTATG
GAATTATGGTATATGCAGATTTTTATGGGTCCCGGCGATCCGGTTCGCGGAACGGGAGTG
TCCTGCCGCGAGAGGTCCTCGCCGGCGATCCTTGTCGCCCGTATTAGGAAAGTAGATCAC
GTTTTTTGTTCCCATTGTGCGCTTTTTTCGCTGCGCTAGTTTTTTTCCCCGAACCCAGCG
AACTGCTCTAATTTTTTAATTCTTCACGGCTTTTCATTGGGCTCCTGGAAAAACGCGGAC
AAGGTTATAACGCTCTACTTACCTGCAATTGTGGCCATAACTCGCACTGCTCTCGTTTTT
AAGATCCGTTTGTTTGTGTTTGTTTGTCCGCGATGGCATTCACGTTTTTACGAGCTCGTT
CCTTCGGGTCCAAAATTATGCCAGTTTGTTTTGTCTCTGGCAATTATTGGAAATTTCATT
GGGTCGATTTCGCTGCCTTCCTTGCTCTTCCCTTGAGAAAAGTGAATAGGTTGTGCCATA
AAAATCGCTGCTCCTGAAGACCAAATGAAATGGATTTGTGTAAGCATTAAAAACGCGAGG
CAAGCCCCAAGATTCCTCCACTGCTTTTTTTATATTGCCCACTGCTAAATGCAGCTAATT
CGTCGATTGTTTAAAAATTAAATTACTTATGTTGCCATTCATACATCCCCTCACATTTTA
TGGCCATTTGAGTGCGGGGTGCACAGTTCTGTCTTAAGTGGCGGATGGAAACCACCACAT
TTACTCGAGGGATGATGTGCTCTAATATCTCCTCATCAAATGGGATGGTTTCTATGGAAA
GGCAAAATCGTTGTAAAGTGAGGCGGAGTTAAAAAATACCTTGTTATAGCCTTTTTAAAA
TAACACAAGATCGTTCGAATTGACTAGAAATATCAAAGTCTTTTTGTATTGAAGCGAGTG
TAGTCTCAATTTATGCTTAATTTTAAGAAATACATCTCTTTATTAGCCCCAAAATGAAAC
AAATGGTCTACTAATTAAGCAAGTCAACAGAATTTTTATGCAATTATTCAAAATGAAATA
ATATATACATAAGATGTTTTTGGGAATCTGTCATGGGGTTTCTGAAATAGGTTTGCCAAA
CAAATTTTAAGTATAAATGTATACATATGTCAACTAATAAATTTAGCAAATAAAATGTAC
CTGCAAGTATCTATAAATTTATTGGACCAATTTTGTGTAAAAAACTGAACTGGCACTCTT
CCCAAGAATGGGACTTCGAGGACTCCTTGCTGAATCACTTACTCAACCCATTCCAACTCA
TCCAATCCGCGCAATCATCATAAATTTTGGCCTTTTTGTTGTAATTGTTTTATGGCAGAA
ATTACTCAATCATCAAGCATAATTCCCTCGTTTTCGCCGTTTTATTGCCAATTTTTGCAC
TGCCTTTGCCTTTTTCCCGCCCTTTCCTCAGCGTTTTGCGAATCTTTGCCGGCATTTCTA
TTGCGCGGACAATCCGGCCAGTGTGTTGGCCATTTACTTGCCATGATGACGGGCATAATC
AGCGAGATCGGCGCTTTGTGAGTGCAGAATGTGCAATAAAGCGGCAACAATCGGCAGGGA
TTCGCCTTCCCATATTCCGGGTATTGCCGGCCCGGGAAAATGCGAAAGTGTTTGCGGATC
GAGATGGAAGATAGAGGATTGAGTATTGAAACGAGGAAGGTACTTCCGCCGGCGGACACT
TTCGCCTAACCAAGCCAATCCAACCCATCCCAATCCAATCCAACCCACCCGATCGCCATA
AAGGGTATTTACTGTCGCTGCCGCAGAGCCTCGCTTGACGACTTAACCCAAGCGGTCGTT
TCGCGTCCATTCTCCGGACGGAGTCAAAGACAAAGGCCGGCGGAGGTGGACAATAGGCAA
GGTTGTTGCTTGTGGGTAGGGTTTGAGCTATGAGCTATGAGCTGTGAGCTGTTAGCCCTG
AACCCCGAACCTCGAGAATTGAACCTTTCCCGGGGCAAGAAGGCTTGCATGTGGGCCTTT
TCCAGGTCGGCCAGTAGGTAGAGTTGTTGCGATGCGGCTATGCCGGGCGAGTTAATGCCA
ATGCAAATTGCGGGCGCAATATAACCCAATAATTTGAAGTAACTGGCAGGAGCGAGGTAT
CCTTCCTGGTTACCCGGTACTGCATAACAATGGAACCCGAACCGTAACTGGGACAGATCG
AAAAGCTGGCCTGGTTTCTCGCTGTGTGTGCCGTGTTAATCCGTTTGCCATCAGCGAGAT
TATTAGTCAATTGCAGTTGCAGCGTTTCGCTTTCGTCCTCGTTTCACTTTCGAGTTAGAC
TTTATTGCAGCATCTTGAACAATCGTCGCAGTTTGGTAACACGCTGTGCCATACTTTCAT
TTAGACGGAATCGAGGGACCCTGGACTATAATCGCACAACGAGACCGGGTTGCGAAGTCA
GGGCATTCCGCCGATCTAGCCATCGCCATCTTCTGCGGGCGTTTGTTTGTTTGTTTGCTG
GGATTAGCCAAGGGCTTGACTTGGAATCCAATCCCGATCCCTAGCCCGATCCCAATCCCA
ATCCCAATCCCTTGTCCTTTTCATTAGAAAGTCATAAAAACACATAATAATGATGTCGAA
GGGATTAGGGGCGCGCAGGTCCAGGCAACGCAATTAACGGACTAGCGAACTGGGTTATTT
TTTTGCGCCGACTTAGCCCTGATCCGCGAGCTTAACCCGTTTTGAGCCGGGCAGCAGGTA
GTTGTGGGTGGACCCCACGATTTTTTTGGCCAAACCTCCAAGCTAACTTGCGCAAGTGGC
AAGTGGCCGGTTTGCTGGCCCAAAAGAGGAGGCACTATCCCGGTCCTGGTACAGTTGGTA
CGCTGGGAATGATTATATCATCATAATAAATGTTTTGCCCAACGAAACCGAAAACTTTTC
AAATTAAGTCCCGGCAACTGGGTTCCCATTTTCCATTTTCCATGTTCTGCGGGCAGGGGC
GGCCATTATCTCGCTACAGCAGTTCCCAAATGGTTATGGCTGGACACCCCTGCCGCCGCT
CCAACGGGGTGGATGAAGCCCCCAAAACCCGAAAGTCATGGCAGCCATGGCAGTGTGGGG
CTGTTAAACGTGCGGCATAATATTAAGACTTCATAAAAGCGCAAATAATTCGCTGGCAGG
CGATCGATAATACATACATACAAATATATAGTGGGATACACACACTCTCTGCCGGCAAAC
ACACACCACCCGACCCGACTGAGCGGCATAATGCCATATCATTCTTGATGAAGCCGATAA
AATCCCATTATTAAGGGGGCCCGCCCGTCCCGCTCGCTCCTGCGGAGCAACCGCCTGCGG
GCGGGCGAGACAAAAGATTCGCTCATCCGCTATGAATACCAAATCGGAACTCTCTCTCTC
TCCAGCTCGGGAGTGCCATGGCCAGCATGGCCAGGACCTCCTCATGGTCCTGCCGAGCAG
AGAACGCGGCTCCATCCCGCTGCTCCGGGTCCTGCTCCTCCGCTTTGTCCCGCCTCGTTA
TCGCCGCTCAGCACCGAGAGCACAGCAGCGCATCCACTCTCAGCACCGCACGATTAGCAC
CGTTCCGCTCAGGCTGTCCCGCTCGCACCTGCCTGGGTCGCTGCGATTGGCCGCTCCCAG
CGACGGCGGCCATTTGCCTGCAGAGCGCAGCGGTATAAAAGGGCGCGGGGTGGCTGAGAG
CAGCACACTCGAGCTGTGACCGCCGCACAGTCAACAACTAACTGCCTTCGTTAATATCCT
CTGAATAAGCCAACTTTGAATCACAAGACGCATACCAAAC";

$demo_matrix = "
; Kr matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       30      37      0       3       1       5       4       2
c       |       4       0       35      37      41      9       1       4
g       |       4       2       3       0       0       11      7       0
t       |       6       5       6       4       2       19      32      38
//
; Med matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       0       0       0       0       0       0
c       |       6       0       0       10      0       0
g       |       1       10      2       0       0       7
t       |       3       0       8       0       10      3
//
; Stat92E matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       0       0       0       0       0       0       0       3       3
c       |       0       0       3       0       2       0       0       0       0
g       |       0       0       0       3       1       3       3       0       0
t       |       3       3       0       0       0       0       0       0       0
//
; bcd matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       8       45      46      1       0       3
c       |       3       1       0       0       45      27
g       |       1       1       2       13      1       3
t       |       36      1       0       34      2       15
//
; eve matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       0       5       3       6       1       2       7       0       0       8       4       2       3       1       0
c       |       8       0       4       2       3       3       0       0       1       1       0       5       0       6       7
g       |       1       4       2       1       0       3       2       0       0       0       4       2       6       2       0
t       |       0       0       0       0       5       1       0       9       8       0       1       0       0       0       2
//
; gt matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       3       0       5       0       0       6       0       0       0       5
c       |       0       0       0       0       0       2       5       0       3       0
g       |       0       1       2       0       5       0       2       4       0       0
t       |       5       7       1       8       3       0       1       4       5       3
//
; hb matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       1       4       0       0       2       53      24      2
c       |       0       4       2       0       2       4       14      20
g       |       3       4       0       0       0       22      10      55
t       |       99      91      101     103     99      24      55      26
//
; kni matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       6       18      7       0       3       15      10      10
c       |       2       7       6       33      4       7       8       2
g       |       22      0       2       0       4       4       10      2
t       |       3       8       18      0       22      7       5       19
//
; pan matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       15      7       1       22      22      21      9       7
c       |       1       2       23      0       1       1       2       6
g       |       3       0       0       0       1       0       10      9
t       |       6       16      1       3       1       3       4       3
//
; prd matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       0       5       0       0       3       0       9       7
c       |       6       0       8       5       2       7       0       2
g       |       0       4       1       4       4       2       0       0
t       |       3       0       0       0       0       0       0       0
//
; tin matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       10      0       0       0       2       10
c       |       1       11      0       0       1       0
g       |       0       0       0       0       8       1
t       |       0       0       11      11      0       0
//
; ttk matrix, derived from OregAnno by J.V. Turatsinze (2007).
a       |       3       7       0       0       8       0       3
c       |       5       0       0       0       0       5       5
g       |       0       1       8       8       0       0       0
t       |       0       0       0       0       0       3       0
";

print "<TD><B>";
print $query->hidden(-name=>'bg_method',-default=>'bgfile');
print $query->hidden(-name=>'uth_pval',-default=>'1e-4');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'organism',-default=>'Drosophila_melanogaster_EnsEMBL');
print $query->hidden(-name=>'return_rank',-default=>'');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->submit(-label=>"DEMO 1");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B>";
print $query->start_multipart_form(-action=>"matrix-scan_form.cgi");
print $query->hidden(-name=>'bg_method',-default=>'bgfile');
print $query->hidden(-name=>'uth_pval',-default=>'1e-4');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'organism',-default=>'Drosophila_melanogaster_EnsEMBL');
print $query->hidden(-name=>'analysis_type',-default=>'analysis_crer');
print $query->hidden(-name=>'return_rank',-default=>'');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->submit(-label=>"DEMO 2");
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
print "<TD><B><A HREF='tutorials/tut_matrix-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

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
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }
  $boxes_matches .= "<BR/>";
   @return_fields_matches = qw( limits normw);
  foreach my $field (@return_fields_matches) {
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }
  $boxes_matches .= "<BR/>";
   @return_fields_matches = qw(weight_limits bg_residues);
  foreach my $field (@return_fields_matches) {
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }
  
  ### Return fields
  my $boxes_occ = "";
  @return_fields_occ = qw(distrib);
  foreach my $field (@return_fields_occ) {
    $boxes_occ .= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }
  $boxes_occ .= "<BR/>";
    @return_fields_occ = qw(occ_proba);
  foreach my $field (@return_fields_occ) {
    $boxes_occ .= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
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
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }
	$boxes_crer .= "<BR/>";
  @return_fields_crer = qw(limits);
  foreach my $field (@return_fields_crer) {
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>' '.$field.' ');
  }
  @return_fields_crer = qw(crer_sites);
  foreach my $field (@return_fields_crer) {
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>' '."sites".' ');
  }
  
  
  
  ### Return fields
  my $boxes_add = "";
  @return_fields_add = qw(matrix freq_matrix weight_matrix bg_model);
  foreach my $field (@return_fields_add) {
    $boxes_add.= $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
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

## Occurences
	 my $thresh_occ = 
     $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
		      $query->Tr({-align=>center,-valign=>MIDDLE},
				 [
				  $query->th(["<A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					      " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					      " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A> "]),

				   $query->th({-colspan=>3,-align=>left},["observed distribution of scores (distrib)"
					     ]),
				  ### Threshold on score
				  $query->td(['score',
					      $query->textfield(-name=>'lth_score',
								-default=>$default{lth_score},
								-size=>5),
					      $query->textfield(-name=>'uth_score',
								-default=>$default{uth_score},
								-size=>5)
					     ]),
				
				### Threshold on occ_inv_cum
				  $query->td(['occ inv_cum',
					      $query->textfield(-name=>'lth_inv_cum',
								-default=>$default{lth_inv_cum},
								-size=>5),
					      $query->textfield(-name=>'uth_inv_cum',
								-default=>$default{uth_inv_cum},
								-size=>5)
					     ]),

				$query->th({-colspan=>3,-align=>left},["over-representation of matches (occ_proba)"
					     ]),
				
				### Threshold on exp_occ
				  $query->td(['exp_occ',
					      $query->textfield(-name=>'lth_exp_occ',
								-default=>$default{lth_exp_occ},
								-size=>5),
					      $query->textfield(-name=>'uth_exp_occ',
								-default=>$default{uth_exp_occ},
								-size=>5)
					     ]),

				  ### Threshold on P-value of the score
				  $query->td(['occ P-value',
					      $query->textfield(-name=>'lth_occ_pval',
								-default=>$default{lth_occ_pval},
								-size=>5),
					      $query->textfield(-name=>'uth_occ_pval',
								-default=>$default{uth_occ_pval},
								-size=>5)
					     ]),
				
				 ### Threshold on P-value of the score
				  $query->td(['occ E-value',
					      $query->textfield(-name=>'lth_occ_eval',
								-default=>$default{lth_occ_eval},
								-size=>5),
					      $query->textfield(-name=>'uth_occ_eval',
								-default=>$default{uth_occ_eval},
								-size=>5)
					     ]),
					     
				### Threshold on Sig of the score
				  $query->td(['occ Sig',
					      $query->textfield(-name=>'lth_occ_sig',
								-default=>$default{lth_occ_sig},
								-size=>5),
					      $query->textfield(-name=>'uth_occ_sig',
								-default=>$default{uth_occ_sig},
								-size=>5)
					     ]),

				  ### Threshold on rank
				  $query->td(['rank',
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
				  $query->td(['crer_size<b>*</b>',
					      $query->textfield(-name=>'lth_crer_size',
								-default=>$default{lth_crer_size},
								-size=>5),
					      $query->textfield(-name=>'uth_crer_size',
								-default=>$default{uth_crer_size},
								-size=>5)
					     ]),
					     
				### Threshold on P-value of the score
				  $query->td(['site P-value<b>*</b>',
					      $query->textfield(-name=>'lth_pval',
								-default=>$default{lth_pval},
								-size=>5),
					      $query->textfield(-name=>'uth_pval',
								-default=>$default{uth_pval},
								-size=>5)
					     ]),
				
				### Threshold on crer_sites
				  $query->td(['crer_sites',
					      $query->textfield(-name=>'lth_crer_sites',
								-default=>$default{lth_crer_sites},
								-size=>5),
					      $query->textfield(-name=>'uth_crer_sites',
								-default=>$default{uth_crer_sites},
								-size=>5)
					     ]),
					     
				### Threshold on crer_pval
				  $query->td(['crer_pval',
					      $query->textfield(-name=>'lth_crer_pval',
								-default=>$default{lth_crer_pval},
								-size=>5),
					      $query->textfield(-name=>'uth_crer_pval',
								-default=>$default{uth_crer_pval},
								-size=>5)
					     ]),
								### Threshold on crer_pval
				  $query->td(['crer_sig',
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
  			<INPUT TYPE='radio' NAME='analysis_type' VALUE='analysis_occ' $checked{analysis_occ}><BR/>
  			<A HREF='help.matrix-scan.html#return_fields'>Enrichment in sites <BR/> (over-representation of matches)</A></th> ",
			"<th bgcolor='#F6E6CA'>
			<INPUT TYPE='radio' NAME='analysis_type' VALUE='analysis_crer' $checked{analysis_crer}><BR/>
			<A HREF='help.matrix-scan.html#return_fields'>CRERs <BR/> (Cis-Regulatory element <BR>Enriched Regions)</A> </th>",
				"</tr>",
			"<tr align='left' valign='top'><td><b>Fields to <BR/> return</b></td>",
  			"<td bgcolor='#CCCCCC'>$boxes_matches</td>",
  			"<td bgcolor='#D6EEFA'>  $boxes_occ</td> ",
			"<td bgcolor='#F6E6CA'>$boxes_crer</td>",
				"</tr>",

			$query->Tr({-align=>middle,-valign=>TOP},
				 [
				  $query->td({-colspan=>4},[ "<b>Other fields to return</b>  $boxes_add"],
				  			)
				 ]),
			"<tr align='left' valign='top'><td><b>Thresholds</b></td>",
  			"<td bgcolor='#CCCCCC'>$thresh_matches</td>",
  			"<td bgcolor='#D6EEFA'>  $thresh_occ</td> ",
			"<td bgcolor='#F6E6CA'>$thresh_crer</td>",
				"</tr>",
		     );

}





