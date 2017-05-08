#!/usr/bin/perl
#### this cgi script fills the HTML form for the program matrix-quality
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "matrix_web_forms.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
use RSAT::matrix;
use RSAT::MatrixReader;

### Read the CGI query
$query = new CGI;

local @supported_input_formats = sort(keys( %RSAT::MatrixReader::supported_input_format));
local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

################################################################
### default values for filling the form
$default{demo_descr}="";
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "transfac";
$default{kfold}="none";
$default{permutation1} = "1";
$default{sep_perm1} = "";
$default{permutation2} = "1";
$default{sep_perm2} = "";
$checked{$default{nwd}} ="";
$default{tag1} = "sequence_set1";
$default{tag2} = "sequence_set2";
$default{pseudo_prior} = "pseudo_prior";
$default{pseudo_counts}="1";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{bg_method}="bgfile";
$checked{$default{bg_method}} = "CHECKED";
$default{organism}="Escherichia_coli_GCF_000005845.2_ASM584v2";
#$default{html_title}="";
$default{markov_order} = "0";
$default{m_sites}="1";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

&ListParameters() if ($ENV{rsat_echo} >=2);

#&ReadMatrixFromFile();

################################################################
### print the form ###


################################################################
### header
&RSA_header("matrix-quality", "form");
print "<center>";
print "Evaluate the quality of a Position-Specific Scoring Matrix (PSSM), by
    comparing score distributions obtained with this matrix in various
    sequence sets.</p>\n";
print "The most classical use of the program is to compare score distributions
    between <em>positive</em> sequences (e.g. true binding sites for the considered
    transcription factor) and <em>negative</em> sequences (e.g. intergenic
    sequences between convergently transcribed genes).<p>\n";
print "<p>Program developed by <a target='_top' href='http://www.ccg.unam.mx/ccg-OrganicG/personalInfo?idPersona=253'>Alejandra Medina Rivera</a>, \n";
print " <a target='_top' href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</A>,\n";
print "and <a target='_top' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</A>.</p>\n";
print "</center>\n";
print "<b>Citation</b>: Medina-Rivera, A., Abreu-Goodger, C., Salgado-Osorio, H., Collado-Vides, J. and van Helden, J. (2010). Empirical and theoretical evaluation of transcription factor binding motifs. Nucleic Acids Res. 2010 Oct 4. [Epub ahead of print] <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/20923783'>[Pubmed 20923783]</a> <a target='_blank' href='http://nar.oxfordjournals.org/content/early/2010/10/04/nar.gkq710.full.pdf'>[Full text]</a>.";


## demo description
print $default{demo_descr};
print $query->start_multipart_form(-action=>"matrix-quality.cgi");


################################################################
#### Matrix specification
print "<hr>";
print "<h2 style='margin-left: 50px;'> Title ";

print $query->textfield(-name=>'html_title',
			 -default=>$default{html_title},
			 -size=>30) ."</h2>";

print "<fieldset>
<legend><b><a href='help.convert-matrix.html#io_format'>1 - Matrix </a></b></legend>";


&GetMatrix();
print "<p></p>";
print $query->checkbox(-name=>'matrix_sites',
  		       -checked=>$default{m_sites},
		       -label=>'');

print "&nbsp;Matrix file includes sites";
print "<p><font color='orange'>Only the first matrix will be taken in acount</font></p>";

print "<\p><b>K fold validation</B>&nbsp;";
print $query->popup_menu(-name=>'kfold',
			 -Values=>["none",0,3,4,5,6,7,8,9,10],
			 -default=>$default{kfold});
print "&nbsp;"x5, "<font color='orange'><b>Note:</b> validation requires a matrix with binding sites, in a suitable format (e.g. transfac, meme).</font>";

print "</fieldset><p/>";


################################################################
#### Sequence specification

print "<fieldset>
<legend><b><a href='help.formats.html'>2 - Sequences </a></b></legend>";


print "<h2> Mandatory Sequence </h2>";

&SeqBoxMQ(1);
print "<hr>";

print "<h2> Optional Sequence </h2>";
&SeqBoxMQ(2);

print "</fieldset><p/>";

################################################################
#### Background specifiaction

print "<fieldset>
<legend><b><a href='help.matrix-scan.html#markov_order'>3 - Background </a></b></legend>";
my %bg_params =(
    "markov" => 1,
    "markov_message" => 1
    );
&GetBackgroundModel(%bg_params);

#print "<br/>Note: Only Bernoulli models are supported. Higher-order Markov models are converted into Markov 0 (Bernoulli).";
print "</fieldset><p/>";



################################################################
## Send results by email or display on the browser
print "<p>\n";
&SelectOutput("email");


################################################################
## Action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"matrix-quality_form.cgi");
$demo_html_title=" LexA and CTCF matrices from RegulonDB 2015";
$demo_matrix="AC  LexA.RegulonDB2012
XX
ID  LexA.RegulonDB2012
XX
DE  waCTGtatawAwahmCAGya
P0       A     C     G     T
1       16     6     1    17
2       26     4     7     3
3        1    39     0     0
4        2     1     0    37
5        0     1    39     0
6        6     3     8    23
7       30     5     0     5
8        1     3     5    31
9       31     1     6     2
10      16     4     6    14
11      33     0     2     5
12      18     7     1    14
13      25     1     5     9
14      16     9     3    12
15      23    16     1     0
16       1    38     1     0
17      37     1     0     2
18       4     0    35     1
19       2    12     7    19
20      21     4     7     8
XX
BA  40 sequences
XX
BS  GACTGTATAAAACCACAGCC; site_0; 1; 20; 0; p
BS  TCCTCTCTGCAATACGCAGA; site_1; 1; 20; 0; p
BS  AGCTGAATAAATATACAGCA; site_2; 1; 20; 0; p
BS  CACTGTATACTTTACCAGTG; site_3; 1; 20; 0; p
BS  ACCAGAATATACATAATAGT; site_4; 1; 20; 0; p
BS  CACTGTATAAATAAACAGCT; site_5; 1; 20; 0; p
BS  AACTGATTAAAAACCCAGCG; site_6; 1; 20; 0; p
BS  AACTGTCGATAAGCGCAGCC; site_7; 1; 20; 0; p
BS  CGCTGCCCTGAAAGCCAGGC; site_8; 1; 20; 0; p
BS  TACTGGATAAAAAAACAGTT; site_9; 1; 20; 0; p
BS  TACTGTATAAACAGCCAATA; site_10; 1; 20; 0; p
BS  TGCTGTATGGATTAACAGGA; site_11; 1; 20; 0; p
BS  TACTGTATGGATGTACAGTA; site_12; 1; 20; 0; p
BS  ACCTGTATAAATAACCAGTA; site_13; 1; 20; 0; p
BS  ATCTGCTGGCAAGAACAGAC; site_14; 1; 20; 0; p
BS  TACTGTATATAAAAACAGTA; site_15; 1; 20; 0; p
BS  CACTGTATAAAAATCCTATA; site_16; 1; 20; 0; p
BS  AACTGTCGATACGTACAGTA; site_17; 1; 20; 0; p
BS  CACTGGATAGATAACCAGCA; site_18; 1; 20; 0; p
BS  TACTGTATAAAATCACAGTT; site_19; 1; 20; 0; p
BS  TGCTGGATAGATATCCAGCG; site_20; 1; 20; 0; p
BS  AACTGAAAAATGGCACAGTA; site_21; 1; 20; 0; p
BS  TACCGGATCGTTTCACAGTA; site_22; 1; 20; 0; p
BS  TACTGTATATAAAAACAGTA; site_23; 1; 20; 0; p
BS  AACTGGATAAAATTACAGGG; site_24; 1; 20; 0; p
BS  TACTGTATATAAAACCAGTT; site_25; 1; 20; 0; p
BS  TACTGTACACAATAACAGTA; site_26; 1; 20; 0; p
BS  TACTGTATGAGCATACAGTA; site_27; 1; 20; 0; p
BS  CTATGTTTATATAACCATCA; site_28; 1; 20; 0; p
BS  AGCTGGCGTTGATGCCAGCG; site_29; 1; 20; 0; p
BS  AACTGGATAATCATACAGTA; site_30; 1; 20; 0; p
BS  TGCAGGTTATAAAACCAGCA; site_31; 1; 20; 0; p
BS  AACTGTATATAAATACAGTT; site_32; 1; 20; 0; p
BS  TACTGTATAAATAAACAGTA; site_33; 1; 20; 0; p
BS  ATCTGTATATATACCCAGCT; site_34; 1; 20; 0; p
BS  AACTGCACAATAAACCAGAG; site_35; 1; 20; 0; p
BS  TGCTGTATATACTCACAGCA; site_36; 1; 20; 0; p
BS  AACTGTATATACACCCAGGG; site_37; 1; 20; 0; p
BS  ACCTGAATGAATATACAGTA; site_38; 1; 20; 0; p
BS  TACTGATGATATATACAGGT; site_39; 1; 20; 0; p
CC  program: consensus
CC  id: LexA.1nt_upstream-noorf-ovlp-2str.20.cons1
CC  matrix.nb: 1
CC  command: consensus -a /state/partition1/space14/PGC/RegulonDB/local/data/Matrices/data/bg_freqs/1nt_upstream-noorf_Escherichia_coli_K12-ovlp-2str.cons_bg -L 20 -d -pt 0 -pf 1 -c2
CC  sites: 40
CC  cons.unadjusted.information: 11.7602
CC  cons.adjusted.information: 10.9922
CC  cons.ln.Pval: -367.453
CC  cons.Pval: 2.6136E-160
CC  cons.ln.exp: -244.504
CC  cons.exp: 6.50303E-107
CC  matrix.nb: 1
CC  sites: 40
CC  consensus.strict: taCTGtataaAaaaaCAGta
CC  consensus.strict.rc: TACTGTTTTTTTATACAGTA
CC  consensus.IUPAC: waCTGtatawAwahmCAGya
CC  consensus.IUPAC.rc: TRCTGKDTWTWTATACAGTW
CC  consensus.regexp: [at]aCTGtata[at]A[at]a[act][ac]CAG[ct]a
CC  consensus.regexp.rc: T[AG]CTG[GT][AGT]T[AT]T[AT]TATACAGT[AT]
XX
//
AC  CRP.RegulonDB2012
XX
ID  CRP.RegulonDB2012
XX
DE  wtGtGatcyrsaTCACahwttt
P0       A     C     G     T
1      110    43    27    80
2       37    39    23   161
3       37    16   146    61
4       26    48    18   168
5       34     5   195    26
6      192    32    18    18
7       63    47    45   105
8       67    83    41    69
9       52    71    47    90
10     112    32    57    59
11      50    77    79    54
12     129    37    41    53
13      18    15    15   212
14      25   221    12     2
15     231     2    15    12
16      22   224     6     8
17     190    13    41    16
18      86    60    20    94
19      90    19     0   151
20      71    22    20   147
21      51    34    16   159
22      57    41    50   112
XX
BA  260 sequences
XX
BS  ATGTGATTCATATCACATATTT; site_0; 1; 22; 0; p
BS  CTGTGATTGGTATCACATTTTT; site_1; 1; 22; 0; p
BS  ATGTGATCCAGATCACATCTAT; site_2; 1; 22; 0; p
BS  ATCTGACTCACATCACACTTTT; site_3; 1; 22; 0; p
BS  ATGTTATCCACATCACAATTTC; site_4; 1; 22; 0; p
BS  ATTTGATACCCATCACACTTTC; site_5; 1; 22; 0; p
BS  ATGTGATTTTCATCACGATTTA; site_6; 1; 22; 0; p
BS  ATGTGATTAACAGCACATTTTT; site_7; 1; 22; 0; p
BS  ATGTGAGTAGTGTCACATTTTT; site_8; 1; 22; 0; p
BS  ATGTGATTTGCTTCACATCTTT; site_9; 1; 22; 0; p
BS  GTGTGAAGTTGATCACAAATTT; site_10; 1; 22; 0; p
BS  ATGTGATACAAATCACATAAAT; site_11; 1; 22; 0; p
BS  GTGTGATCGTCATCACAATTCG; site_12; 1; 22; 0; p
BS  ATGAGATTTTCATCACACATTT; site_13; 1; 22; 0; p
BS  AAGTGATTTAGATCACATAATA; site_14; 1; 22; 0; p
BS  GTGTGATCTGCATCACGCATTA; site_15; 1; 22; 0; p
BS  ATGTGTACGAAATCACATTTTT; site_16; 1; 22; 0; p
BS  TTTTGACATGTATCACAAATTT; site_17; 1; 22; 0; p
BS  TTGTGATTCAGATCACAAAGAT; site_18; 1; 22; 0; p
BS  ACGTGATCTTCATCACAAATAA; site_19; 1; 22; 0; p
BS  TCGTGATCAAGATCACATTCTC; site_20; 1; 22; 0; p
BS  ACGTGATGATGTTCACAATTTG; site_21; 1; 22; 0; p
BS  ATTTGAAGTAGCTCACACTTAT; site_22; 1; 22; 0; p
BS  ATTTGATTTAGATCGCAATTTG; site_23; 1; 22; 0; p
BS  TTGTGATTCGATTCACATTTAA; site_24; 1; 22; 0; p
BS  ATGTGATTGATATCACACAAAA; site_25; 1; 22; 0; p
BS  ATGCGATTCCACTCACAATATT; site_26; 1; 22; 0; p
BS  AAGTGACCGAAATCACACTTAA; site_27; 1; 22; 0; p
BS  ATGAGATTCAGATCACATATAA; site_28; 1; 22; 0; p
BS  GAGTGACGTAGATCACACTTAT; site_29; 1; 22; 0; p
BS  CTTCGATACACATCACAATTAA; site_30; 1; 22; 0; p
BS  ATTTGAAGTTCATCACACTTCA; site_31; 1; 22; 0; p
BS  CGGTGATCTATTTCACAAATTA; site_32; 1; 22; 0; p
BS  ATGCGATCTATATCACGCTGTG; site_33; 1; 22; 0; p
BS  ATTTGTTCCTCTTCACATTTTT; site_34; 1; 22; 0; p
BS  ATTTGCACGGCGTCACACTTTG; site_35; 1; 22; 0; p
BS  ATGTGAAATAAATCAAAATTTC; site_36; 1; 22; 0; p
BS  TAGTGATCCACGCCACATTTTG; site_37; 1; 22; 0; p
BS  TTGCGATCAAAATAACACTTTT; site_38; 1; 22; 0; p
BS  TTGTGAATCTTTTCACAGTTTA; site_39; 1; 22; 0; p
BS  ATTCAATATTCATCACACTTTT; site_40; 1; 22; 0; p
BS  AAGTGATGCAAATCACATAAAT; site_41; 1; 22; 0; p
BS  TCGTGAACTACGGCACACTTTG; site_42; 1; 22; 0; p
BS  TTATGACCCTCTTCACATTTCG; site_43; 1; 22; 0; p
BS  GTGTGATGCAAGCCACATTTTT; site_44; 1; 22; 0; p
BS  ATGTGAGCGAGATCAAATTCTA; site_45; 1; 22; 0; p
BS  ACGTGACGTTCATCACAAAACG; site_46; 1; 22; 0; p
BS  ATGTGAATTGCCGCACACATTA; site_47; 1; 22; 0; p
BS  ATGTGAGTTAGCTCACTCATTA; site_48; 1; 22; 0; p
BS  AGGTGAGATGCATCACGCTTCG; site_49; 1; 22; 0; p
BS  AGGTGACCGGTTTCACAAATAT; site_50; 1; 22; 0; p
BS  CTGTGAGTGATTTCACAGTATC; site_51; 1; 22; 0; p
BS  TTGTTATTAGTCTCACACTTTT; site_52; 1; 22; 0; p
BS  GTGCGAGCCAGCTCAAACTTTT; site_53; 1; 22; 0; p
BS  ATATGACGGCGGTCACACTTAT; site_54; 1; 22; 0; p
BS  ATGCGATCTGGTTCAAATAATT; site_55; 1; 22; 0; p
BS  TTTCGAGGTTGATCACATTTCC; site_56; 1; 22; 0; p
BS  ATGTGCGACCACTCACAAATTA; site_57; 1; 22; 0; p
BS  AGGTGAGAGCCATCACAAATGT; site_58; 1; 22; 0; p
BS  ATGTGAGCCAGCTCACCATAAA; site_59; 1; 22; 0; p
BS  ATGAGATCGAGCACACATTTTA; site_60; 1; 22; 0; p
BS  ATGTGAAGCAAATCACCCACTT; site_61; 1; 22; 0; p
BS  ATTTGAGTAAGTTCTCAATTTT; site_62; 1; 22; 0; p
BS  TTGTGAGCTTGCTCGCACTTCG; site_63; 1; 22; 0; p
BS  TTATGACGCTCTTCACACTCTG; site_64; 1; 22; 0; p
BS  ATGCGCGACGCATCGCAAATTT; site_65; 1; 22; 0; p
BS  ATATGACAACCATCACAAAAAT; site_66; 1; 22; 0; p
BS  CTGTAACAGAGATCACACAAAG; site_67; 1; 22; 0; p
BS  TAGTGAAGCAGATCGCATTATA; site_68; 1; 22; 0; p
BS  TATTGTCCCCGATCACACTTTT; site_69; 1; 22; 0; p
BS  CCATGATCCGCGCCACACTTTT; site_70; 1; 22; 0; p
BS  CAGTGATCCAGGTCACGATAAC; site_71; 1; 22; 0; p
BS  GTTAAATATAGATCACAATTTT; site_72; 1; 22; 0; p
BS  TATAGATCTCCGTCACATTTTT; site_73; 1; 22; 0; p
BS  ATGTGACTACCATCACTTTAAT; site_74; 1; 22; 0; p
BS  TTGCGAAGCGCGTCACTATTTA; site_75; 1; 22; 0; p
BS  ATATGACGGTGTTCACAAAGTT; site_76; 1; 22; 0; p
BS  ATTTGCACTGTGTCACAATTCC; site_77; 1; 22; 0; p
BS  TTGTGATCGTTATCTCGATATT; site_78; 1; 22; 0; p
BS  TTGCAAGCAACATCACGAAATT; site_79; 1; 22; 0; p
BS  GTGTTAAATTGATCACGTTTTA; site_80; 1; 22; 0; p
BS  TTGCGAGCGAGCGCACACTTGT; site_81; 1; 22; 0; p
BS  ATTTATTCCATGTCACACTTTT; site_82; 1; 22; 0; p
BS  TAGAGATCTACTTCACAAATCA; site_83; 1; 22; 0; p
BS  TTTTGCGCGAGGTCACTATTTT; site_84; 1; 22; 0; p
BS  TTGCGGGTCGCGTCACATTTAA; site_85; 1; 22; 0; p
BS  AGGTGATTTTGATCACGGAATA; site_86; 1; 22; 0; p
BS  TTTTGAATCCCATCACAAACCC; site_87; 1; 22; 0; p
BS  CTGTGGTTGCCATCACAGATAT; site_88; 1; 22; 0; p
BS  TTGTGGCCTGCTTCAAACTTTC; site_89; 1; 22; 0; p
BS  AAGTGAACCATATCTCAATTCA; site_90; 1; 22; 0; p
BS  TTATGAAGCCCTTCACAGAATT; site_91; 1; 22; 0; p
BS  GTGCGGGCGTGATCACAATTAC; site_92; 1; 22; 0; p
BS  CTGTGACTCGATTCACGAAGTC; site_93; 1; 22; 0; p
BS  ACGCGACTTTTATCACTTTTTA; site_94; 1; 22; 0; p
BS  CGGTGACGGAGTTCACCCTTTA; site_95; 1; 22; 0; p
BS  TTGCAATTCGTGTCACAAAATA; site_96; 1; 22; 0; p
BS  TTACGACAGCTATCACGAATTT; site_97; 1; 22; 0; p
BS  ATTTGAAGCAGTTAACGCTATT; site_98; 1; 22; 0; p
BS  ACATGAGCAACCGCACATATTT; site_99; 1; 22; 0; p
BS  CTGCGAGTGGGAGCACGGTTTT; site_100; 1; 22; 0; p
BS  TTTTGCGCTAAAGCACATTTCT; site_101; 1; 22; 0; p
BS  GTGCGAAATCCGTCACAGTTCA; site_102; 1; 22; 0; p
BS  TTGTGATGTGGTTAACCAATTT; site_103; 1; 22; 0; p
BS  TCGTGAACGATCCCACGAATTT; site_104; 1; 22; 0; p
BS  CTGCGAGCATGGTCATATTTTT; site_105; 1; 22; 0; p
BS  GAGTGATATGTATAACATTATG; site_106; 1; 22; 0; p
BS  ATATGACCAACCTCTCATAATT; site_107; 1; 22; 0; p
BS  ATGCAAAGGACGTCACATTACC; site_108; 1; 22; 0; p
BS  GAGTGATCGAGTTAACATTGTT; site_109; 1; 22; 0; p
BS  TGGTGAGGAACTTAACAATATT; site_110; 1; 22; 0; p
BS  TTGAGATTCAACTCTCAAATTT; site_111; 1; 22; 0; p
BS  TTGTGAATCAGATCAGAAAACC; site_112; 1; 22; 0; p
BS  ATGTATGACAGATCACTATTTT; site_113; 1; 22; 0; p
BS  ATATGATAAATATCAAACAATG; site_114; 1; 22; 0; p
BS  ATGTAAGCTGTGCCACGTTTTT; site_115; 1; 22; 0; p
BS  TTGTGAGTTTTGTCACCAAATA; site_116; 1; 22; 0; p
BS  TTGCGATGAATGTCACATCCTC; site_117; 1; 22; 0; p
BS  TTATGACGAGGCACACACATTT; site_118; 1; 22; 0; p
BS  TTGAGGGGTTGATCACGTTTTG; site_119; 1; 22; 0; p
BS  TGTTGCTTTTGATCACAATAAG; site_120; 1; 22; 0; p
BS  TCGTGACAGGAATCACGGAGTT; site_121; 1; 22; 0; p
BS  GCTTGAGCCGCAGCACAATGTG; site_122; 1; 22; 0; p
BS  GTTTGACGGCTATCACGTTTCA; site_123; 1; 22; 0; p
BS  ATTTAATTCGTATCGCAAATTA; site_124; 1; 22; 0; p
BS  TGGTGATCCATAAAACAATATT; site_125; 1; 22; 0; p
BS  CCGTGAAAGCGATCACAAAGGG; site_126; 1; 22; 0; p
BS  TCGTGATACTCATCACCATGAC; site_127; 1; 22; 0; p
BS  ATGTGATCTACAGCATGTTATG; site_128; 1; 22; 0; p
BS  GCGTGCCAGTTTTCACATTCTT; site_129; 1; 22; 0; p
BS  TCTTGCTTACCGTCACATTCTT; site_130; 1; 22; 0; p
BS  CAGTGACAGATTTCACGAAAAT; site_131; 1; 22; 0; p
BS  GGTTGCACTCTCTCACATTTTT; site_132; 1; 22; 0; p
BS  GCCTGACGGAGTTCACACTTGT; site_133; 1; 22; 0; p
BS  TCATGAAACTGTGCACATTTTA; site_134; 1; 22; 0; p
BS  GCAGGATTTAGCTCACACTTAT; site_135; 1; 22; 0; p
BS  ACGCAGCGAAGATCACAATTTA; site_136; 1; 22; 0; p
BS  CATTTAAACAGATCACAAAATC; site_137; 1; 22; 0; p
BS  TATTGATCTAACTCACGAAAAT; site_138; 1; 22; 0; p
BS  TTACAAGGCACATCACGTTATG; site_139; 1; 22; 0; p
BS  CTGTGACCGTGGTCGCAGTTGG; site_140; 1; 22; 0; p
BS  TTGGAATATCCATCACATAACG; site_141; 1; 22; 0; p
BS  ATATGCGCGAAATCAAACAATT; site_142; 1; 22; 0; p
BS  GTGCGAGTCTGCTCGCATAATC; site_143; 1; 22; 0; p
BS  ATGTGACAGATAAAACGTTTTA; site_144; 1; 22; 0; p
BS  AATTGACCGATGCCACGTTTTG; site_145; 1; 22; 0; p
BS  ATATGTTTCGTTTCACAGTTCT; site_146; 1; 22; 0; p
BS  ATGTGCGCATCTCCACATTACC; site_147; 1; 22; 0; p
BS  ATGTTAATTTCCTCACATCGTG; site_148; 1; 22; 0; p
BS  ATGCGGTGAGCATCACATCACC; site_149; 1; 22; 0; p
BS  ACGAGATCTGACACACACTATA; site_150; 1; 22; 0; p
BS  CAGCGACATCTGTCACATTCCT; site_151; 1; 22; 0; p
BS  TTGGGCGACAGATCACGCAAAA; site_152; 1; 22; 0; p
BS  TTTTGAATCGTGTCTCATTCTG; site_153; 1; 22; 0; p
BS  GTGAGGCATAAATCACATTACG; site_154; 1; 22; 0; p
BS  TCGTGATATTGCTCACGCCAAA; site_155; 1; 22; 0; p
BS  TTCTAATAGCCATCACAAAACG; site_156; 1; 22; 0; p
BS  TTGTGGATAAAATCACGGTCTG; site_157; 1; 22; 0; p
BS  CTCTGAGATGGATCAAAGAATT; site_158; 1; 22; 0; p
BS  CTTAGAAACCGATCACATACAG; site_159; 1; 22; 0; p
BS  AATCGATTGCGTTCACGTTTAC; site_160; 1; 22; 0; p
BS  TATTGATTTAAATCAAAGATTC; site_161; 1; 22; 0; p
BS  AATTGAACCAAATCATAAAATC; site_162; 1; 22; 0; p
BS  ATTGAACCCCGATCACACCATA; site_163; 1; 22; 0; p
BS  TAATGAAAAGGATGACATATTC; site_164; 1; 22; 0; p
BS  ACTTGCGTGACTACACATTCTT; site_165; 1; 22; 0; p
BS  AGCAGATACAACTCACACAATG; site_166; 1; 22; 0; p
BS  CTGTGCTGCGCATAATACTTTG; site_167; 1; 22; 0; p
BS  CCATGCTCAATCTCACAAAGTG; site_168; 1; 22; 0; p
BS  ATTCGACAAAGCGCACAATCCG; site_169; 1; 22; 0; p
BS  ACAGGAAGCACATCACAAAGAC; site_170; 1; 22; 0; p
BS  ATTTGCCACAGGTAACAAAAAA; site_171; 1; 22; 0; p
BS  TACGGATCTTCATCACATAAAA; site_172; 1; 22; 0; p
BS  CTGTGACAAGCTCCGCAAATCG; site_173; 1; 22; 0; p
BS  CAACGCTTTGGCTCACAGTTTA; site_174; 1; 22; 0; p
BS  ACAAAACATATGTCACAATATT; site_175; 1; 22; 0; p
BS  ATATGATCTATATCAATTTCTC; site_176; 1; 22; 0; p
BS  CTGTTGATATGATCACGTTATA; site_177; 1; 22; 0; p
BS  TTTAGATGTAAATCACTCCATT; site_178; 1; 22; 0; p
BS  ATGTGAAATAAAACAATTATTT; site_179; 1; 22; 0; p
BS  ACCTGTCACAAATCACAAAAAG; site_180; 1; 22; 0; p
BS  TAGTGCTCAGCGACACTATTTT; site_181; 1; 22; 0; p
BS  TTCTTATCTACCTCACAAAGGT; site_182; 1; 22; 0; p
BS  CTGTGCGCGCAACGACATTTTG; site_183; 1; 22; 0; p
BS  TTACAAAATTGTTAACAATTTT; site_184; 1; 22; 0; p
BS  ATTTTAACAACTTGACATATAT; site_185; 1; 22; 0; p
BS  CTTCGTAACGCCTCGCAAATTT; site_186; 1; 22; 0; p
BS  AACCTCTTTGCGTCACATTTTT; site_187; 1; 22; 0; p
BS  TTTTGTAACAATTCAAACTTCT; site_188; 1; 22; 0; p
BS  TCAGGCGTTAAATCACGTTTTC; site_189; 1; 22; 0; p
BS  TTCATCTCTATGTCACATTTTG; site_190; 1; 22; 0; p
BS  TTGTTACCTGCCTCTAACTTTG; site_191; 1; 22; 0; p
BS  ATGTGATATTGCTCTCCTATGG; site_192; 1; 22; 0; p
BS  ACGTTAACTGAAACGCATATTT; site_193; 1; 22; 0; p
BS  TATTGATGTAAATCAAATTCAC; site_194; 1; 22; 0; p
BS  CCATGAAACGGAACACGAAAAT; site_195; 1; 22; 0; p
BS  ATAATAATCTAATCACATCTTG; site_196; 1; 22; 0; p
BS  TGGTGCGCATGATAACGCCTTT; site_197; 1; 22; 0; p
BS  GTGCAGGCTTGATCACAACTCC; site_198; 1; 22; 0; p
BS  AGGTAATTCGAATGACATTGCT; site_199; 1; 22; 0; p
BS  CTGGGATGAAAGTGACATTTGA; site_200; 1; 22; 0; p
BS  CCGACGAATAGATCACAATTTA; site_201; 1; 22; 0; p
BS  ATTAGTAAGTTATCACCATTTG; site_202; 1; 22; 0; p
BS  CTGTGATAGTGTCATCATTTTC; site_203; 1; 22; 0; p
BS  GCTTAATTTAAATAACAAAATC; site_204; 1; 22; 0; p
BS  ATGTGATGAGAAAGTCAATTTG; site_205; 1; 22; 0; p
BS  TGACGCATGAAATCACGTTTCA; site_206; 1; 22; 0; p
BS  ATGCAAAATCAAAAACAATTTC; site_207; 1; 22; 0; p
BS  AATCGATTTAACACACCATTTA; site_208; 1; 22; 0; p
BS  TTGTAAACAGATTAACACCTCG; site_209; 1; 22; 0; p
BS  CAGTTATTTTTAACAAATTTTT; site_210; 1; 22; 0; p
BS  TGTTGCTGACCTTCAAAAATTA; site_211; 1; 22; 0; p
BS  TTGTCATCTTTCTGACACCTTA; site_212; 1; 22; 0; p
BS  AAGTTCGATATTTCTCGTTTTT; site_213; 1; 22; 0; p
BS  CTGGAACGCTTTTCGCATTCTG; site_214; 1; 22; 0; p
BS  ACCCGCTTTAAAACACGCTATC; site_215; 1; 22; 0; p
BS  GTTTGACATTGTTCTCTCACTT; site_216; 1; 22; 0; p
BS  CCTCGGTTTAGTTCACAGAAGC; site_217; 1; 22; 0; p
BS  TAGTTACATGTTTAACACTTGA; site_218; 1; 22; 0; p
BS  ATCTGGGTAGCATCACAGCAGA; site_219; 1; 22; 0; p
BS  AGGCGACCTGGGTCATGCTGAA; site_220; 1; 22; 0; p
BS  ATCTTGAAATAATCACATTGAT; site_221; 1; 22; 0; p
BS  ATTGTTATCCGCTCACAATTCC; site_222; 1; 22; 0; p
BS  TGTGGAAATTAATCCCACTATT; site_223; 1; 22; 0; p
BS  CTGCATATTAATTGACATTTCT; site_224; 1; 22; 0; p
BS  CGGGTATTAGCACCACATATAA; site_225; 1; 22; 0; p
BS  ATGTAATAAAATTCATGGTAAT; site_226; 1; 22; 0; p
BS  TTTTACTTTTGGTTACATATTT; site_227; 1; 22; 0; p
BS  TCTTTATCTTTGTAGCACTTTC; site_228; 1; 22; 0; p
BS  ACGCGCACTATGTCAACTCTTG; site_229; 1; 22; 0; p
BS  CAATGAAAAATTGCACAGTAAA; site_230; 1; 22; 0; p
BS  GTTAAATTGATGTAACATAATC; site_231; 1; 22; 0; p
BS  TTGAGTCATAAATAACCTTTAG; site_232; 1; 22; 0; p
BS  AGGTGAATCGCGCCAGCAAATT; site_233; 1; 22; 0; p
BS  CTGAGACTAGTACGACTTTTTG; site_234; 1; 22; 0; p
BS  AGTTGTAACTATGCACAAATGT; site_235; 1; 22; 0; p
BS  CCCCATGGCAGATGACATTTTT; site_236; 1; 22; 0; p
BS  AAAATATCCTTGTCACATTCGT; site_237; 1; 22; 0; p
BS  ATGTAAATTGGTCAACCATTGT; site_238; 1; 22; 0; p
BS  ACACTATGAGCAAAACACATTT; site_239; 1; 22; 0; p
BS  ATTTAATCATGTTTACAGTAAT; site_240; 1; 22; 0; p
BS  CTGGGAGCAGGCTCGATTTATG; site_241; 1; 22; 0; p
BS  ACATGATAAATTTGACGAAGAA; site_242; 1; 22; 0; p
BS  CAGTGAACCCCTTCCCAACCAC; site_243; 1; 22; 0; p
BS  GAATGCTAATCATCAGAAATGT; site_244; 1; 22; 0; p
BS  TTTGGGTTGTTATCAAATCGTT; site_245; 1; 22; 0; p
BS  TGTAAATCTGCATCGGAATTTG; site_246; 1; 22; 0; p
BS  TTCATACCACCATCACAACCAC; site_247; 1; 22; 0; p
BS  TATTTAACGCCGTCAGAAATGT; site_248; 1; 22; 0; p
BS  AAGGTGACTATACCACACTCAT; site_249; 1; 22; 0; p
BS  TTAGCAAATACCTCACAGTGAA; site_250; 1; 22; 0; p
BS  CTGCAACCATCTACAAATAACC; site_251; 1; 22; 0; p
BS  CGGGTAATAACAACACTCATAT; site_252; 1; 22; 0; p
BS  TTTCCTGAAAATTCACGCTGTA; site_253; 1; 22; 0; p
BS  ATTTAATGAAGAGAATTTTTTT; site_254; 1; 22; 0; p
BS  CTGCGTGAAGCAGCAGTAAATC; site_255; 1; 22; 0; p
BS  ACAGGCACCTGATAAAGCCATT; site_256; 1; 22; 0; p
BS  GGCAAGTGGTATTCGCACTTTT; site_257; 1; 22; 0; p
BS  GTACGTTTGCAGTGAAATAACT; site_258; 1; 22; 0; p
BS  ATAACAAAATCCTAATGTTATT; site_259; 1; 22; 0; p
CC  program: meme
CC  matrix.nb: 1
CC  id: CRP.2nt_upstream-noorf-ovlp-2str.22.meme_m1
CC  id: CRP.2nt_upstream-noorf-ovlp-2str.22.meme_m1
CC  ac: 
CC  name: CRP.2nt_upstream-noorf-ovlp-2str.22.meme_m1
CC  command: meme data/Sites_FNA_NR//CRP.21.fna -dna -mod oops -revcomp -text -nostatus -minw 22 -maxw 22 -bfile /state/partition1/space14/PGC/RegulonDB/local/data/Matrices/data/bg_freqs/2nt_upstream-noorf_Escherichia_coli_K12-ovlp-2str.meme_bg
CC  sites: 260
CC  meme.llr: 1864
CC  meme.E-value: 3.9e-382
CC  matrix.nb: 1
CC  sites: 260
CC  consensus.strict: atGtGatccagaTCACattttt
CC  consensus.strict.rc: AAAAATGTGATCTGGATCACAT
CC  consensus.IUPAC: wtGtGatcyrsaTCACahwttt
CC  consensus.IUPAC.rc: AAAWDTGTGATSYRGATCACAW
CC  consensus.regexp: [at]tGtGatc[ct][ag][cg]aTCACa[act][at]ttt
CC  consensus.regexp.rc: AAA[AT][AGT]TGTGAT[CG][CT][AG]GATCACA[AT]
XX
//
";

$demo_seq1=">sulA|2.75|1020250|I|-74|17.3
ccggctgtagtgttttccgtagagacacgcgcaattttacttgctgcggatgagaacgacgaagaacgatgtgcatagcctgaagtgtacataatcaatccagcccctgtgagttactgtatggatgtacagtacatccagtgacaacaaagatcaaccctattttcggaaagagcctcgcaaattttgtcgttggtgacgggaaaacataaattaatcttgccccttaagaataagttgcctattttcgtagttaacggatccgttaatgtgaatcattcttttatgttatgattttaaaaggaattttatgaaaagcctctcctataagcggatctataaatcacaagaatacctggcaacgttgggcacaattgaataccgatcattgtttggcagt
>ysdA|2.62|3851691|I|-358|n.d.
tgagataccgcgtgcacaacagctggcaacaggcagcggaaaggtacgtcagctggcagtgctcctgaaccacaggagacgcgtatgaacctggtggatatcgccattcttatcctcaaactcattgttgcagcactgcaactgcttgatgctgttctgaaatacctgaagtaattcagattcaagtcgcaccaaaggggagcgggaaaccgctccccttttatatttagcgtgcgggttggtgtcggatgcgatgctgacgcatcttatccgccctaccatctctcccggcaacatttattgccgcttttgtttacatattctgccgctaaacaattccccattcctggcgtatatctggctaacattcatcaatgtgatagattcctctcccgcattt
>sbmC|2.48|2079201|I|111|5.1
cgcatactgaccacctgtaatttctgtcagaatgacgccctcactgttttcgggaagcgtaaagtaacccggcaccgtcacgacggtgtcgcagcgtaatttttcggcgggtgtttcatctggattgtcgtaatagacagcaacccactccttcggcacaatatttttgctatctacccacatcatcaactgctcaaagcctttctttaccgtctgttcccacgggccaacgagatggaaacctgcaacggtacgtttctcttcctgcttaatctcgtagttcatgacgcctccattgatactgtttttatatacagtatagttgcaaaattaaaaccacaaggaatgagtgttgattatgcgagcagactcgcactcctgccagtctgctgcaaaagaa
>dinD|2.46|3815928|I|-200|10.6
taaatacagttacagatttactttctttgcaattgatatcacatggagtgggcaatgaacgaacatcatcaaccttttgaagagataaaactgattaatgcaaacggagcagaacaatggtcagcaagacaacttgggaaactactgggttattcagagtatcgtcactttatacctgtattaacgcgcgccaaagaagcctgtgagaacagtggtcacacaattgatgaccatttcgaggagatcctcgatatggtcaaaattggctcaaatgccaaaagagcattaaaagacatcgtactctcccgctatgcctgttacctggtagtacaaaacggcgaccctgcgaaaccggtcattgcggcagggcagacttattttgctatccagacccgacggc
>recN|2.38|2749644|I|113|31.4
ttctgaccccctctctggatgcgattaccctggtgcccatgttcccgcatacgttgtcagcacgaccactggtcataaacagcagcagcacgatccgtctgcgtttttcgcatcgccgtaacgacctggaaatcagttgcgacagccagatagcactgccgattcaggaaggtgaagatgtcctgattcgtcgctgtgattaccatctgaatctgattcatccgaaagattacagttatttcaacacattaagcaccaagctcggctggtcaaaaaaattattctaattttacgccagcctctttactgtatataaaaccagtttatactgtacacaataacagtaatggtttttcatacaggaaaacgactatgttggcacaactgaccatcagcaact
>lexA|2.26|4255589|I|-469|4.60
ggtgaaccacttctggcgcaacagcatattgaaggtcattatcaggtcgatccttccttattcaagccgaatgctgatttcctgctgcgcgtcagcgggatgtcgatgaaagatatcggcattatggatggtgacttgctggcagtgcataaaactcaggatgtacgtaacggtcaggtcgttgtcgcacgtattgatgacgaagttaccgttaagcgcctgaaaaaacagggcaataaagtcgaactgttgccagaaaatagcgagtttaaaccaattgtcgttgaccttcgtcagcagagcttcaccattgaagggctggcggttggggttattcgcaacggcgactggctgtaacatatctctgagaccgcgatgccgcctggcgtcgcggtttgtt
>dinI|1.93|1120728|I|13|3.1
tgtttatcttcttttgtcgcgccaataaccgataaattattcgctgcggcataacgtaccgatacgtggccttcattatcaggaaacgcatactgaatacggcgggaaagttcgccagccagggcgtcaatagccccagctggcaatggagaagttttcgctatggtgacttcaattcgcataatagccccctgttgaatatactggttatttatacaggtaaaataacctaatgacaacaggaagctacgatttttattgtttaacggaccagcgtaccgtttccccggcgaggaatggcaccagcgtgtcatcagtcagtgcgatgctttcagcaacctgttgctcttcacgtaccagttcgatgaatgtgtcgttgaccggcaacccatagaactgc
>minC|1.85|1225340|II|232|0.6
gtgcactgacgttgagtacaacgggggcatgttttaaaaatgcgggggcctgagcgattttgtcttccagcgcctgatggataaccttaggttctgcctcatgcagatgaaccacagataaagtgaagctactgcctttaagctcgattggcgtgtttgacatcctggccttactcaattagctattaatcatcgccagcgcgcgatgatgttccgaagactataaggcatgttatagtctggattatattgaggcaagtcaccctcccatttattcagagtaaaagtctattctgtgataaatggcgctgattcatagcttaaaaaatacccttgtcaatcaacccattgccgtcgtacttttgattgttcttatttacgcttcctttttccgcaccct
>yjiW|1.83|4577915|I|32|1.8
cggcttccagccactgacctttcagggtgatggcgggaatacggctgtaatccgggtagcgactcgcataaccgacggtgacatgacggttatttgccggggagacttctgcttcgaacggttgtgcaatagaatgcgtgtcagtcataactgctattctccaggaatagtgattgtgattagcgatgcgggtgtgttggcgcacatccgcaccgcgctaaatacctgtatatatcatcagtaaatatggggaaagtccagctaaaaatagaataaaatgggcaatttctggaatgatttaaatatatttatgtgggttatgattggcgtgaaataataaaaagcgcaccggaaaggtgcgccagaaaataatgttcaggattttttacgtgaggctttt
>umuD|1.81|1230172|I|-213|25.7
acttcaggcagattattatgttgtttatcaagcctgcggatctccgcgaaattgtgacttttccgctatttagcgatcttgttcagtgtggctttccttcaccggcagcagattacgttgaacagcgcatcgatctgaatcaactgttgatccagcatcccagcgcgacttacttcgtcaaagcaagtggtgattctatgattgatggtggaattagtgacggtgatttactgattgtcgatagcgctattaccgccagccatggtgatattgtcatcgctgctgttgacggcgagtttacggtgaaaaaattgcaactacgcccgacggtacagcttattcccatgaacagcgcgtactcgcccattaccatcagtagtgaagatacgctggatgtctt
>recA|1.75|2821814|I|45|10.0
catcggcagaccacctgccccaagcgcgatatccagtgaaagcgaaccggtagagatggtttccacatccatggaacggtcttcacccaggcgcatgatggagcctttaccaaattgtttctcaatctggcccagtgctgccgccaacgctttctgtttgttttcgtcgatagccatttttactcctgtcatgccgggtaataccggatagtcaatatgttctgttgaagcaattatactgtatgctcatacagtatcaagtgttttgtagaaattgttgccacaaggtctgcaatgcatacgcagtagcctgacgacgcaccgcatcacggtcgccgctgaagcattcccgccgggtaatgccttcaccgcgggcagtggcaaaagcaaaccagacggt
>uvrA,ssb|1.50|4272039|I|-50|1.9,1.2
ggttgatgtttttgagattatgggtgcgggcgccccgaacttcgatcttatccattcacctttcccggattaaacgcttttttgcccggtggcatggtgctaccggcgatcacaaacggttaattatgacacaaattgacctgaatgaatatacagtattggaatgcattacccggagtgttgtgtaacaatgtctggccaggtttgtttcccggaaccgaggtcacaacatagtaaaagcgctattggtaatggtacaatcgcgcgtttacacttattcagaacgatttttttcaggagacacgaacatggccagcagaggcgtaaacaaggttattctcgttggtaatctgggtcaggacccggaagtacgctacatgccaaatggtggcgcagttgc
>otsB,otsA|1.42|1979908|III|NA|0.9,1.4
ccggcacacccgccagtcgccatgatgcctgagttgcacctgtgccaatttttactgacattccgcccagtcggttaacgactgcgaagccagattcatcggttaaatcatcgcccagaaatacgggcgttcgcccgataaagggagcttcctgcataaaagctgcaattgcctcacctttactggtacctctcggtttgatctcgacaacacactttccctgctgtaacgccatttgtggccagatctgagtaatacgttgcgctaatgtcattaatgcgtcttcatgctgcggagcctgacgataatgcagcgcaaaagccatccctttcgcctccagctccgcgccgggatactgagcgatgactgtatgcagttgcacgctaatatcacgcgcaat
>yfaX,yfaW|1.41|2359548|III|NA|n.d,n.d
cgtgataatcaccgccgccagcgcctttttctgctgtcgcaccgccagtaaaccaggcgcgaacctgtttaatttttggtagggtcatgatgttctccattgttatgaggcttgtaagtcaaagggacttttccatcccaacagacgtgaaatatccctggcgcaggcaatggccttgcccgccagataatcacggtattcttcattgatttgtaagcgggtaccgaccaccgagatcgcagcggtaagctcgttattggcgttaaacaccggcgcagcgacacaacggacatcggcgtaatcttcgccgttgtcatagctccagccctgacggcgaatacgcgccagttcttcgtgaagttgctgtggatgagtaatcgttgtgggtgtcgcctgctcc
>polA|1.39|4047394|III|NA|0.8
caacggcggcagaagtgtttggtttgccactggaaaccgtcaccagcgagcaacgccgtagcgcgaaagcgatcaactttggtctgatttatggcatgagtgctttcggtctggcgcggcaattgaacattccacgtaaagaagcgcagaagtacatggacctttacttcgaacgctaccctggcgtgctggagtatatggaacgcacccgtgctcaggcgaaagagcagggctacgttgaaacgctggacggacgccgtctgtatctgccggatatcaaatccagcaatggtgctcgtcgtgcagcggctgaacgtgcagccattaacgcgccaatgcagggaaccgccgccgacattatcaaacgggcgatgattgccgttgatgcgtggttacaggc
>uvrB|1.37|812662|I|NA|3.0
agataaatgcaatggcagtcactgaacaggcatctcttgccataaaactgtcatcactcatcttgacaaatgttaaaaaagccgttgctttggggataacccggtaaggccggagttttatctcgccacagagtaaattttgctcatgattgacagcggagtttacgctgtatcagaaatattatggtgatgaactgtttttttatccagtataatttgttggcataattaagtacgacgagtaaaattacatacctgcccgcccaactccttcaggtagcgactcatgagtaaaccgttcaaactgaattccgcttttaaaccttctggcgatcagccagaggcgattcgacgtctcgaagaggggctggaagatggcctggcgcaccagacgttactt
>molR|1.33|2194300|I|177|1.4
ggttctgacaaactgctgcgcctgacgctggatctcggcggtgaaaaacgcaatgtcttctccggtattcgttctgcttacccggatccgcaggcactgattggtcgtcacaccattatggtggctaacctggcaccacgtaaaatgcgcttcggtatctctgaaggcatggtgatggctgccggtcctggcgggaaagatattttcctgctaagcccggatgccggtgctaaaccgggtcatcaggtgaaataatcccccttcaaggcgctgcatcgacagcgccttttctttataaattcctaaagttgttttcttgcgattttgtctctctctaacccgcataaatactggtagcatctgcattcaactggataaaattacagggatgcagaatgag
>ydjF|1.30|1852700|III|NA|n.d.
ttccgccgagttggttagcaacgtcaggccactacggtcctgtaacaatttgagcaattccattacggtactactggaatcggctgccatggtggttttattgtcgataaagggtagtgccttgcgtgcaataagctgcttctcttcataaaacgatgaagcgcgcttataaaaatggatattctccgtcaacatcgctgtatttaaaacagcaccaccataggttctggtcaaaaagccttcatcttccagcttctcaagatcgcggcgaatggtttcttcggttacctgaaaaatcccactcaaatttgagactgtcacctttttatcgttggcaaccatttgcttaattgcctgaatcctgtcttttgccgccacgattacacccctgtatcttttt
>dinQ|1.29|3646206|I|-192|n.d.
ggataatcatacagtacatgcaggttataaaaccagcacgtccttgcaatagtttcagtatggtattagcattgatgcgttagatgatggctatctcactccagtcagagccaccaactcagggctggaaagtaaaaaaccgacgcaaagtcggtttttttacatccggattcggacaaggcttaatatgacgatgacccagtgaaagtatataaatcgtcactgcgatatataccgaagtgctccctccgccagctgaagaaatcgctaattcttgcaatgttagccactggctaatagtattgagctgttagataagaactctctcactccagccagagccaccaactcagggctggaaagtaaaaaaccgacgcaaagtcggtttttttacgtcctg
";
$demo_seq2=">CRP|34166-34187
ttgatgacatAAGCAGGATTTAGCTCACACTTatcgacggtg
>CRP|42068-42089
ataagctgtaTTCTGTGATTGGTATCACATTTttgtttcggg
>CRP|989834-989855
aaaagatgctAAAGGTTATTTATGACTCAACAgccacaagcc
>CRP|1019435-1019456
acgtagttgaAAACTTACAAGTGTGAACTCCGtcaggcatat
>CRP|1078316-1078337
aggtgcaaccGCAAAAAATGTGAGAGAGTGCAacctgatgaa
>CRP|1078346-1078367
caacctgatgAAAAATAGTGTCGCTGAGCACTaaaatttaat
>CRP|1078378-1078399
aaaatttaatGTAAATGGTGTGTTAAATCGATtgtgaataac
>CRP|1102599-1102620
atcttaaatcAAGTGTTAAACATGTAACTAAAtgtaactcgt
>CRP|1102727-1102748
aaatcatacaAATGGTGATAACTTACTAATAAtgcatataaa
>CRP|1156884-1156905
acgtttgaaaTATTGTGACATATGTTTTGTCAaaatgtgcaa
>CRP|1156938-1156959
tctgaagttgAAACGTGATAGCCGTCAAACAAattggcactg
>CRP|1229740-1229761
ttatttctatTAATATGATAAATATCAAACAAtgtttaatgt
>CRP|1236677-1236698
cagagtcaggGAGATGTGAGCCAGCTCACCATaaaaaagccg
>CRP|1236742-1236763
gagttatcaaGATGTGATTAGATTATTATTCTtttactgtat
>CRP|1327314-1327335
atgcctgttgTAAACTGTGAGCCAAAGCGTTGtttaaccaag
>CRP|1333754-1333775
ctcttttatcAATTTGGGTTGTTATCAAATCGttacgcgatg
>CRP|1445434-1445455
tattgtcacgATTTGCGGAGCTTGTCACAGCTgacaaagcga
>CRP|1451820-1451841
tataaaaataGGGTGCGAAATCCGTCACAGTTcaaacataca
>CRP|1486144-1486165
ttacacttgtTTTTATGAAGCCCTTCACAGAAttgtcctttc
>CRP|1570148-1570169
aaggactcgtGTTTAAATAACAAAATCCTAATgttatttatc
>CRP|1617054-1617075
cgatttagcaAAACGTGGCATCGGTCAATTCAttcatttgac
>CRP|1666675-1666696
ttttacttagGCATGTGATTAACAGCACATTTttcgggcttt
>CRP|1686856-1686877
ttcctgttcaAAGTATTATGCGCAGCACAGCCactctccatt
>CRP|1694255-1694276
gcaatacccaCAGCGTGATATAGATCGCATTAatctttaaaa
>CRP|1697225-1697246
tactccctgaTTATGTGACAGATAAAACGTTTtaccttttat
>CRP|1697283-1697304
attatcgttgCGTAATGTGATTTATGCCTCACtaaaatttga
>CRP|1819823-1819844
cttgacctgtGGTTATGACCCTCTTCACATTTcgggcaaata
>CRP|1840289-1840310
gggtcattttTTTCTTGCTTACCGTCACATTCttgatggtat
>CRP|1860569-1860590
gctgcacctaAATCGTGATGAAAATCACATTTttatcgtaat
>CRP|1887880-1887901
ccaacgaaaaGGTTGCGAAGCGCGTCACTATTtatttttatc
>CRP|1899854-1899875
ctttgcaaacGAATGTGACAAGGATATTTTACctttcgaaat
>CRP|1899906-1899927
cgaaagttaaATTACGGATCTTCATCACATAAaataattttt
>CRP|1976480-1976501
attttcaataATGCGTGATGCAGATCACACAAaacactcaat
>CRP|1984293-1984314
attaattctcCATAGGAGAGCAATATCACATCgcagaattac
>CRP|2175286-2175307
tgtttttaaaTATCGAGATAACGATCACAAAAacgacaatat
>CRP|2229737-2229758
aatgggctaaAATTTGCGATGCGTCGCGCATTtttgatgtat
>CRP|2229789-2229810
ttgcataattAATGAGATTCAGATCACATATAaagccacaac
>CRP|2238619-2238640
attgttaagaTACTGTGAAATCACTCACAGATtgaaagcggt
>CRP|2239764-2239785
gaatacaggaCTTCGTGAATCGAGTCACAGCAatggaaacgg
>CRP|2342630-2342651
gaaatcgccgAACAGTTATTTTTAACAAATTTttctcttccc
>CRP|2350504-2350525
ttcataaattAAATGTGAATTGCCGCACACATtattaaataa
>CRP|2350554-2350575
aaatgttcaaAATGACGCATGAAATCACGTTTcactttcgaa
>CRP|2350594-2350615
aattatgagcGAATATGCGCGAAATCAAACAAttcatgtttt
>CRP|2459272-2459293
cacttcgcgcTCCTGTTACAGCACGTAACATAgtttgtataa
>CRP|2475715-2475736
tgagatgagcTAAAGTGAACCATATCTCAATTcaccttcatt
>CRP|2475745-2475766
ttcaccttcaTTTTTAGATGTAAATCACTCCAttgatgcaat
>CRP|2510931-2510952
aacatagcagAAATGTATGACAGATCACTATTtttgaagcct
>CRP|2510980-2511001
gacgtcattaTAGTGTGTGTCAGATCTCGTTTtccttaacca
>CRP|2531450-2531471
taggtgctttTTTGTGGCCTGCTTCAAACTTTcgcccctcct
>CRP|2531577-2531598
gctgaatcgaTTTTATGATTTGGTTCAATTCTtcctttagcg
>CRP|2599022-2599043
atacctcactTCTCGTGATCAAGATCACATTCtcgctttccc
>CRP|2632236-2632257
tgctgatttaGAATTTGATCTCGCTCACATGTtaccttctca
>CRP|2663370-2663391
acgcttgccaACATTTCTGATGATTAGCATTCccttcgccat
>CRP|2714575-2714596
acacccttgaATCTTTGATTTAAATCAATAAAaaccacacat
>CRP|2749002-2749023
tgcaactgaaGAATGTGAAAACTGGCACGCTCgcggagattg
>CRP|2786870-2786891
atgttactaaTTTGTTGCTTTTGATCACAATAagaaaacaat
>CRP|2823460-2823481
atctgcatcaCATTGTGCTGCGGCTCAAGCAGgaagccgcca
>CRP|2823761-2823782
ttatctttcaTTTTGCGATCAAAATAACACTTttaaatcttt
>CRP|2866075-2866096
cgttacaggcGCTGTGACCGTGGTCGCAGTTGgcttgttgtt
>CRP|2866192-2866213
ggtggaatttTGTGCAGGCTTGATCACAACTCcttgctctgc
>CRP|2932050-2932071
tttcgattatTAAAGTGATGGTAGTCACATAAagtcaccttc
>CRP|2932090-2932111
tctagctaatAAGTGTGACCGCCGTCATATTAcagagcgttt
>CRP|2980314-2980335
cgacatgtcgTTATGTGATGGATATTCCAATTttcaaattaa
>CRP|3056640-3056661
acatgtcaccAAATTTAATGAAGAGAATTTTTttaacggggg
>CRP|3071905-3071926
aagaagacatTTATCTGACTCACATCACACTTttatcccctt
>CRP|3084528-3084549
gtttttaacaGAATGAGACACGATTCAAAAAAaagtggaaat
>CRP|3086224-3086245
gcgattacacTGATGTGATTTGCTTCACATCTttttacgtcg
>CRP|3098854-3098875
gagatctacaAAGTTAGAGGCAGGTAACAAAAcgaagaatta
>CRP|3103570-3103591
tcaggggcaaAAATGTTATCCACATCACAATTtcgttttgca
>CRP|3103620-3103641
atgtttgcaaTTATTTGCCACAGGTAACAAAAaaccagtccg
>CRP|3103718-3103739
atccgcatcaCGATGTGAGGAAATTAACATGAatcttaagct
>CRP|3126199-3126220
ggttaactcaATGTTAAATTGATGTAACATAAtcacttacgt
>CRP|3190064-3190085
cgtgcgtaaaTATATTGTCCCCGATCACACTTtttagtgaaa
>CRP|3265148-3265169
aatgaacaggATATGTGCGACCACTCACAAATtaactttcaa
>CRP|3276817-3276838
attatgtttcTTTTGTGAATCAGATCAGAAAAccattatctt
>CRP|3316406-3316427
ctgcaggcatTATAGTGATCCACGCCACATTTtgtcaacgtt
>CRP|3371692-3371713
agaaagcttaTAATGCGATCTGCTTCACTAAAgtggcattat
>CRP|3382385-3382406
ccacatctcaAGAATGTGTAGTCACGCAAGTTtagcgtttat
>CRP|3408005-3408026
acgcggggcgAAGTGCGAGCAAGCTCACAAAAggcacgtaaa
>CRP|3408090-3408111
cataaagaaaAATTGAGAACTTACTCAAATTTctttgagtgt
>CRP|3408194-3408215
tttaaatgcaATTCTTTGATCCATCTCAGAGGattggtcaaa
>CRP|3483905-3483926
cgcgcaacggAAGGCGACCTGGGTCATGCTGAagcgagacac
>CRP|3484005-3484026
atgtactgcaTGTATGCAAAGGACGTCACATTaccgtgcagt
>CRP|3490449-3490470
tcactttttaTTCCGTGATCAAAATCACCTCTtaaaatgcaa
>CRP|3490502-3490523
attgcaataaAACATTTAAACAGATCACAAAAtcacctaaaa
>CRP|3491702-3491723
acgctgtcgtCTTTGTGATGTGCTTCCTGTTAggtttcgtca
>CRP|3491957-3491978
atataaaggtGAATTTGATTTACATCAATAAGcggggttgct
>CRP|3534771-3534792
ttattgattcTGTTGATATGATCACGTTATACccaatgtgcg
>CRP|3542637-3542658
gcctaaccaaATGGCTTTATCAGGTGCCTGTTgcagcacggc
>CRP|3544345-3544366
ttccccgcgaTAATATGACCAACCTCTCATAAtttaaattta
>CRP|3550615-3550636
atgagcaaggAAATGTGATCTCAACCACTTAAagctagtgca
>CRP|3550965-3550986
aaagatttggAATTGTGACACAGTGCAAATTCagacacataa
>CRP|3559920-3559941
ataaacgccaTAATGTTATACATATCACTCTAaaatgttttt
>CRP|3575726-3575747
ttaagtgtatAAGTGTGAGCTACTTCAAATTTgtgggcttaa
>CRP|3590465-3590486
taagatagtaAGGTGTCAGAAAGATGACAAGGcggtgacggc
>CRP|3598903-3598924
tgttttatcaGACCGTGATTTTATCCACAAGTtcaatgcaag
>CRP|3740582-3740603
ggtagcgcaaAGTGTGCCGTAGTTCACGATCTcgacagataa
>CRP|3740748-3740769
aaaacaaggaAGCCTGGGATGAAAGTGACATTtgagcagtta
>CRP|3754608-3754629
cggtaacagcTTTACGACAGCTATCACGAATTtacgggcaag
>CRP|3769939-3769960
acgaaggcatAACATGCTGTAGATCACATCAGgtgaacgccg
>CRP|3769981-3770002
taagaaaataTCTTGTGATTCAGATCACAAAGattcaacaaa
>CRP|3770025-3770046
atcaaaacaaAAATGTGACACTACTCACATTTaaatgccatt
>CRP|3770098-3770119
tacacaagcgTTTTGTGATGAACGTCACGTCAattacctctc
>CRP|3770142-3770163
ccccctatatTTATGTGATTGATATCACACAAaaggccgtcg
>CRP|3845314-3845335
tatgaagtgaAAAGGTGAGATGCATCACGCTTcgcgcggtgt
>CRP|3851099-3851120
agcaagggaaAATTGAGGGGTTGATCACGTTTtgtactgaat
>CRP|3886364-3886385
ctccccgaacGATTGTGATTCGATTCACATTTaaacaatttc
>CRP|3904771-3904792
ctattgataaAAATATGACCATGCTCGCAGTTattaactttg
>CRP|3931272-3931293
ccatgtaaaaCGTTTCGAGGTTGATCACATTTccgtaacgtc
>CRP|3963636-3963657
tctgttctgtTAAATGTGTTTTGCTCATAGTGtggtagaata
>CRP|3989009-3989030
cagcaaggtgTTAAATTGATCACGTTTTAGACcattttttcg
>CRP|4014310-4014331
gaagggtgatTTATGTGATTTGCATCACTTTTggtgggtaaa
>CRP|4014362-4014383
gcatttgcgtCATGGTGATGAGTATCACGAAAaaatgttaaa
>CRP|4047827-4047848
agctgtcacgTTTTGTGATGGCTATTAGAAATtcctatgcaa
>CRP|4056175-4056196
attgcccctaAAAGGCGTTATCATGCGCACCAtcgtgcaaaa
>CRP|4056306-4056327
ctttttatgcTCCGTGAAAGCGATCACAAAGGgactctgcaa
>CRP|4095577-4095598
agggaaagatGAACGTGATGATGTTCACAATTtgctgaattg
>CRP|4095607-4095628
tttgctgaatTGTGGTGATGTGATGCTCACCGcatttcctga
>CRP|4095631-4095652
gctcaccgcaTTTCCTGAAAATTCACGCTGTAtcttgaaaaa
>CRP|4098672-4098693
acggcattaaGTGGGTGATTTGCTTCACATCTcgggcatttt
>CRP|4116211-4116232
taacgagcaaAAACGAGAAATATCGAACTTAAaatgtgtgtg
>CRP|4116234-4116255
cgaacttaaaATGTGTGTGCCTCGTCATAAAAtgagcgttat
>CRP|4122586-4122607
ttttcatgaaAAGTGTGATGAATATTGAATTTttcgatccgc
>CRP|4198125-4198146
taagaccagaAAACGTGATTTAACGCCTGATTtgtcgtacct
>CRP|4199694-4199715
ccagggcaatTTTCGTGTTCCGTTTCATGGTTaatcctccag
>CRP|4213384-4213405
ttaaaatggaAATTGTTTTTGATTTTGCATTTtaaatgagta
>CRP|4238287-4238308
acaaaaaataTAGATCTCCGTCACATTTTTGCgttatacagg
>CRP|4244553-4244574
caccgtcgctTTGTGTGATCTCTGTTACAGAAttggcggtaa
>CRP|4244582-4244603
gaattggcggTAATGTGGAGATGCGCACATAAaatcgccacg
>CRP|4244616-4244637
tcgccacgatTTTTGCAAGCAACATCACGAAAttccttacat
>CRP|4244648-4244669
ttccttacatGACCTCGGTTTAGTTCACAGAAgccgtgttct
>CRP|4285473-4285494
tttgtttagtATTTGGGCGACAGATCACGCAAaagtagaatt
>CRP|4285526-4285547
cggcagggtaATTTTTGAAGGTCAGCAACAAAagttgattaa
>CRP|4328297-4328318
aaacgaaatcCATGTGTGAAGTTGATCACAAAtttaaacact
>CRP|4339704-4339725
gagggaggatGACTGCGAGTGGGAGCACGGTTttcaccctct
>CRP|4339818-4339839
taaactcagaTTTACTGCTGCTTCACGCAGGAtctgagttta
>CRP|4346883-4346904
tcgttaccggCTTTAGCAAATACCTCACAGTGaatattggct
>CRP|4366535-4366556
ttaattattaATTTGTGAAATAGATCACCGCTttgggattac
>CRP|4417879-4417900
atgtggaattATTTGCGGGTCGCGTCACATTTaatcataaat
>CRP|4434685-4434706
tacgctttgaAAATGATGACACTATCACAGTTggcgcattca
>CRP|4453703-4453724
caatgttgcgCTCAGGTGAATCGCGCCAGCAAattacggatt
>CRP|4464211-4464232
cataaagcccCATGGCAGATGACATTTTTGGTtggctgcaga
>CRP|4464283-4464304
aacgttcccgAAACGCAGCGAAGATCACAATTtatcgttcag
>CRP|4492489-4492510
aactgtcattATTTGTGATGAAGATCACGTCAgaaaattgtt
>CRP|4492539-4492560
atgttacgcaTAACGTGATGTGCCTTGTAATTcttatcagta
>CRP|4538002-4538023
acaaaaccagATTTGCAATTCGTGTCACAAAAtatgtcgatc
>CRP|4549367-4549388
gggcgatatgTTATGTAAATTGGTCAACCATTgttgcgatga
>CRP|4549388-4549409
ggtcaaccatTGTTGCGATGAATGTCACATCCtctgatcaat
>CRP|4549473-4549494
agagtgaaatTCTTGTGATGTGGTTAACCAATttcagaattc
>CRP|4609153-4609174
gcggtttcaaAATTGTGATCTATATTTAACAAagtgatgaca
>CRP|4609177-4609198
tttaacaaagTGATGACATTTCTGACGGCGTTaaataccgtt
>CRP|4615197-4615218
tgaaagtgaaTTATTTGAACCAGATCGCATTAcagtgatgca
>CRP|4615250-4615271
agatttccttAATTGTGATGTGTATCGAAGTGtgttgcggag
";

$demo_descr="<H4>Comment on the demonstration example : </H4><blockquote class ='demo'>In this demonstration, we will analyse the PSSM of the Transcription Factor LexA available in RegulonDB. </p>
As a positive set we will use the obtained sequences from the ChIP-chip experiment (Wade et al. Genes Dev. 2005) of transcription factor LexA in the Escherichia coli K12 genome.</p>
As a negative sequence set we will use the reported binding sites of CRP in the Escherichia coli K12 Genome annotated at RegulonDB. </p>";
$demo_markov=1;

print "<td><b>";
print $query->hidden(-name=>'demo_descr',-default=>$demo_descr."</blockquote>");
print $query->hidden(-name=>'html_title',-default=>$demo_html_title);
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'kfold',-default=>$default{kfold});

print $query->hidden(-name=>'tag1',-default=>'positive_set');
print $query->hidden(-name=>'sequence1',-default=>$demo_seq1);
print $query->hidden(-name=>'permutation1',-default=>1);
print $query->hidden(-name=>'scanopt1',-default=>'');

print $query->hidden(-name=>'tag2',-default=>'negative_set');
print $query->hidden(-name=>'sequence2',-default=>$demo_seq2);
print $query->hidden(-name=>'markov_order',-default=>$demo_markov);
print $query->hidden(-name=>'scanopt2',-default=>'');
print $query->hidden(-name=>'nwd',-default=>'CHECKED');

print $query->submit(-label=>"DEMO");
print "</b></td>\n";
print $query->end_form;


print "<td><b><a href='help.matrix-quality.html'>MANUAL</A></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

