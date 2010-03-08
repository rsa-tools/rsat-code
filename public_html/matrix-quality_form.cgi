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
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
use RSAT::matrix;
use RSAT::MatrixReader;

### Read the CGI query
$query = new CGI;

local @supported_input_formats = sort(keys( %RSAT::MatrixReader::supported_input_format));
local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

################################################################
### default values for filling the form
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "consensus";
$default{permutation1} = "1";
$default{permutation2} = "1";
$default{tag1} = "sequence_set1";
$default{tag2} = "sequence_set2";
$default{pseudo_prior} = "pseudo_prior";
$default{pseudo_counts}="1";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{bg_method}="bgfile";
$checked{$default{bg_method}} = "CHECKED";
$default{organism}="Escherichia_coli_K12";

&ReadMatrixFromFile();

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 


################################################################
### print the form ###


################################################################
### header
&RSA_header("matrix-quality", "form");
print "<CENTER>";
print "description position-specific scoring matrices (PSSM), and calculate statistical parameters.<P>\n";
print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"matrix-quality.cgi");

#print "<FONT FACE='Helvetica'>";

################################################################
#### Matrix specification
print "<hr>";
&GetMatrix();
print "<p><font color=red>Only the first matrix will be taken in acount</font></p>";
print "<hr>";

################################################################
#### Sequence specification

print "<h2> Mandatory Sequence </h2>";

&SeqBoxMQ(1);
print "<hr>";

print "<h2> Optional Sequence </h2>";
&SeqBoxMQ(2);

print "<hr>";

################################################################
#### Background specifiaction
my %bg_params =(
    "markov" => 1,
    "markov_message" => 1
    );
&GetBackgroundModel(\%bg_params);

print "<br/>Note: Only Bernoulli models are supported. Higher-order Markov models are converted into Markov 0 (Bernoulli).";
print "<hr>";



################################################################
### send results by email or display on the browser
print "<p>\n";
&SelectOutput("email", email_only=>1);



################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"matrix-quality_form.cgi");
$demo_matrix="PRIOR FREQUENCIES DETERMINED BY OBSERVED FREQUENCIES.
letter   1: A  prior frequency = 0.296739
letter   2: T  prior frequency = 0.327174
letter   3: C  prior frequency = 0.202174
letter   4: G  prior frequency = 0.173913

THE LIST OF MATRICES FROM FINAL CYCLE

MATRIX 1
number of sequences = 23
unadjusted information = 11.7434
sample size adjusted information = 10.3038
ln(p-value) = -177.341   p-value = 9.58579E-78
ln(expected frequency) = -108.44   expected frequency = 8.04114E-48
A	|  12   0   0   0   1  12   1  12   6  10   7  13   4  12   0  23   0   1  12   6  11
T	|   3   1  23   0  14   5  15   5  15   6  11   5  12   2   0   0   0  13   6  13   8
C	|   3  22   0   0   2   3   5   2   2   5   5   2   4   7  23   0   0   8   2   2   3
G	|   5   0   0  23   6   3   2   4   0   2   0   3   3   2   0   0  23   1   3   2   1
   1|1   :    1/1     ACTGTATAAAACCACAGCCAA
   2|2   :    2/1     GCTGCGCTTATCGACAGTTAT
   3|3   :    3/1     CCTGGCTTTCAGGGCAGCGTT
   4|4   :    4/1     ACTGTTTTTTTATCCAGTATA
   5|5   :    5/1     ATTGGCTGTTTATACAGTATT
   6|6   :    6/1     CCTGTTAATCCATACAGCAAC
   7|7   :    7/1     ACTGTACATCCATACAGTAAC
   8|8   :    8/1     TCTGCTGGCAAGAACAGACTA
   9|9   :    9/1     ACTGTATATAAAAACAGTATA
  10|10  :   10/1     GCTGGATATCTATCCAGCATT
  11|11  :   11/1     GCTGGATATCTATCCAGCATT
  12|12  :   12/1     ACTGTGCCATTTTTCAGTTCA
  13|13  :   13/1     ACTGTGCCATTTTTCAGTTCA
  14|14  :   14/1     ACTGTATATAAAACCAGTTTA
  15|15  :   15/1     ACTGTACACAATAACAGTAAT
  16|16  :   16/1     ACTGTATGAGCATACAGTATA
  17|17  :   17/1     GCTGGCGTTGATGCCAGCGGC
  18|18  :   18/1     ACTGTTTATTTATACAGTAAA
  19|19  :   19/1     TCTGTATATATACCCAGCTTT
  20|20  :   20/1     TCTGGTTTATTGTGCAGTTTA
  21|21  :   21/1     GCTGTATATACTCACAGCATA
  22|22  :   22/1     ACTGTATATACACCCAGGGGG
  23|23  :   23/1     CCTGAATGAATATACAGTATT
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


print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'consensus');

print $query->hidden(-name=>'tag1',-default=>'positive_set');
print $query->hidden(-name=>'sequence1',-default=>$demo_seq1);
print $query->hidden(-name=>'permutation1',-default=>1);

print $query->hidden(-name=>'tag2',-default=>'negative_set');
print $query->hidden(-name=>'sequence2',-default=>$demo_seq2);

print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.convert-matrix.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_PSSM.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

