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
$default{demo_descr}="";
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "meme";
$default{kfold}="0";
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
$default{organism}="Escherichia_coli_K_12_substr__MG1655_uid57779";
#$default{html_title}="";
$default{markov_order} = "0";
$default{nwd}="";


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
 print "<\p><b>K fold validation</B>&nbsp;";
 print $query->popup_menu(-name=>'kfold',
			   -Values=>[0,3,4,5,6,7,8,9,10],
			   -default=>$default{kfold});
print "<p><font color=red>Only the first matrix will be taken in acount</font></p>";

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
$demo_html_title=" LexA matrix from RegulonDB";
$demo_matrix="********************************************************************************
MEME - Motif discovery tool
********************************************************************************
MEME version 3.5.4 (Release date: 3.5.4)

For further information on how to interpret these results or to get
a copy of the MEME software please access http://meme.nbcr.net.

This file may be used as input to the MAST algorithm for searching
sequence databases for matches to groups of motifs.  MAST is available
for interactive use and downloading at http://meme.nbcr.net.
********************************************************************************


********************************************************************************
REFERENCE
********************************************************************************
If you use this program in your research, please cite:

Timothy L. Bailey and Charles Elkan,
\"Fitting a mixture model by expectation maximization to discover
motifs in biopolymers\", Proceedings of the Second International
Conference on Intelligent Systems for Molecular Biology, pp. 28-36,
AAAI Press, Menlo Park, California, 1994.
********************************************************************************


********************************************************************************
TRAINING SET
********************************************************************************
DATAFILE= /bio/jvanheld2/matrix_eval/data/Sites_FNA_NR/LexA.fna
ALPHABET= ACGT
Sequence name            Weight Length  Sequence name            Weight Length  
-------------            ------ ------  -------------            ------ ------  
LexA|65834-65853         1.0000     40  LexA|738568-738587       1.0000     40  
LexA|738659-738678       1.0000     40  LexA|812655-812674       1.0000     40  
LexA|832259-832278       1.0000     40  LexA|932351-932370       1.0000     40  
LexA|1020162-1020181     1.0000     40  LexA|1229931-1229950     1.0000     40  
LexA|1229951-1229970     1.0000     40  LexA|1944050-1944069     1.0000     40  
LexA|1944102-1944121     1.0000     40  LexA|2749749-2749768     1.0000     40  
LexA|2749771-2749790     1.0000     40  LexA|2821851-2821870     1.0000     40  
LexA|3208751-3208770     1.0000     40  LexA|3851322-3851341     1.0000     40  
LexA|3995930-3995949     1.0000     40  LexA|4255050-4255069     1.0000     40  
LexA|4255091-4255110     1.0000     40  LexA|4255112-4255131     1.0000     40  
LexA|4271978-4271997     1.0000     40  
********************************************************************************

********************************************************************************
COMMAND LINE SUMMARY
********************************************************************************
This information can also be useful in the event you wish to report a
problem with the MEME software.

command: meme /bio/jvanheld2/matrix_eval/data/Sites_FNA_NR/LexA.fna -dna -mod oops -revcomp -nostatus -minw 22 -maxw 22 -bfile /bio/jvanheld2/matrix_eval/data/bg_freqs/2nt_upstream-noorf_Escherichia_coli_K12-ovlp-2str.meme_bg -dir /home/scmbb/installations 

model:  mod=          oops    nmotifs=         1    evt=           inf
object function=  E-value of product of p-values
width:  minw=           22    maxw=           22    minic=        0.00
width:  wg=             11    ws=              1    endgaps=       yes
nsites: minsites=       21    maxsites=       21    wnsites=       0.8
theta:  prob=            1    spmap=         uni    spfuzz=        0.5
em:     prior=   dirichlet    b=            0.01    maxiter=        50
        distance=    1e-05
data:   n=             840    N=              21
strands: + -
sample: seed=            0    seqfrac=         1
Letter frequencies in dataset:
A 0.310 C 0.190 G 0.190 T 0.310 
Background letter frequencies (from /bio/jvanheld2/matrix_eval/data/bg_freqs/2nt_upstream-noorf_Escherichia_coli_K12-ovlp-2str.meme_bg):
A 0.294 C 0.206 G 0.206 T 0.294 
********************************************************************************


********************************************************************************
MOTIF  1	width =   22   sites =  21   llr = 280   E-value = 2.6e-046
********************************************************************************
--------------------------------------------------------------------------------
	Motif 1 Description
--------------------------------------------------------------------------------
Simplified        A  4627::::8:7374644:a:14
pos.-specific     C  ::::a::113:1:3124a::31
probability       G  ::12::a1:1111::11::a22
matrix            T  5361:a:8:61522221:::43

         bits    2.3     * *          *    
                 2.0     * *          * *  
                 1.8     ***          ***  
                 1.6     ***          ***  
Information      1.4     ***          ***  
content          1.1     ***          ***  
(19.2 bits)      0.9     ****         ***  
                 0.7 *  *******       ***  
                 0.5 *********** ***  ***  
                 0.2 *************** ***** 
                 0.0 ----------------------

Multilevel           TATACTGTATATAAAAACAGTA
consensus            ATAG     C A CTCC   CT
sequence                          T T      
                                           
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
	Motif 1 sites sorted by position p-value
--------------------------------------------------------------------------------
Sequence name            Strand  Start   P-value                     Site       
-------------            ------  ----- ---------            ----------------------
LexA|2749749-2749768         +      9  5.40e-11   CCAGCCTC TTTACTGTATATAAAACCAGTT TATACTGTAC
LexA|3851322-3851341         -     11  8.07e-10   TATAGAAG TTTACTGTATAAATAAACAGTA ATATTTGGAC
LexA|812655-812674           -     11  2.09e-09   CAACAAAT TATACTGGATAAAAAAACAGTT CATCACCATA
LexA|4255112-4255131         +      9  3.54e-09   CTCACAGC ATAACTGTATATACACCCAGGG GGCGGAATGA
LexA|1229951-1229970         +      9  3.54e-09   AAGAACAG ACTACTGTATATAAAAACAGTA TAACTTCAGG
LexA|4255091-4255110         +      9  4.19e-09   AATCGCCT TTTGCTGTATATACTCACAGCA TAACTGTATA
LexA|2749771-2749790         +      9  5.80e-09   AACCAGTT TATACTGTACACAATAACAGTA ATGGTTTTTC
LexA|3995930-3995949         +      9  4.18e-08   TAATCAGC AAATCTGTATATATACCCAGCT TTTTGGCGGA
LexA|4271978-4271997         -     11  5.35e-08   TGCATTCC AATACTGTATATTCATTCAGGT CAATTTGTGT
LexA|1944050-1944069         +      9  6.79e-08   GATAAAAA AATGCTGGATAGATATCCAGCG AAGGATGAAG
LexA|832259-832278           -     11  6.79e-08   AAACCTGA AATACTGTATAAACAGCCAATA TTGTGGCATT
LexA|2821851-2821870         +      9  9.59e-08   GAAGCAAT TATACTGTATGCTCATACAGTA TCAAGTGTTT
LexA|65834-65853             -     11  1.83e-07   GGGCAGTA ATGACTGTATAAAACCACAGCC AATCAAACGA
LexA|1020162-1020181         -     14  2.24e-07      CTGGA TGTACTGTACATCCATACAGTA ACTCACAGGG
LexA|4255050-4255069         -     11  6.76e-07   TGGAACCA TAAACTGCACAATAAACCAGAG ATTTATCGAA
LexA|932351-932370           -     11  2.13e-06   CCAGTACT GTTGCTGTATGGATTAACAGGA GTGTAATCAA
LexA|738568-738587           -     11  4.05e-06   CCACGGCG ATAACTGTCGATAAGCGCAGCC AGCTGCTGGC
LexA|738659-738678           -     11  1.16e-05   TTCGAAAT AACGCTGCCCTGAAAGCCAGGC GTCAGGATAA
LexA|3208751-3208770         +     10  1.26e-05  ATTTTGAAA TAAGCTGGCGTTGATGCCAGCG GCAAACCGA 
LexA|1944102-1944121         -     12  1.26e-05    AATAAAT TATACTGTGCCATTTTTCAGTT CATCGAGACA
LexA|1229931-1229950         -     11  4.68e-05   ATATACAG TAGTCTGTTCTTGCCAGCAGAT CAATACTGAT
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
	Motif 1 block diagrams
--------------------------------------------------------------------------------
SEQUENCE NAME            POSITION P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
LexA|2749749-2749768              5.4e-11  8_[+1]_10
LexA|3851322-3851341              8.1e-10  10_[-1]_8
LexA|812655-812674                2.1e-09  10_[-1]_8
LexA|4255112-4255131              3.5e-09  8_[+1]_10
LexA|1229951-1229970              3.5e-09  8_[+1]_10
LexA|4255091-4255110              4.2e-09  8_[+1]_10
LexA|2749771-2749790              5.8e-09  8_[+1]_10
LexA|3995930-3995949              4.2e-08  8_[+1]_10
LexA|4271978-4271997              5.3e-08  10_[-1]_8
LexA|1944050-1944069              6.8e-08  8_[+1]_10
LexA|832259-832278                6.8e-08  10_[-1]_8
LexA|2821851-2821870              9.6e-08  8_[+1]_10
LexA|65834-65853                  1.8e-07  10_[-1]_8
LexA|1020162-1020181              2.2e-07  13_[-1]_5
LexA|4255050-4255069              6.8e-07  10_[-1]_8
LexA|932351-932370                2.1e-06  10_[-1]_8
LexA|738568-738587                4.1e-06  10_[-1]_8
LexA|738659-738678                1.2e-05  10_[-1]_8
LexA|3208751-3208770              1.3e-05  9_[+1]_9
LexA|1944102-1944121              1.3e-05  11_[-1]_7
LexA|1229931-1229950              4.7e-05  10_[-1]_8
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
	Motif 1 in BLOCKS format
--------------------------------------------------------------------------------
BL   MOTIF 1 width=22 seqs=21
LexA|2749749-2749768     (    9) TTTACTGTATATAAAACCAGTT  1 
LexA|3851322-3851341     (   11) TTTACTGTATAAATAAACAGTA  1 
LexA|812655-812674       (   11) TATACTGGATAAAAAAACAGTT  1 
LexA|4255112-4255131     (    9) ATAACTGTATATACACCCAGGG  1 
LexA|1229951-1229970     (    9) ACTACTGTATATAAAAACAGTA  1 
LexA|4255091-4255110     (    9) TTTGCTGTATATACTCACAGCA  1 
LexA|2749771-2749790     (    9) TATACTGTACACAATAACAGTA  1 
LexA|3995930-3995949     (    9) AAATCTGTATATATACCCAGCT  1 
LexA|4271978-4271997     (   11) AATACTGTATATTCATTCAGGT  1 
LexA|1944050-1944069     (    9) AATGCTGGATAGATATCCAGCG  1 
LexA|832259-832278       (   11) AATACTGTATAAACAGCCAATA  1 
LexA|2821851-2821870     (    9) TATACTGTATGCTCATACAGTA  1 
LexA|65834-65853         (   11) ATGACTGTATAAAACCACAGCC  1 
LexA|1020162-1020181     (   14) TGTACTGTACATCCATACAGTA  1 
LexA|4255050-4255069     (   11) TAAACTGCACAATAAACCAGAG  1 
LexA|932351-932370       (   11) GTTGCTGTATGGATTAACAGGA  1 
LexA|738568-738587       (   11) ATAACTGTCGATAAGCGCAGCC  1 
LexA|738659-738678       (   11) AACGCTGCCCTGAAAGCCAGGC  1 
LexA|3208751-3208770     (   10) TAAGCTGGCGTTGATGCCAGCG  1 
LexA|1944102-1944121     (   12) TATACTGTGCCATTTTTCAGTT  1 
LexA|1229931-1229950     (   11) TAGTCTGTTCTTGCCAGCAGAT  1 
//

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
	Motif 1 position-specific scoring matrix
--------------------------------------------------------------------------------
log-odds matrix: alength= 4 w= 22 n= 399 bayes= 4.16993 E= 2.6e-046 
    55  -1104   -211     83 
    96   -211   -211     18 
   -30   -211   -111    108 
   118  -1104     21   -162 
 -1104    228  -1104  -1104 
 -1104  -1104  -1104    177 
 -1104  -1104    228  -1104 
 -1104   -111    -53    137 
   137    -53   -211   -262 
 -1104     47   -111    108 
   128   -211   -111   -104 
    -4   -111    -53     70 
   118   -211   -111    -62 
    55     69  -1104    -30 
   108   -111   -211    -30 
    38     21    -53    -30 
    55     88   -111   -162 
 -1104    228  -1104  -1104 
   177  -1104  -1104  -1104 
  -262  -1104    221  -1104 
  -162     47    -12     55 
    38    -53    -12     -4 
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
	Motif 1 position-specific probability matrix
--------------------------------------------------------------------------------
letter-probability matrix: alength= 4 w= 22 nsites= 21 E= 2.6e-046 
 0.428571  0.000000  0.047619  0.523810 
 0.571429  0.047619  0.047619  0.333333 
 0.238095  0.047619  0.095238  0.619048 
 0.666667  0.000000  0.238095  0.095238 
 0.000000  1.000000  0.000000  0.000000 
 0.000000  0.000000  0.000000  1.000000 
 0.000000  0.000000  1.000000  0.000000 
 0.000000  0.095238  0.142857  0.761905 
 0.761905  0.142857  0.047619  0.047619 
 0.000000  0.285714  0.095238  0.619048 
 0.714286  0.047619  0.095238  0.142857 
 0.285714  0.095238  0.142857  0.476190 
 0.666667  0.047619  0.095238  0.190476 
 0.428571  0.333333  0.000000  0.238095 
 0.619048  0.095238  0.047619  0.238095 
 0.380952  0.238095  0.142857  0.238095 
 0.428571  0.380952  0.095238  0.095238 
 0.000000  1.000000  0.000000  0.000000 
 1.000000  0.000000  0.000000  0.000000 
 0.047619  0.000000  0.952381  0.000000 
 0.095238  0.285714  0.190476  0.428571 
 0.380952  0.142857  0.190476  0.285714 
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
	Motif 1 regular expression
--------------------------------------------------------------------------------
[TA][AT][TA][AG]CTGTA[TC]A[TA]A[ACT][AT][ACT][AC]CAG[TC][AT]
--------------------------------------------------------------------------------




Time  0.06 secs.

********************************************************************************


********************************************************************************
SUMMARY OF MOTIFS
********************************************************************************

--------------------------------------------------------------------------------
	Combined block diagrams: non-overlapping sites with p-value < 0.0001
--------------------------------------------------------------------------------
SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
LexA|65834-65853                 6.95e-06  10_[-1(1.83e-07)]_8
LexA|738568-738587               1.54e-04  10_[-1(4.05e-06)]_8
LexA|738659-738678               4.39e-04  10_[-1(1.16e-05)]_8
LexA|812655-812674               7.96e-08  10_[-1(2.09e-09)]_8
LexA|832259-832278               2.58e-06  10_[-1(6.79e-08)]_8
LexA|932351-932370               8.09e-05  10_[-1(2.13e-06)]_8
LexA|1020162-1020181             8.51e-06  13_[-1(2.24e-07)]_5
LexA|1229931-1229950             1.78e-03  10_[-1(4.68e-05)]_8
LexA|1229951-1229970             1.35e-07  8_[+1(3.54e-09)]_10
LexA|1944050-1944069             2.58e-06  8_[+1(6.79e-08)]_10
LexA|1944102-1944121             4.80e-04  11_[-1(1.26e-05)]_7
LexA|2749749-2749768             2.05e-09  8_[+1(5.40e-11)]_10
LexA|2749771-2749790             2.21e-07  8_[+1(5.80e-09)]_10
LexA|2821851-2821870             3.64e-06  8_[+1(9.59e-08)]_10
LexA|3208751-3208770             4.80e-04  9_[+1(1.26e-05)]_9
LexA|3851322-3851341             3.07e-08  10_[-1(8.07e-10)]_8
LexA|3995930-3995949             1.59e-06  8_[+1(4.18e-08)]_10
LexA|4255050-4255069             2.57e-05  10_[-1(6.76e-07)]_8
LexA|4255091-4255110             1.59e-07  8_[+1(4.19e-09)]_10
LexA|4255112-4255131             1.35e-07  8_[+1(3.54e-09)]_10
LexA|4271978-4271997             2.03e-06  10_[-1(5.35e-08)]_8
--------------------------------------------------------------------------------

********************************************************************************


********************************************************************************
Stopped because nmotifs = 1 reached.
********************************************************************************

CPU: genomix

********************************************************************************
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
print $query->hidden(-name=>'matrix_format',-default=>'meme');
print $query->hidden(-name=>'kfold',-default=>$default{kfold});

print $query->hidden(-name=>'tag1',-default=>'positive_set');
print $query->hidden(-name=>'sequence1',-default=>$demo_seq1);
print $query->hidden(-name=>'permutation1',-default=>1);
print $query->hidden(-name=>'scanopt1',-default=>'');

print $query->hidden(-name=>'tag2',-default=>'negative_set');
print $query->hidden(-name=>'sequence2',-default=>$demo_seq2);
print $query->hidden(-name=>'markov_order',-default=>$demo_markov);
print $query->hidden(-name=>'scanopt2',-default=>'');

print $query->submit(-label=>"DEMO");
print "</b></td>\n";
print $query->end_form;


print "<td><b><a href='help.matrix-quality.html'>MANUAL</A></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

