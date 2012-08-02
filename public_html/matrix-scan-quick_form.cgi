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

$default{sequence_file} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence} = ""; ### [-f <Name of sequence file---default: standard input>]
$default{sequence_format} = "fasta"; ### automatic conversion from any format to wc

$default{origin}="end";

$default{bg_method}="bginput";
$checked{$default{bg_method}} = "CHECKED";
$default{markov_order} = "0";
$default{organism} = "Saccharomyces cerevisiae";
$default{matrix_format} = "tab";
$default{pseudo_counts} = 1;
$default{pseudo_distribution} = "pseudo_prior";
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{decimals} = "1";
$default{analysis_type} = "analysis_sites";
$checked{$default{analysis_type}} = "CHECKED";

## Return fields
$default{return_field} = "sites";
$default{return_site_limits} = "on";

## Threshold values for site detection
$default{thresh_field} = "weight";
$default{thresh_value} = "1";

### print the form ###
&RSA_header("matrix-scan QUICK and SIMPLE");
&ListParameters() if ($ENV{rsat_echo} >= 2);

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
}

&ReadMatrixFromFile();

### head
print "<center>";
print "Scan a DNA sequence with a profile matrix<br>\n";
print "This quick version was programmed by <a href='mailto:defrance@bigre.ulb.ac.be'>Matthieu Defrance</a>, Web interface by <A HREF='mailto:morgane\@bigre.ulb.ac.be (Morgane Thomas-Chollier)'>Morgane Thomas-Chollier</A><br>\n";
print "</CENTER>";

print "<div align=center>";
print "<b>Citation</b>: <a href='mailto:jturatsi\@bigre.ulb.ac.be (Jean Valery Turatsinze)'>Jean Val&eacute;ry Turatsinze</A>, <A HREF='mailto:morgane\@bigre.ulb.ac.be (Morgane Thomas-Chollier)'>Morgane Thomas-Chollier</A>, <a href='mailto:defrance@bigre.ulb.ac.be'>Matthieu Defrance</a> and <A HREF='mailto:jvanheld\@bigre.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</a> (2008).<br> Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nat Protoc, 3, 1578-1588. <a href='http://www.ncbi.nlm.nih.gov/pubmed/18802439'>Pubmed 18802439</a>";
print "</p>";

print "<a href='matrix-scan_form.cgi'><b><font color=red>--> Click here to access the ADVANCED form <-- </font></b></a> <i>(custom background, CRER detection, overrepresentation of sites,...)</i>";
print "</div>";

## demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"matrix-scan.cgi");

################################################################
#### sequence
print "<fieldset>
<legend><b><a href='help.formats.html'>Sequences </a></b></legend>";
&MultiSequenceChoice("",1);
print "</fieldset><p/>";

################################################################
#### Matrix specification
print "<fieldset>
<legend><b><a href='help.convert-matrix.html#io_format'>Matrix </a></b></legend>";

&GetMatrix("consensus"=>0,"no_pseudo"=>1);
print "</fieldset><p/>";

################################################################
## Background model
print "<fieldset>
<legend><b><a href='help.convert-matrix.html#io_format'>Background </a></b></legend>";

my %bg_params =("markov" => 1,
		"bg_input" => 1,
		"bg_window" => 0,
		"markov_message" => 0,
		"simple"=>1,
	       );
&GetBackgroundModel(%bg_params);

print "</fieldset><p/>";

################################################################
#### scanning options

print "<fieldset>
<legend><b>Scanning options</b></legend>";

################################################################
#### origin for calculating positions
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#origin'><B>Sequence Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['start',
				   'center',
				   'end'],
			 -default=>$default{origin});
print "<br/>";

################################################################
## Fields to return + thresholds
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#return_fields'><B>Return</B></A>\n";

my %returns = ("sites" => "sites only",
				"pval" => "sites + pval");

my $Popup = "";
    $Popup .=  "<SELECT NAME='return_field' onChange=\"toggle(this.options[this.selectedIndex].value)\">";
     foreach my $f (keys %returns) {
     	if ($f eq $default{return_field}){
			$Popup .=  "<OPTION  SELECTED VALUE=$f>$returns{$f}</option>\n";
     	} else {
     		$Popup .=  "<OPTION VALUE=$f  >$returns{$f}</option>\n";
     	}
     }
    $Popup .=  "</SELECT>";
    print $Popup;

print "&nbsp;"x2, "<i>Calculating pval is slighly slower</i>";
print "<br/>";

## thresholds	 
print "&nbsp;"x4,  "<A HREF='help.matrix-scan.html#thresholds'><B>Threshold</B></A>\n";
print "<br/>";

print "<table style='padding-left:50px;'>";
print "<tr>";
print "<td style='text-align:right;'>weight score >=</td>";
print "<td rowspan=2>";
print $query->textfield(-name=>'thresh_value',
							    -default=>$default{thresh_value},
							    -size=>5);
print "</td>";
print "<td><i>if return is '<b>sites only</b>', the threshold is set on the weight score</i></td>";
print "</tr>";

print "<tr>";
print "<td style='text-align:right;'>pval <=</td>";
print "<td><i>if return is '<b>sites + pval </b>', the threshold is set on the pval</i></td>";
print "</tr>";
print "</table>";


print "</fieldset>";

## add the -quick option
print $query->hidden(-name=>'quick',-default=>'CHECKED');
print $query->hidden(-name=>'analysis_type',-default=>$default{analysis_type});
print $query->hidden(-name=>'pseudo_counts',-default=>$default{pseudo_counts});
print $query->hidden(-name=>'decimals',-default=>$default{decimals});
print $query->hidden(-name=>'strands',-default=>"both");
print $query->hidden(-name=>'bg_pseudo',-default=>$default{bg_pseudo});
print $query->hidden(-name=>'return_site_limits',-default=>$default{return_site_limits});
################################################################
### send results by email or display on the browser
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
print $query->start_multipart_form(-action=>"matrix-scan-quick_form.cgi");

$demo_sequence=">dm3_flyBaseGene_CG2328-RA
GCCACTGAAGCTGTGGGTTTGCTCCTGCCGAGCGAATCCAACGCGAGTAG
GGTCCCATTCGGGGCCCGAGTAGCCAGAGTCCTGCAGCTCACTCGAAACC
GCCACTCACCGTGGCTAATTGCCCATCAATAAAGGGCCCGGGCAGTGAGG
AATTCCTCCGAAAGTCGGGTCCTCCGTTCTCCAGCCGAAGATTTTTTCGA
GCAACCAAAATATTATGGTGTGCCCCGCTGTTCTCGCACAGTCAGCGCGA
ATTTGCTGCGGTGAGTCGATGCTGTTTCGCAGGACCTTCTTCCATTTTCG
TCTCCCTCTTGCTCAGCCTGTCCCTGTTCCTCTGCAGTTCCCTATCTCCT
GATGCCTGTGCTCCTTTGGCGGCACTGTGTCCTGTCGTCGTTGTTTTCCT
GTGATTTGACATGTCTGTTAGCAGGATGCCTGACCCTGAGGCCGAGCCCT
GGTCTCAGTGTCCACTGTTCCACTTTGATGTGATTCGTCAGTGCGGTGGA
CTACTGCTACTGCTCTCTTGCTGGACTGCGTCTTGAGTCCTGTTCGGCTG
CCCCCTCCCGTGACCTCTGACCCTGCACTCTGCGGCTTTCCAGCGGCGTT
TGTTGGCGAATCTGACCCCGAGCTCCTGCTGCTCCTTCGCTCCTTCGCTC
CTTCTCCGCATCTCCGCTCTTTGGACTTCGTACGAATCAAAATTGGTCAC
AGCACCGAGTGAATTGCCCCGGAGACCGCAATGCGCTGTATTTATAGTAA
ACGTGTCCGATTGATTTGGCCACCCGTGGCGGCTCTGTCACAGATGCCTC
AATTTGCATCTATCGAATGGTTTACATGGCTCTAAAAAGGTACCTCGATG
GGTTGGTCACAATGTGGTGGCCTCTCAACATTGCAAGGCTCTTACTTGTG
AATTATTAAGTTATTAACTGCTGCGATGTAAGTCATGGCAGTTTCTGTTT
TCTTTATAGGATATATATAGGAAGGATTAAAGGAGGCATGTACAATAATA
TGAGTATGATTTAGCTCAAATTCCAAATATGATAAAAGTACAAAGCATAC
GATAATATAATCAAATTACGCTGACAATCACGATAATGTTCTTGTAGTAG
TATTTGTGTAATATTTATGTTTTTTTAAGATAAGAAACGGTAATAAAATC
CACGTAAGTGTAAAAAATGGATGCCCTAATCTATGCCATGATGTGTTCTA
CTTTCGAGATTTCGCCTCTGCCCTCATTGATGGTTTCCCGGGGCTACTTG
GCCCAAAAATCCCGGCCGTCCAAAAAGACGATCCTTAAAAAAGAAACCGC
TAATCATTGGGCCGCACAAAGAGCGGACAATCGCTCACCTAATTATTTGG
CCCGATTGTGAGGAGCGGACAGTCGGCTCGTGGACGCTTTTTGTGGCCTC
TTTTTGTTTCGACAAAAAGCGAGCCAATTTTTTTTCTTTCTGGGCCACTT
TGTTGCTCTTTTTATGAGTTTTTTCCATTGTCAGTTTTTCCGGGCCTGTC
TCGCAGCCCTCGATTCCCGCGATGCCTGCCCTACAAACCTCCTAATTACG
GCAGTTAGTCGTTGTCCGGGACAGGAGAGTATGCGGAAGGACATGCGTGA
GTTTATTGCCCGCTCGAATTTCCACTAAAAATTGGGCCGAAAAAAAAACA
ACTAGGTAGGACTAGGAACTGCAAACTAGCAAAGCGGACGCGCCTTTTTA
TTGGTGCACCTTCGGCGGAACCGCAGGATAACAGCAGTAAAAGCGACGAC
GAGGACACAAGGATCCTCGAAATCGAGAGCGACCTCGCTGCATTAGAAAA
CTAGATCAGTTTTTTGTTTTGGCCGACCGATTTTTGTGCCCGGTGCTCTC
TTTACGGTTTATGGCCGCGTTCCCATTTCCCAGCTTCTTTGTTCCGGGCT
CAGAAATCTGTATGGAATTATGGTATATGCAGATTTTTATGGGTCCCGGC
GATCCGGTTCGCGGAACGGGAGTGTCCTGCCGCGAGAGGTCCTCGCCGGC
GATCCTTGTCGCCCGTATTAGGAAAGTAGATCACGTTTTTTGTTCCCATT
GTGCGCTTTTTTCGCTGCGCTAGTTTTTTTCCCCGAACCCAGCGAACTGC
TCTAATTTTTTAATTCTTCACGGCTTTTCATTGGGCTCCTGGAAAAACGC
GGACAAGGTTATAACGCTCTACTTACCTGCAATTGTGGCCATAACTCGCA
CTGCTCTCGTTTTTAAGATCCGTTTGTTTGTGTTTGTTTGTCCGCGATGG
CATTCACGTTTTTACGAGCTCGTTCCTTCGGGTCCAAAATTATGCCAGTT
TGTTTTGTCTCTGGCAATTATTGGAAATTTCATTGGGTCGATTTCGCTGC
CTTCCTTGCTCTTCCCTTGAGAAAAGTGAATAGGTTGTGCCATAAAAATC
GCTGCTCCTGAAGACCAAATGAAATGGATTTGTGTAAGCATTAAAAACGC
GAGGCAAGCCCCAAGATTCCTCCACTGCTTTTTTTATATTGCCCACTGCT
AAATGCAGCTAATTCGTCGATTGTTTAAAAATTAAATTACTTATGTTGCC
ATTCATACATCCCCTCACATTTTATGGCCATTTGAGTGCGGGGTGCACAG
TTCTGTCTTAAGTGGCGGATGGAAACCACCACATTTACTCGAGGGATGAT
GTGCTCTAATATCTCCTCATCAAATGGGATGGTTTCTATGGAAAGGCAAA
ATCGTTGTAAAGTGAGGCGGAGTTAAAAAATACCTTGTTATAGCCTTTTT
AAAATAACACAAGATCGTTCGAATTGACTAGAAATATCAAAGTCTTTTTG
TATTGAAGCGAGTGTAGTCTCAATTTATGCTTAATTTTAAGAAATACATC
TCTTTATTAGCCCCAAAATGAAACAAATGGTCTACTAATTAAGCAAGTCA
ACAGAATTTTTATGCAATTATTCAAAATGAAATAATATATACATAAGATG
TTTTTGGGAATCTGTCATGGGGTTTCTGAAATAGGTTTGCCAAACAAATT
TTAAGTATAAATGTATACATATGTCAACTAATAAATTTAGCAAATAAAAT
GTACCTGCAAGTATCTATAAATTTATTGGACCAATTTTGTGTAAAAAACT
GAACTGGCACTCTTCCCAAGAATGGGACTTCGAGGACTCCTTGCTGAATC
ACTTACTCAACCCATTCCAACTCATCCAATCCGCGCAATCATCATAAATT
TTGGCCTTTTTGTTGTAATTGTTTTATGGCAGAAATTACTCAATCATCAA
GCATAATTCCCTCGTTTTCGCCGTTTTATTGCCAATTTTTGCACTGCCTT
TGCCTTTTTCCCGCCCTTTCCTCAGCGTTTTGCGAATCTTTGCCGGCATT
TCTATTGCGCGGACAATCCGGCCAGTGTGTTGGCCATTTACTTGCCATGA
TGACGGGCATAATCAGCGAGATCGGCGCTTTGTGAGTGCAGAATGTGCAA
TAAAGCGGCAACAATCGGCAGGGATTCGCCTTCCCATATTCCGGGTATTG
CCGGCCCGGGAAAATGCGAAAGTGTTTGCGGATCGAGATGGAAGATAGAG
GATTGAGTATTGAAACGAGGAAGGTACTTCCGCCGGCGGACACTTTCGCC
TAACCAAGCCAATCCAACCCATCCCAATCCAATCCAACCCACCCGATCGC
CATAAAGGGTATTTACTGTCGCTGCCGCAGAGCCTCGCTTGACGACTTAA
CCCAAGCGGTCGTTTCGCGTCCATTCTCCGGACGGAGTCAAAGACAAAGG
CCGGCGGAGGTGGACAATAGGCAAGGTTGTTGCTTGTGGGTAGGGTTTGA
GCTATGAGCTATGAGCTGTGAGCTGTTAGCCCTGAACCCCGAACCTCGAG
AATTGAACCTTTCCCGGGGCAAGAAGGCTTGCATGTGGGCCTTTTCCAGG
TCGGCCAGTAGGTAGAGTTGTTGCGATGCGGCTATGCCGGGCGAGTTAAT
GCCAATGCAAATTGCGGGCGCAATATAACCCAATAATTTGAAGTAACTGG
CAGGAGCGAGGTATCCTTCCTGGTTACCCGGTACTGCATAACAATGGAAC
CCGAACCGTAACTGGGACAGATCGAAAAGCTGGCCTGGTTTCTCGCTGTG
TGTGCCGTGTTAATCCGTTTGCCATCAGCGAGATTATTAGTCAATTGCAG
TTGCAGCGTTTCGCTTTCGTCCTCGTTTCACTTTCGAGTTAGACTTTATT
GCAGCATCTTGAACAATCGTCGCAGTTTGGTAACACGCTGTGCCATACTT
TCATTTAGACGGAATCGAGGGACCCTGGACTATAATCGCACAACGAGACC
GGGTTGCGAAGTCAGGGCATTCCGCCGATCTAGCCATCGCCATCTTCTGC
GGGCGTTTGTTTGTTTGTTTGCTGGGATTAGCCAAGGGCTTGACTTGGAA
TCCAATCCCGATCCCTAGCCCGATCCCAATCCCAATCCCAATCCCTTGTC
CTTTTCATTAGAAAGTCATAAAAACACATAATAATGATGTCGAAGGGATT
AGGGGCGCGCAGGTCCAGGCAACGCAATTAACGGACTAGCGAACTGGGTT
ATTTTTTTGCGCCGACTTAGCCCTGATCCGCGAGCTTAACCCGTTTTGAG
CCGGGCAGCAGGTAGTTGTGGGTGGACCCCACGATTTTTTTGGCCAAACC
TCCAAGCTAACTTGCGCAAGTGGCAAGTGGCCGGTTTGCTGGCCCAAAAG
AGGAGGCACTATCCCGGTCCTGGTACAGTTGGTACGCTGGGAATGATTAT
ATCATCATAATAAATGTTTTGCCCAACGAAACCGAAAACTTTTCAAATTA
AGTCCCGGCAACTGGGTTCCCATTTTCCATTTTCCATGTTCTGCGGGCAG
GGGCGGCCATTATCTCGCTACAGCAGTTCCCAAATGGTTATGGCTGGACA
CCCCTGCCGCCGCTCCAACGGGGTGGATGAAGCCCCCAAAACCCGAAAGT
CATGGCAGCCATGGCAGTGTGGGGCTGTTAAACGTGCGGCATAATATTAA
GACTTCATAAAAGCGCAAATAATTCGCTGGCAGGCGATCGATAATACATA
CATACAAATATATAGTGGGATACACACACTCTCTGCCGGCAAACACACAC
CACCCGACCCGACTGAGCGGCATAATGCCATATCATTCTTGATGAAGCCG
ATAAAATCCCATTATTAAGGGGGCCCGCCCGTCCCGCTCGCTCCTGCGGA
GCAACCGCCTGCGGGCGGGCGAGACAAAAGATTCGCTCATCCGCTATGAA
TACCAAATCGGAACTCTCTCTCTCTCCAGCTCGGGAGTGCCATGGCCAGC
ATGGCCAGGACCTCCTCATGGTCCTGCCGAGCAGAGAACGCGGCTCCATC
CCGCTGCTCCGGGTCCTGCTCCTCCGCTTTGTCCCGCCTCGTTATCGCCG
CTCAGCACCGAGAGCACAGCAGCGCATCCACTCTCAGCACCGCACGATTA
GCACCGTTCCGCTCAGGCTGTCCCGCTCGCACCTGCCTGGGTCGCTGCGA";

$demo_matrix = "
; Matrices derived from OregAnno by J.V. Turatsinze (2007).
;
; MATRIX 1/12 : Kr
;
AC  Kr
XX
ID  Kr
XX
PO       A     C     G     T
1       30     4     4     6
2       37     0     2     5
3        0    35     3     6
4        3    37     0     4
5        1    41     0     2
6        5     9    11    19
7        4     1     7    32
8        2     4     0    38
XX
//
;
; MATRIX 5/12 : eve
;
AC  eve
XX
ID  eve
XX
PO       A     C     G     T
1        0     8     1     0
2        5     0     4     0
3        3     4     2     0
4        6     2     1     0
5        1     3     0     5
6        2     3     3     1
7        7     0     2     0
8        0     0     0     9
9        0     1     0     8
10       8     1     0     0
11       4     0     4     1
12       2     5     2     0
13       3     0     6     0
14       1     6     2     0
15       0     7     0     2
XX
//
";

$descr="<H4>Comment on the demonstration example : </H4><blockquote class ='demo'>In this demonstration, we will analyse
the promoter of Drosophila melanogaster even-skipped gene (eve). We will scan the 5500 bp sequence upstream the transcription start site with
matrices representing the binding specificity of 2 transcription factors known to regulate eve. These matrices were built from
binding sites annotated in the <a target=_blank href='http://www.oreganno.org'>ORegAnno</a> database by Jean-Valery Turatsinze.<p/>";

## demo 1
print "<TD><B>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr."The program will return individual matches, i.e. sequence segments scoring above the predefined threshold. In this example, threshold is set on the Pval.
</blockquote>");
print $query->hidden(-name=>'bg_method',-default=>'bgfile');
print $query->hidden(-name=>'thresh_field',-default=>'pval');
print $query->hidden(-name=>'thresh_value',-default=>'1e-4');
print $query->hidden(-name=>'bgfile',-default=>'CHECKED');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'markov_order',-default=>'0');
print $query->hidden(-name=>'organism',-default=>'Drosophila_melanogaster');
print $query->hidden(-name=>'return_field',-default=>'pval');
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'origin',-default=>'end');
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->submit(-label=>"DEMO");
print "</B></TD>";
print $query->end_form;


print "<TD><B><A HREF='help.matrix-scan.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_matrix-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

&ListParameters() if ($ENV{rsat_echo} >= 2);
&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->end_html;

exit(0);
