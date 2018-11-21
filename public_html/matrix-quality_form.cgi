#!/usr/bin/env perl
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
$main::quality=1;
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
$default{organism}="";
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
print "<p>Program developed by <a target='_top' href='http://liigh.unam.mx/amedina/index.html'>Alejandra Medina Rivera</a>, \n";
print " <a target='_top' href='http://morgane.bardiaux.fr/'>Morgane Thomas-Chollier</A>,\n";
print "and <a target='_top' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</A>.</p>\n";
print "</center>\n";
print "<b>Citation</b>: Medina-Rivera, A., Abreu-Goodger, C., Salgado-Osorio, H., Collado-Vides, J. and van Helden, J. (2010). Empirical and theoretical evaluation of transcription factor binding motifs. Nucleic Acids Res. 2010 Oct 4. [Epub ahead of print] <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/20923783'>[Pubmed 20923783]</a> <a target='_blank' href='http://nar.oxfordjournals.org/content/early/2010/10/04/nar.gkq710.full.pdf'>[Full text]</a>.";


## demo description
#print $default{demo_descr};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"matrix-quality.cgi", -onreset=>"resetHandler()");


################################################################
#### Matrix specification
print "<hr>";
print "<h2 style='margin-left: 50px;'> Title ";

print $query->textfield(-name=>'html_title', -id=>'html_title',
			 -default=>$default{html_title},
			 -size=>30) ."</h2>";

print "<fieldset> <legend><b><a class='iframe' href='help.convert-matrix.html#io_format'>1 - Matrix </a></b></legend>";


&GetMatrix();
print "<p></p>";
print $query->checkbox(-name=>'matrix_sites',
  		       -checked=>$default{m_sites},
		       -label=>'');

print "&nbsp;Matrix file includes sites";
print "<p><font color='orange'>Only the first matrix will be taken in acount</font></p>";

print "<\p><b>K fold validation</B>&nbsp;";
print $query->popup_menu(-name=>'kfold', -id=>'kfold',
			 -Values=>["none",0,3,4,5,6,7,8,9,10],
			 -default=>$default{kfold});
print "&nbsp;"x5, "<font color='orange'><b>Note:</b> validation requires a matrix with binding sites, in a suitable format (e.g. transfac, meme).</font>";

print "</fieldset><p/>";


################################################################
#### Sequence specification

print "<fieldset>
<legend><b><a class='iframe' href='help.formats.html'>2 - Sequences </a></b></legend>";


print "<h2> Mandatory Sequence </h2>";

&SeqBoxMQ(1);
print "<hr>";

print "<h2> Optional Sequence </h2>";
&SeqBoxMQ(2);

print "</fieldset><p/>";

################################################################
#### Background specifiaction

print "<fieldset>
<legend><b><a class='iframe' href='help.matrix-scan.html#markov_order'>3 - Background </a></b></legend>";
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
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 

$demo_matrix = "";
$demo_seq1 = "";
$demo_seq2 = "";

open(my $fh, "demo_files/matrix_quality_demo_matrix.tf");
while(my $row = <$fh>){
    chomp $row;
    $demo_matrix .= $row;
    $demo_matrix .= "\\n";
}

open($fh, "demo_files/matrix_quality_demo_seq1.fa");
while(my $row = <$fh>){
    chomp $row;
    $demo_seq1 .= $row;
    $demo_seq1 .= "\\n";
}

open(my $fh, "demo_files/matrix_quality_demo_seq2.fa");
while(my $row = <$fh>){
    chomp $row;
    $demo_seq2 .= $row;
    $demo_seq2 .= "\\n";
}

print '<script>
function setDemo(demo_matrix, demo_seq1, demo_seq2){
    $("#reset").trigger("click");
    
    descr = "<H4>Comment on the demonstration example : </H4><blockquote class =\'demo\'>In this demonstration, we will analyse the PSSM of the Transcription Factors LexA and CRP available in RegulonDB. </p> \
    As first sequence set we will input the obtained sequences from the ChIP-chip experiment (Wade et al. Genes Dev. 2005) of \ transcription factor LexA in Escherichia coli K12 .</p>\
    As second sequence set we will use the reported CRP binding sites in the Escherichia coli K12 annotated in \RegulonDB. </p>\ In the results you will observe there is an enrichment of LexA binding sites in the LexA ChIP-chip reported sequence set, \ and since LexA does not usually binds in the same sequences as CRP, you will notice there is no enrichment of LexA binding sites in CRP sequences.\ Hence, you will observe a reciprocal behavior of CRP predicted binding sites enrcihments, enriched in CRP reported bindings sequences, and not enriched in LexA ChIP-chip sequences. </p>    </blockquote> \ <p> \.";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;

    
    html_title.value = "LexA and CTCF matrices from RegulonDB 2015";
    matrix.value = demo_matrix;
    matrix_format.value = "transfac";
    kfold.value = "none";
    tag1.value = "LexA_peaks";
    sequence1.value = demo_seq1;
    permutation1.value = 1;
    
    tag2.value = "CRP_binding_sites";
    sequence2.value = demo_seq2;
    permutation2.value = 1;

    markov_order.value = 1;
    $("#nwd").prop("checked",true);
    
    $("#organism_bg_name").val("Escherichia_coli_GCF_000005845.2_ASM584v2");
    $("#organism_bg").val("Escherichia_coli_GCF_000005845.2_ASM584v2");
}
function resetHandler(){
    $("#db_choice").val("").change();
}
</script>';


$demo_markov=1;

print "<td><b>";
print '<button type="button" onclick="setDemo('. "'$demo_matrix'" .',' . "'$demo_seq1'" . ','. "'$demo_seq2'" .')">DEMO</button>';
print "</b></td>\n";

print "<td><b><a href='sample_outputs/matrix-quality_demo_output/matrix-quality_2018-03-21.175002_synthesis.html'>[Sample Output]</a></B></TD>\n";
print "<td><b><a class='iframe' href='help.matrix-quality.html'>MANUAL</A></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

