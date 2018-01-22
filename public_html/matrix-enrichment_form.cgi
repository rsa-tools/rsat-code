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
$main::quality=0; 
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
$default{tag3} = "sequence_set3";
$default{tag4} = "sequence_set4";
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
&RSA_header("matrix-enrichment", "form");
print "<center>";
print "Evaluate the enrichment of a set of motifs in one or several sequence sets. (The web version only allows for four sequence sets)</p>\n";
print "The most classical use of the program is to indentify transcription factor binding sites that could be enriched or depleted in one or several sequence sets.<p>\n";
print "<p>Program developed by <a target='_top' href='http://folk.uio.no/jamondra/'>Jaime A Castro-Mondragon</a>, \n";
print " <a target='_top' href='http://www.ibens.ens.fr/spip.php?article94&lang=en'>Samuel Collombet</A>,\n";
print " <a target='_top' href='http://liigh.unam.mx/amedina/'>Alejandra Medina-Rivera</A>,\n";
print " <a target='_top' href='http://morgane.bardiaux.fr/'>Morgane Thomas-Chollier</A>,\n";
print "and <a target='_top' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/ '>Jacques van Helden</A>.</p>\n";
print "</center>\n";


## demo description
#print $default{demo_descr};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"matrix-enrichment.cgi", -onreset=>"resetHandler()");


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
print "<p><font color='orange'>Only the first 50 matrices will be taken into acount</font></p>";

print "</fieldset><p/>";


################################################################
#### Sequence specification

print "<fieldset>
<legend><b><a class='iframe' href='help.formats.html'>2 - Sequences </a></b></legend>";


print "<h2> Mandatory Sequence Set </h2>";

&SeqBoxMQ(1);
print "<hr>";

print "<h2> Optional Sequence Set 1 </h2>";
&SeqBoxMQ(2);

print "<h2> Optional Sequence Set 2 </h2>";
&SeqBoxMQ(3);

print "<h2> Optional Sequence Set 3 </h2>";
&SeqBoxMQ(4);

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

open(my $fh, "demo_files/Ballester_etal_elife_2014_4TFs_motifs.tf");
while(my $row = <$fh>){
    chomp $row;
    $demo_matrix .= $row;
    $demo_matrix .= "\\n";
}
close ($fh);

$demo_seq1_url= $ENV{rsat_www}."/demo_files/Ballester_etal_elife_2014_hg18_cebpa_singletons.fa";
$demo_seq2_url= $ENV{rsat_www}."/demo_files/Ballester_etal_elife_2014_hg18_foxa1_singletons.fa";
$demo_seq3_url= $ENV{rsat_www}."/demo_files/Ballester_etal_elife_2014_hg18_hnf4a_singletons.fa";
$demo_seq4_url= $ENV{rsat_www}."/demo_files/Ballester_etal_elife_2014_hg18_hnf6_singletons.fa";

$demo_bg_url= $ENV{rsat_www}."/demo_files/all_human_ENCODE_DNAse_mk1_bg.ol";

print '<script>
function setDemo(demo_matrix, demo_seq1_url, demo_seq2_url, demo_seq3_url, demo_seq4_url, demo_bg_url ){
    $("#reset").trigger("click");
    
    descr = "<H4>Comment on the demonstration example : </H4><blockquote class =\'demo\'>In this demonstration, we will assess the enrichment of four liver Transcription Factors CEBP-alpha, FOXA1, HNF4 and HNF6, in the reported singleton sites of each TF. </p> \
   Singleton sequences, are ChIP-seq assayed regions where only one of the four TFs had signal and there was no other TF in a surrounding 300bp window.</p>\
    These data was published in the Ballestar et al, eLife, 20015 article <a target=\'_top\' href=\'https://www.ncbi.nlm.nih.gov/pubmed/25279814\'>[Pubmed]</A>. </p>    </blockquote> \ <p> \.";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;

    html_title.value = "Zoo-ChIP_liver_Transcription_Factors_eLife_2015";
    matrix.value = demo_matrix;
    matrix_format.value = "transfac";

    tag1.value = "CEBP-alpha_singleton_sites";
    sequence_url1.value  = demo_seq1_url ;

    tag2.value = "FOXA1_singleton_sites";
    sequence_url2.value = demo_seq2_url ;

    tag3.value = "HNF4_singleton_sites";
    sequence_url3.value  = demo_seq3_url ;

    tag4.value = "HNF6_singleton_sites";
    sequence_url4.value  = demo_seq4_url ;
    
    $("#url").prop("checked",true);
    bgmodel_url.value = demo_bg_url ;
    
}
function resetHandler(){
    $("#db_choice").val("").change();
}
</script>';

$demo_markov=1;

print "<td><b>";

#print '<button type="button" onclick="setDemo('. "'$demo_matrix'" .',' . "'$demo_seq1_url'" .',' . "'$demo_seq2_url'" .',' . "'$demo_seq3_url'" .',' . "'$demo_seq4_url'"  .')">DEMO</button>';
print '<button type="button" onclick="setDemo('. "'$demo_matrix'" .',' . "'$demo_seq1_url'" .',' . "'$demo_seq2_url'" .',' . "'$demo_seq3_url'" .',' . "'$demo_seq4_url'" . ',' ."'$demo_bg_url'" .')">DEMO</button>';


################################################################
## Epromoters demo
################################################################

open(my $fh, "demo_files/Capstarr_selected_TFs.tf");
while(my $row = <$fh>){
    chomp $row;
    $demo_2_matrix .= $row;
    $demo_2_matrix .= "\\n";
}
close ($fh);

$demo_2_seq1_url= $ENV{rsat_www}."/demo_files/CapStarrseq_Active_Prom_HELA_merge_IP.fasta";
$demo_2_seq2_url= $ENV{rsat_www}."/demo_files/CapStarrseq_Active_Prom_K562_merge_IP.fasta";
$demo_2_seq3_url= $ENV{rsat_www}."/demo_files/CapStarrseq_InactiveProm_FDR95_All_samples.fasta";
$demo_2_seq4_url="";
$demo_2_bg_url= $ENV{rsat_www}."/demo_files/All_hpromoterRegions_BG_model_mkv_2.oligos";
$demo_markov=2;


print '<script>
function setDemo2(demo_matrix, demo_seq1_url, demo_seq2_url, demo_seq3_url, demo_seq4_url, demo_bg_url ){
    $("#reset").trigger("click");
      
    demo_descr.innerHTML = "<H4>Comment on the demonstration example : </H4><blockquote class =\'demo\'>In this demonstration, we will assess the enrichment of seven Transcription Factors (JUN-FOS, YY, SP1, USF2, REST, DDIT, FCP2), in human promoters displaying enhancer activity (Epromoters), the enrichment is measured in three sets of sequences: HeLa cell, K562 and a set of promoters without enhancer activity. </p> \
    These data was published in the Dao et al, Nat Gen, 2017 article <a target=\'_top\' href=\'doi:10.1038/ng.3884\'>[DOI]</A>. </p>    </blockquote> \ <p> \.";
    demo.value = descr;

    html_title.value = "Enrichment in Epromoters";
    matrix.value = demo_matrix;
    matrix_format.value = "transfac";

    tag1.value = "Hela";
    sequence_url1.value  = demo_seq1_url ;

    tag2.value = "K562";
    sequence_url2.value = demo_seq2_url ;

    tag3.value = "Inactive promoters";
    sequence_url3.value  = demo_seq3_url ;

    tag4.value = "";
    sequence_url4.value  = demo_seq4_url ;
    
    $("#url").prop("checked",true);
    bgmodel_url.value = demo_bg_url ;
    
}
function resetHandler(){
    $("#db_choice").val("").change();
}
</script>';


print "<td><b>";

print '<button type="button" onclick="setDemo2('. "'$demo_2_matrix'" .',' . "'$demo_2_seq1_url'" .',' . "'$demo_2_seq2_url'" .',' . "'$demo_2_seq3_url'" .',' . "'$demo_2_seq4_url'" . ',' ."'$demo_2_bg_url'" .')">DEMO 2</button>';

################################################################



print "</b></td>\n";


print "<td><b><a class='iframe' href='help.matrix-enrichment.html'>MANUAL</A></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

