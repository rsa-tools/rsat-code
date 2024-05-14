#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program footprint-scan
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
use RSAT::MatrixReader;
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### Read the CGI query
$default{demo_descr1} = "";

################################################################
## Default values for filling the form

## matrix-scan
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "transfac";
$default{pseudo_distribution} = "pseudo_prior";
$checked{$default{pseudo_distribution}} = "CHECKED";

## Background model
$default{markov_order} = "1";

#$default{leaders} = 'checked';
$default{bg_method}="bgfile";
$checked{$default{bg_method}} = "CHECKED";
$default{organism}="";
$default{uth_pvalue} = "1e-4";
$default{taxon} = "Gammaproteobacteria";
$default{uth_occ_th} = "5";
$default{img_format}="jpg";
$default{info_lines}="CHECKED";
$default{pseudo_freq} = "0.01";



### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
   if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
  if ($key eq "visualize"){
  	$checked{$query->param($key)} = "CHECKED";
  }
}



################################################################
### print the form ###

################################################################
### header
&RSA_header("footprint-scan", "form");
print "<center>";
print "Scan promoters of orthologous genes with one or more user-provided motifs to detect enriched (over-represented) motifs.<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print "</center>";
print "</blockquote>\n";

################################################################
## Display the form only if it is relveant for the organisms supported
## on this RSAT instance.
&check_phylo_tools();

################################################################
## Form header


print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'105\')\" id=\"heading105\"><font color='#0D73A7'>Information about footprint-scan</font> </div>\n";
 print "<div id=\"menu105\" class=\"menu_collapsible\">\n";
print "<BLOCKQUOTE>\n";
print "<p>In the search of new putative Transcription Factor Binding Sites one common aproach is scanning a set of regulator regions from one organism with a PSSM, althogth, this search rases a problem of high False Positive Rate (FPR) when this method is not used carefully, and even when an expert uses it, the FPR is not as small as we would like, and an experimental validation must be done in order to decide if a site is in fact a real TFBS, to circumvent the problem we propose, under the asumption of conservation of regulation, to take advantage of the availability of many bacteria genomes.</p>

<p>Nowadays, the goal in bioinformatics is to increase the statistical power when scanning a genome sequences with regulatory motifs, we propose in <i>footprint-scan</i> the use of additional sequence data from related species in order to achieve this goal.</p> ";
print "</BLOCKQUOTE>\n";
print "</div></p>\n";

&ListParameters() if ($ENV{rsat_echo} >=2);



&ReadMatrixFromFile() ;

## demo description
#print $default{demo_descr1};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"footprint-scan.cgi");

################# Matrix input
 &Panel1();


################# Select reference organism, taxon and query genes
 &Panel2();

################# Scanning Parameters
 &Panel3();

################# Drawing parameters
 &Panel4();



################################################################
## Select output mode. Email is preferred since footprint discovery
## may take a while.
print "<p>\n";
&SelectOutput('email');

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;
print "<TD><b>";

################################################################
### data for the demo
$demo_queries = "lexA\\nrecA\\n";

$demo_matrix_file=$ENV{RSAT}."/public_html/demo_files/LexA.2nt_upstream-noorf-ovlp-2str.20.tf";
$demo_matrix = "";
open(my $fh, $demo_matrix_file);
while (my $row = <$fh>){
    chomp $row;
    $demo_matrix .= $row . "\\n";
}

my $demo_org = "Escherichia coli_str._K-12_substr._MG1655 GCF_000005845.2 ASM584v2";
my $org = $demo_org;
$org =~ s/\ /_/g;
print '<script>
function setDemo(demo_matrix){
    $("#reset").trigger("click");
    
    var descr1 = "<H4>Comment on the demonstration example 1 :</H4> \
    <blockquote class =\"demo\"> \
    <p>In this demonstration, we apply <i>footprint-scan<\i> to evaluate the enrichment of LexA binding site in the upstream sequences of two of its target genes: lexA (the factor is auto-regulated) and recA.</p>\n \
    <p> For each query gene, the orthologs are collected at the level of Gammaproteobacteria, their upstream sequences are scanned with the matrix, and the number of observed sites is compared to the random expectation.</p>\n \
    </blockquote>";
    
    demo_descr.innerHTML = descr1;
    
    $("#queries").val("'. $demo_queries . '");
    $("#organism_name").val("Escherichia coli_str._K-12_substr._MG1655 GCF_000005845.2 ASM584v2").trigger("change");
    $("#organism").val("Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2");
    $("#taxon").val("Gammaproteobacteria");
    $("#taxon_name").val("Gammaproteobacteria").trigger("change");
    $("#matrix").val(demo_matrix);
    $("#matrix_format").val("transfac");
    $("#bg_method_bginput").prop("checked", true);
    $("#markov_order").val("0");
}
</script>';
print '<button type="button" onclick="setDemo('. "'$demo_matrix'" .')">DEMO</button>';

print "<td><b><a href='sample_outputs/footprint-scan_demo_output/footprint-scan_Gammaproteobacteria__lexA_recA_2018-03-21.183943_2018-03-21.183943_J8t3qJ/lexA_recA/Gammaproteobacteria/Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2_ASM584v2/all_matrices_report.html'>[Sample Output]</a></B></TD>\n";
print "<td><b><a href='help.footprint-scan.html'>[MANUAL]</a></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



#################### Internal functions #######################

sub Panel1 {

  print "<fieldset>\n<legend><b><a href='help.formats.html'>Matrix </a></b></legend>\n";

  &GetMatrix();

  print '<b><font style="font-size:80%"><a href="tutorials/tut_peak-motifs.html#seq" target="_blank"></a></font></b></br>';
  print "</fieldset><p/>";
}

##########################################
sub Panel2 {
  print "<p class=\"clear\"></p>\n";
  print "<div class=\"menu_heading_open\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>Select reference organism, query genes and taxon</b> </div>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible_display\">\n";

  print "<p/><fieldset>\n";
  print "<legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Select reference organism, query genes and taxon</a></b></legend>\n";

  &PrintOrthoSelectionSection();

### use predicted leader genes
  print "<br>";
  print $query->checkbox(-name=>'leaders',
			 -checked=>$default{leaders},
			 -label=>'');
  print "<A class='iframe' HREF='help.footprint-scan.html#leader'><B>\n";
  print "predict operon leader genes";
  print "</B></A>\n";

  print "<br/>";

  print "</fieldset><p/>";
  print '</div>
</div>
<p class="clear"></p>';
}

##########################################
sub Panel3 {
################################################################
## Background model
print "<hr>";

my %bg_params =("markov" => 1,
		"bg_input" => 1,
		"bg_window" => 1,
		"markov_message" => 1
    );
&GetBackgroundModel(%bg_params);

print '
<div class="menu_heading_closed" onclick="toggleMenu(\'102\')" id="heading102"><b>Scanning Parameters</b> </div>
<div id="menu102" class="menu_collapsible">';



print "<p/><fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Scanning Parameters </a></b></legend>";


  ## Occurrences
  my $thresh_occ =
    $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
		  $query->Tr({-align=>center,-valign=>MIDDLE},
			     [
			      $query->th(["<A HREF='help.footprint-scan.html#return_fields'>Field</A> ",
					  " <A HREF='help.footprint-scan.html#thresholds'>Threshold</A> "]),

			      ### Threshold on score
			      $query->td(['Max site Pvalue',
			      		  $query->textfield(-name=>'uth_pvalue',
							    -default=>$default{uth_pvalue},
							    -size=>5)
					 ]),

				### Threshold on occ_inv_cum
			      $query->td(['Min number of sites',
					  $query->textfield(-name=>'uth_occ_th',
							    -default=>$default{uth_occ_th},
							    -size=>5)
					 ]),

			     ]
			    )
		 );
print "<br/> <b>Threshold </b> <br/>";

print "<td bgcolor='#F6E6CA'>$thresh_occ</td>";


print "</p>";


print "</fieldset><p/>";

print '</div>';

}

################################################################
## Comparisons with motif databases
sub Panel4 {
  print '
<br/>
<div>';
  print '
<div class="menu_heading_closed" onclick="toggleMenu(\'103\')" id="heading103">
<b>Drawing parameters</b></div>
<div id="menu103" class="menu_collapsible">';

  ## Tasks
  print "<fieldset><legend><b><a href='help.peak-motifs.html#tasks'>Drawing parameters</a></b></legend>";

print "<br>";
## Image format
print "Format ";
print $query->popup_menu(-name=>'img_format',
			 -Values=>["jpg","png","eps","pdf"],
			 -default=>$default{img_format});
print "<br>";
## draw lines to join points
print $query->checkbox(-name=>'info_lines',
		       -checked=>$default{info_lines},
		       -label=>' Informative lines');


  print "<p/> ";


  print "</fieldset><p/>";

  print '
</div>
</div>
<p class="clear"></p>';
}
