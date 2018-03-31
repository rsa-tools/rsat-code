#!/usr/bin/perl
#### this cgi script fills the HTML form for the program matrix-scan
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
# if ($0 =~ /([^(\/)]+)$/) {
#     push (@INC, "$`lib/");
# }
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "matrix_web_forms.lib.pl";
use RSAT::MatrixReader;
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
$default{markov_order} = "1";
$default{organism} = "";
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
print "<b>Citation</b>: <a href='https://www.linkedin.com/in/jean-valery-turatsinze-23a6aa41/'>Jean Val&eacute;ry Turatsinze</A>, <A HREF='http://morgane.bardiaux.fr/'>Morgane Thomas-Chollier</A>, <a href='https://www.ulb.ac.be/rech/inventaire/chercheurs/4/CH10464.html'>Matthieu Defrance</a> and <A HREF='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a> (2008). Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nat Protoc, 3, 1578-1588. <a href='http://www.ncbi.nlm.nih.gov/pubmed/18802439'>Pubmed 18802439</a>";


print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"matrix-scan.cgi", -id=>"form", -onreset=>"resetHandler()");

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
print "<A class='iframe' HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['single',
				   'both'],
			 -default=>$default{strands});

################################################################
#### origin for calculating positions
print "&nbsp;"x4,  "<A class='iframe' HREF='help.matrix-scan.html#origin'><B>Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',-id=>'origin',
			 -Values=>['start',
				   'center',
				   'end',
				   'genomic'],
			 -default=>$default{origin});

################################################################
#### Offset for calculating positions
print "&nbsp;"x4,  "<A class='iframe' HREF='help.matrix-scan.html#offset'><B>Offset</B></A>\n";
print $query->textfield(-name=>'offset',
			-default=>$default{offset},
			-size=>8);

################################################################
#### decimals
print "&nbsp;"x2,  "<A class='iframe' HREF='help.matrix-scan.html#decimals'><B>score decimals</B></A>\n";
print $query->popup_menu(-name=>'decimals',
			 -Values=>['0',
				   '1','2'],
			 -default=>$default{decimals});

################################################################
#### decimals
print "&nbsp;"x2,  "<A class='iframe' HREF='help.matrix-scan.html#skipscore'><B>handling of N characters</B></A>\n";
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
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
### data for the demo

## Load demo sequences
$demo_sequence_file = $ENV{RSAT}."/public_html/demo_files/Dmelanogaster_eve_up5000.fasta";
$demo_sequence = "";
open($fh, $demo_sequence_file);
while(my $row = <$fh>){
    chomp $row;
    if($row =~ /\'/){
    $row =~ s/\'/\\'/g;
    }
    $demo_sequence .= $row;
    $demo_sequence .= "\\n";
}

## Load demo matrices
$demo_matrix_file = $ENV{RSAT}."/public_html/demo_files/Dmelanogaster_segmentation_12matrices.tf";
$demo_matrix = "";
open($fh, $demo_matrix_file);
while(my $row = <$fh>){
    chomp $row;
    $demo_matrix .= $row;
    $demo_matrix .= "\\n";
}

## demo 1

print '<script>

descr = "<H4>Comment on the demonstration example : </H4><blockquote class =\'demo\'>In this demonstration, we will analyse the promoter of \
<i>Drosophila melanogaster</i> even-skipped gene (<i>eve</i>). We will scan the 5500 \
bp sequence upstream the transcription start site with matrices \
representing the binding specificity of 12 transcription factors known \
to regulate <i>eve</i>. These matrices were built from binding sites \
annotated in the <a target=_blank \
href=\'http://www.oreganno.org\'>ORegAnno</a> database by Jean-Valery \
Turatsinze.<p/>";

function setDemo1(demo_matrix, demo_sequence){
    $("#reset").trigger("click");
    $("#db_choice").val("").change();
    descr_1 = descr + "The program will return individual matches, i.e. sequence segments scoring above the predefined threshold. In this example, threshold is set on the P-value.</blockquote>";
    
    demo_descr.innerHTML = descr_1;
    demo.value = descr_1;
    $("#bg_method_bginput").prop("checked",true);
    $("#uth_pval").val("1e-4");
    background.value = "upstream-noorf";
    markov_order.value = "1";
    
    $("#analysis_type_sites").prop("checked",true);
    $("#return_rank").prop("checked",false);
    matrix.value = demo_matrix;
    matrix_format.value = "transfac";
    $("#consensus_as_name").prop("checked", false);
    sequence.value = demo_sequence;
    origin.value = "genomic";
}

function setDemo2(demo_matrix, demo_sequence){
    $("#reset").trigger("click");
    $("#db_choice").val("").change();
    descr_2 = descr + "The program will return CRERs: regions of a few hundreds residues that have a higher density of matches than expected by chance.</blockquote>";
    
    demo_descr.innerHTML = descr_2;
    demo.value = descr_2;
    $("#bg_method_bgfile").prop("checked",true);
    $("#uth_site_pval").val("1e-4");
    background.value = "upstream-noorf";
    markov_order.value = "1";
    $("#organism_bg_name").val("Drosophila melanogaster")
    $("#organism_bg").val("Drosophila_melanogaster");
    $("#analysis_type_crer").prop("checked",true);
    $("#return_rank").prop("checked",false);
    matrix.value = demo_matrix;
    matrix_format.value = "transfac";
    $("#consensus_as_name").prop("checked", false);
    $("#crer_ids").prop("checked", false);
    sequence.value = demo_sequence;
    origin.value = "genomic";
}

function setDemo3(demo_matrix, demo_sequence){
    $("#reset").trigger("click");
    $("#db_choice").val("").change();
    descr_3 = descr + "The program will return matrices for which the total number of hits in the input sequences is higher than expected by chance.</blockquote>";
    
    demo_descr.innerHTML = descr_3;
    demo.value = descr_3;
    $("#bg_method_input").prop("checked",true);
    $("#uth_site_pval").val("1e-4");
    background.value = "upstream-noorf";
    markov_order.value = "1";
   
    $("#analysis_type_occ").prop("checked",true);
    $("#return_rank").prop("checked",false);
    matrix.value = demo_matrix;
    matrix_format.value = "transfac";
    $("#consensus_as_name").prop("checked", false);
    sequence.value = demo_sequence;
    origin.value = "genomic";
    $("#uth_occ_sig_rank").val("1");
    $("#lth_occ_score").val("5");
}

function resetHandler(){
    $("#db_choice").val("").change();
}

</script>';

print "<TD><B>";
print '<button type="button" onclick="setDemo1('. "'$demo_matrix'" .','."'$demo_sequence'".')">DEMO 1 (sites)</button>';
print "</B></TD>";

## demo2
print "<TD><B>";
print '<button type="button" onclick="setDemo2('. "'$demo_matrix'" .','."'$demo_sequence'".')">DEMO 2 (CRERs)</button>';
print "</B></TD>\n";


## demo3: detect enrichment of hits for PSSMs

print "<TD><B>";
print '<button type="button" onclick="setDemo3('. "'$demo_matrix'" .','."'$demo_sequence'".')">DEMO 3 (enrichment)</button>';
print "</B></TD>\n";

print "<TD><B><A class='iframe'  HREF='help.matrix-scan.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_matrix-scan.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
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
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				       -CHECKED=>$default{'return_'.$field},
				       -label=>' '.$field.' ');
  }
  $boxes_matches .= "<BR/>";
  @return_fields_matches = qw( site_limits normw);
  foreach my $field (@return_fields_matches) {
  		my $display_field = $field;
  		$display_field =~ s/site_//;
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				       -CHECKED=>$default{'return_'.$field},
				       -label=>' '.$display_field.' ');
  }
  $boxes_matches .= "<BR/>";
  @return_fields_matches = qw(weight_limits bg_residues);
  foreach my $field (@return_fields_matches) {
    $boxes_matches .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				       -CHECKED=>$default{'return_'.$field},
				       -label=>' '.$field.' ');
    $boxes_matches .= "<BR/>";
  }


  ### Return fields
  my $boxes_occ = "";
  @return_fields_occ = qw(distrib);
  foreach my $field (@return_fields_occ) {
    $boxes_occ .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				   -CHECKED=>$default{'return_'.$field},
				   -label=>' '.$field.' ');
  }
  $boxes_occ .= "<BR/>";
  @return_fields_occ = qw(occ_proba);
  foreach my $field (@return_fields_occ) {
    $boxes_occ .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				   -CHECKED=>$default{'return_'.$field},
				   -label=>' '.$field.' ');
  }

$boxes_occ .= "&nbsp;&nbsp; <b> sort by </b> &nbsp;&nbsp;".$query->popup_menu(-name=>'sort_distrib',-id=>'sort_distrib',
										-Values=>['scores',
											  'occ_sig'],
										-default=>$default{sort_distrib});

  ### Return fields
  my $boxes_crer = "";
  @return_fields_crer = qw(crer normw);
  foreach my $field (@return_fields_crer) {
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				    -CHECKED=>$default{'return_'.$field},
				    -label=>' '.$field.' ');
  }
  $boxes_crer .= "<BR/>";
  @return_fields_crer = qw(crer_limits);
  foreach my $field (@return_fields_crer) {
  	  	my $display_field = $field;
  		$display_field =~ s/crer_//;
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				    -CHECKED=>$default{'return_'.$field},
				    -label=>' '.$display_field.' ');
  }
  @return_fields_crer = qw(crer_sites);
  foreach my $field (@return_fields_crer) {
    $boxes_crer .= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
				    -CHECKED=>$default{'return_'.$field},
				    -label=>' '."sites".' ');
  }
  $boxes_crer .= "<BR/>";
  $boxes_crer .= $query->checkbox(-name=>'crer_ids',-id=>'crer_ids',
				  -CHECKED=>$default{crer_ids},
				  -label=>' '."crer-specific identifiers".' ');

  ### Return fields
  my $boxes_add = "";
  @return_fields_add = qw(matrix freq_matrix weight_matrix bg_model);
  foreach my $field (@return_fields_add) {
    $boxes_add.= $query->checkbox(-name=>'return_'.$field,-id=>'return_'.$field,
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
			      $query->th([" <A class='iframe' HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					  " <A class='iframe' HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					  " <A class='iframe' HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A>"]),

			      ### Threshold on score
			      $query->td(['Weight<br>score',
					  $query->textfield(-name=>'lth_score',-id=>'lth_score',
							    -default=>$default{lth_score},
							    -size=>5),
					  $query->textfield(-name=>'uth_score',-id=>'uth_score',
							    -default=>$default{uth_score},
							    -size=>5)
					 ]),

			      ### Threshold on P-value of the score
			      $query->td(['P-value',
					  $query->textfield(-name=>'lth_pval',-id=>'lth_pval',
							    -default=>$default{lth_pval},
							    -size=>5),
					  $query->textfield(-name=>'uth_pval',-id=>'uth_pval',
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
			      $query->th(["<A class='iframe' HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A> "]),

			      $query->th({-colspan=>3,-align=>left},["Occurrences"
								    ]),
			      ### Threshold on score
			      $query->td(['Weight<br>score',
					  $query->textfield(-name=>'lth_occ_score',-id=>'lth_occ_score',
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
					  $query->textfield(-name=>'uth_occ_sig_rank',-id=>'uth_occ_sig_rank',
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
			      $query->th([" <A class='iframe' HREF='help.matrix-scan.html#return_fields'>Field</A> ",
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
					  $query->textfield(-name=>'uth_site_pval',-id=>'uth_site_pval',
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
  			<INPUT TYPE='radio' NAME='analysis_type' id='analysis_type_sites' VALUE='analysis_sites' $checked{'analysis_sites'}><BR/>
  			<A class='iframe' HREF='help.matrix-scan.html#return_fields'>Individual matches</A></th>",
		      "<th bgcolor='#D6EEFA'>
  			<INPUT TYPE='radio' NAME='analysis_type' 'analysis_type_crer' VALUE='analysis_crer' $checked{analysis_crer}><BR/>
			<A class='iframe' HREF='help.matrix-scan.html#return_fields'>CRERs <BR/> (Cis-Regulatory element <BR>Enriched Regions)</A> </th>",
		      "<th bgcolor='#F6E6CA'>
  			<INPUT TYPE='radio' NAME='analysis_type' 'analysis_type_occ' VALUE='analysis_occ' $checked{analysis_occ}><BR/>
  			<A class='iframe' HREF='help.matrix-scan.html#return_fields'>Enrichment of hits<br>in the whole input sequence set</A></th> ",
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





