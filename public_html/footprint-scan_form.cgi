#!/usr/bin/perl
################################################################
## this cgi script fills the HTML form for the program footprint-scan
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

# BEGIN {
#     if ($0 =~ /([^(\/)]+)$/) {
# 	push @INC, "$`lib/";
#     }
# }
# use CGI;
# use CGI::Carp qw/fatalsToBrowser/;
# require "RSA.lib";
# require "RSA2.cgi.lib";
# use RSAT::Tree;


# $ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

local @supported_input_formats = sort(keys( %RSAT::MatrixReader::supported_input_format));
local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));
################################################################
## Initialize parameters


# ## Output fields
# my @output_fields = qw(ident
# 			   ali_len
# 			   mismat
# 			   gap_open
# 			   e_value
# 			   bit_sc
# 			   rank);
# my %field_description = ();
# $field_description{ident} = "Percentage of identity";
# $field_description{ali_len} = "Alignment length";
# $field_description{mismat} = "Number of mismatches";
# $field_description{gap_open} = "Number of gap openings";
# $field_description{e_value} = "E-value";
# $field_description{bit_sc} = "Bit score";
# $field_description{rank} = "Rank";
# #$field_description{s_rank} = "target rank";


################################################################
### default values for filling the form

## Default values for dyad-analysis

$default{demo_descr}="";
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "meme";
$default{pseudo_prior} = "pseudo_prior";
$default{pseudo_counts}="1";
$checked{'tf'}="";
$default{tf}="";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{bg_method}="bginput";
$checked{$default{bg_method}} = "CHECKED";



$default{organism}="Escherichia_coli_K12";
$default{background} = "upstream-noorf";
$default{bg_level} = "organism";
$default{markov_order} = "1";
$default{bg_method2}="background";
$default{filter_pval}="1e-4";
## Default parameters for get-orthologs

## Threshold values for occurrence statistics
$default{lth_occ_score} = "none";
$default{uth_occ_score} = "none";
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
$default{uth_occ_sig_rank} = "none";

&LoadGetOrthoDefault(\%default);

## Other default parameters


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
&RSA_header("footprint-scan", "form");
print "<CENTER>";
print "Given one or several genes from a query organism and a PSSM, collect all orthologous genes for a given taxonomical level, <br> scan the promotes and  detect conserved sites.<br>\n";
print "(Program developed by <A href='mailto:amedina\@lcg.unam.mx'>Alejandra Medina-Rivera</A>\n";
print "and <a href='mailto:jvanheld\@bigre.ulb.ac.be'>Jacques van Helden</A>).\n";
#print "<br>Reference: <a target='_blank' href=\"http://www.biomedcentral.com/1471-2105/9/37\">Janky & van Helden, BMC Bioinformatics 2008, 9:37.</a>";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $default{demo_descr};
&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->start_multipart_form(-action=>"footprint-discovery.cgi");

################################################################
## Matrix
print "<hr>";
print "<fieldset>
<legend><b><a href='help.convert-matrix.html#io_format'>1 - Matrix </a></b></legend>";
&GetMatrix();
## TF
    print ("<br><INPUT TYPE='radio' NAME='tf' VALUE='tf' $checked{'tf'}>","<b>TF gene name</b> &nbsp;");
    print $query->textfield(-name=>'tf',
			    -default=>$default{tf},
			    -size=>5);
print "</fieldset><p/>";
################################################################
## Print the options for the selection of orthologs
print "<hr>";
print "<fieldset>
<legend><b><a href=''>2 - Organism & Taxon </a></b></legend>";

&PrintOrthoSelectionSection();

### use predicted leader genes
print "<br>";
print $query->checkbox(-name=>'leaders',
		       -checked=>$default{leaders},
		       -label=>'');
print "<A HREF='help.footprint-discovery.html#leader'><B>\n";
print "predict operon leader genes";
print "</B></A>\n";
print "</fieldset><p/>";
################################################################
#### Options for scanning
print "<hr>";
print "<p><b>Options for <i>scanning</i></b></p>";

#print "<ul>";

print "<fieldset>
<legend><b><a href='help.matrix-scan.html#markov_order'>Background </a></b></legend>";
my %bg_params =("markov" => 1,
		"bg_input" => 1,
		"bg_window" => 1,
		"markov_message" => 1,
		"taxon" => 1
	       );
&GetBackgroundModel(\%bg_params);

#&PrintScanBackgroundOptions(1);

print "<br/>Note: Only Bernoulli models are supported. Higher-order Markov models are converted into Markov 0 (Bernoulli).";
print "</fieldset><p/>";



################################################################
## Print dyad return fields
print "<fieldset>
<legend><b><a href=''>Thresholds </a></b></legend>";
&ReturnTable ;
print "</fieldset><p/>";




################################################################
#### Filtering
print "<hr>";
print "<p><b>Options for <i>filtering</i></b></p>";

#print "<ul>";

### filtering scandyads
print $query->checkbox(-name=>'scan_filter',
		       -checked=>$default{dyads_filter},
		       -label=>'');
print "<a href='help.footprint-scan.html#filtering'><b>\n";
print "Scan filtering</b></a>\n";
print "(only accept genes having at least one site occurrence in the promoter of the query gene)";
   #### Background specifiaction

#### filter threshold
    print "<br>\n";
    print "<A HREF='help.patser.html#uthreshold'><B>Filter Pval</B></A>", "&nbsp"x6;
    print $query->textfield(-name=>'filter_pval',
			    -default=>$default{filter_pval},
			    -size=>6);
print "<fieldset>
<legend><b><a href='help.matrix-scan.html#markov_order'>Background </a></b></legend>";
my %bg_params =(
    "markov" => 1,
    "markov_message" => 1
    );

&PrintScanBackgroundOptions(2);

print "<br/>Note: Only Bernoulli models are supported. Higher-order Markov models are converted into Markov 0 (Bernoulli).";
print "</fieldset><p/>";


################################################################
### send results by email or display on the browser
print "<hr>";
&SelectOutput('email', email_only=>1);

################################################################
### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"footprint-scan_form.cgi");
$demo_queries = "lexA\n";
#$demo_queries .= "recA\n";
#$demo_queries .= "uvrB\n";
print "<TD><B>";
print $query->hidden(-name=>'queries',-default=>$demo_queries);
print $query->hidden(-name=>'organism',-default=>"Escherichia_coli_K12");
print $query->hidden(-name=>'taxon',-default=>"Enterobacteriales");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.footprint-discovery.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_footprint-discovery.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</BLOCKQUOTE>\n";

print "<hr class=solid>";

print $query->end_html;

exit(0);

################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {

 #  my $thresh_matches =
#     $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
# 		  $query->Tr({-align=>center,-valign=>MIDDLE},
# 			     [
# 			      $query->th([" <A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
# 					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
# 					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A>"]),

# 			      ### Threshold on score
# 			      $query->td(['Weight<br>score',
# 					  $query->textfield(-name=>'lth_score',
# 							    -default=>$default{lth_score},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_score',
# 							    -default=>$default{uth_score},
# 							    -size=>5)
# 					 ]),

# 			      ### Threshold on P-value of the score
# 			      $query->td(['P-value',
# 					  $query->textfield(-name=>'lth_pval',
# 							    -default=>$default{lth_pval},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_pval',
# 							    -default=>$default{uth_pval},
# 							    -size=>5)
# 					 ]),
					  
# 				### Threshold on Sig of the score
# 			      $query->td(['Sig',
# 					  $query->textfield(-name=>'lth_sig',
# 							    -default=>$default{lth_sig},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_sig',
# 							    -default=>$default{uth_sig},
# 							    -size=>5)
# 					 ]),

# 			      ### Threshold on proba_M
# 			      $query->td(['P(S|M)<br>proba_M',
# 					  $query->textfield(-name=>'lth_proba_M',
# 							    -default=>$default{lth_proba_M},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_proba_M',
# 							    -default=>$default{uth_proba_M},
# 							    -size=>5)
# 					 ]),

# 			      ### Threshold on proba_B
# 			      $query->td(['P(S|B)<br>proba_B',
# 					  $query->textfield(-name=>'lth_proba_B',
# 							    -default=>$default{lth_proba_B},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_proba_B',
# 							    -default=>$default{uth_proba_B},
# 							    -size=>5)
# 					 ]),

# 			      ### Threshold on normw
# 			      $query->td(['Normalized<br>weight',
# 					  $query->textfield(-name=>'lth_normw',
# 							    -default=>$default{lth_normw},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_normw',
# 							    -default=>$default{uth_normw},
# 							    -size=>5)
# 					 ]),

# 			      ### Threshold on rank
# 			      $query->td(['Rank',
# 					  $query->textfield(-name=>'lth_rank',
# 							    -default=>$default{lth_rank},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_rank',
# 							    -default=>$default{uth_rank},
# 							    -size=>5)
# 					 ]),

# 			     ]
# 			    )
# 		 );

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

			      $query->th({-colspan=>3,-align=>left},["Over-representation"
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
#   my $thresh_crer =
#     $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
# 		  $query->Tr({-align=>center,-valign=>MIDDLE},
# 			     [
# 			      $query->th([" <A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
# 					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
# 					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A> "]),


# 			      ### Threshold on score
# 			      $query->td(['CRER size<b>*</b>',
# 					  $query->textfield(-name=>'lth_crer_size',
# 							    -default=>$default{lth_crer_size},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_crer_size',
# 							    -default=>$default{uth_crer_size},
# 							    -size=>5)
# 					 ]),

# 				### Threshold on P-value of the score
# 			      $query->td(['site P-value<b>*</b>',
# 					  $query->textfield(-name=>'lth_site_pval',
# 							    -default=>$default{lth_site_pval},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_site_pval',
# 							    -default=>$default{uth_site_pval},
# 							    -size=>5)
# 					 ]),

# 				### Threshold on crer_sites
# 			      $query->td(['CRER sites',
# 					  $query->textfield(-name=>'lth_crer_sites',
# 							    -default=>$default{lth_crer_sites},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_crer_sites',
# 							    -default=>$default{uth_crer_sites},
# 							    -size=>5)
# 					 ]),

# 				### Threshold on crer_pval
# 			      $query->td(['CRER pval',
# 					  $query->textfield(-name=>'lth_crer_pval',
# 							    -default=>$default{lth_crer_pval},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_crer_pval',
# 							    -default=>$default{uth_crer_pval},
# 							    -size=>5)
# 					 ]),
# 			      ### Threshold on crer_pval
# 			      $query->td(['CRER sig',
# 					  $query->textfield(-name=>'lth_crer_sig',
# 							    -default=>$default{lth_crer_sig},
# 							    -size=>5),
# 					  $query->textfield(-name=>'uth_crer_sig',
# 							    -default=>$default{uth_crer_sig},
# 							    -size=>5)
# 					 ]),
# 			      $query->Tr({-align=>middle,-valign=>TOP},
# 					 [
# 					  $query->td({-colspan=>4},[ "<b>*</b> =mandatory field"],
# 						    )
# 					 ]),
# 			     ]
# 			    )
# 		 );


  #############################################
  ## Table

  print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>3},
		      "<tr>",
		      "<th bgcolor='#D6EEFA'>
  			Thresholds</th> ",
		      "</tr>",		      		   
		      "<td bgcolor='#D6EEFA'>$thresh_occ</td>",
		      "</tr>",
		     );

}

sub PrintScanBackgroundOptions {
   my $bgtag=shift(@_);
    $checked{$default{'bg_method_'.$bgtag}} = "CHECKED";

    print "<a href='help.oligo-analysis.html#exp_freq'><B>Background model</b></a>&nbsp;";

    &PrintGenomeSubsetBgOptionsFS($bgtag);

  #### Estimated from the input sequence set

  print ("<br><b>Estimate from input sequence</b>");

  ## Markov model
  #$Freq_name="bg_method_".$bgtag;
  print ("<br><input type='radio' NAME='bg_method' VALUE='Markov model (higher order dependencies)' $checked{'Markov model (higher order dependencies)'}>", "Markov model (higher order dependencies)");
   # print ("<br><input type='radio' NAME='bg_method' VALUE='Markov model (higher order dependencies)' CHECKED >", "Markov model (higher order dependencies)");

  print "&nbsp;&nbsp;order &nbsp;";
  print $query->textfield(-name=>'markov_order',
			  -default=>$default{markov_order},
			  -size=>5);

  #### Lexicon partitioning
  #print "<br><input type='radio' NAME='bg_method' VALUE='Lexicon partitioning' $checked{'Lexicon partitioning'}>Lexicon partitioning<p>";

  #### Bernouilli model
  #print "<br><input type='radio' NAME='bg_method' VALUE='Residue frequencies from input sequence' $checked{'Residue frequencies from input sequence'}>Residue frequencies from input sequence<p>";

  #### equiprobable residues
  print "<br><input type='radio' NAME='bg_method' VALUE='Equiprobable residues' $checked{'Equiprobable residues'}>Equiprobable residues (<A HREF='help.oligo-analysis.html#equiprobable'>usually NOT recommended</a>)<p>";

  #### custom expected frequency file
  print ("<a href='help.oligo-analysis.html#upload_freq_file'><b>Upload your own expected frequency file</b></a><BR>");
  print ("<br><input type='radio' NAME='$Freq_name' VALUE='file_upload' $checked{'file_upload'}>");
  print ($query->filefield(-name=>'upload_freq_file',-default=>'starting value',-size=>30,-maxlength=>200));
  print "<p>";
}

sub PrintGenomeSubsetBgOptionsFS {
    my $bgtag=shift(@_);
  $checked{$default{bg_level}} = "CHECKED";
  $checked{$default{'bg_method'.$bgtag}} = "CHECKED";

  #### Calibrated on genome subsets

  print( "<br><input type='radio' NAME='bg_method' VALUE='background' $checked{background} >");

# print( "<br><input type='radio' NAME='bg_method' VALUE='background'  >");
  print ("<b>Genome subset</b>");
#  print "<ul>";

  print ( "&nbsp;&nbsp;<a href='help.oligo-analysis.html#background'>Sequence type</a> &nbsp;&nbsp;&nbsp;&nbsp;",$query->popup_menu(-name=>'background',-Values=>["upstream","upstream-noorf","protein"],-default=>$default{background}));


  print ("<ul>");
  print( "<input type='radio' NAME='bg_level' VALUE='organism' $checked{organism}>", &OrganismPopUpString());
  print( "<br><input type='radio' NAME='bg_level' VALUE='taxon' $checked{taxon}>", &TaxonPopUpString("node"));
  print ("</ul>");
#  print "</ul>";


  print "<p>";
}
