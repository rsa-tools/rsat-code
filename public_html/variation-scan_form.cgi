#!/usr/bin/perl
#### this cgi script fills the HTML form for the program varition-scan
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
$default{mml}=30;



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

## Background model
$default{markov_order} = "1";

$default{leaders} = 'checked';
$default{bg_method}="bgfile";
$checked{$default{bg_method}} = "CHECKED";
$default{organism}="Homo_sapiens_GRCh37";
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
&RSA_header("variation-scan", "form");
print "<CENTER>";
print " variation-scan.<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.epernicus.com/am27'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", <a target='_blank' href=''>Jeremy Delerce</a><sup>ci</sup>\n";
print "</CENTER>";
print "</BLOCKQUOTE>\n";
print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'105\')\" id=\"heading105\"><font color='#0D73A7'>Information about variation-scan</font> </div>\n";
 print "<div id=\"menu105\" class=\"menu_collapsible\">\n";
print "<BLOCKQUOTE>\n";
print "<p>Scan variant sequences with position specific scoring matrices (PSSM)
    and report variations that affect the binding score, in order to predict
    regulatory variants.</p> ";
print "</BLOCKQUOTE>\n";
print "</div></p>\n";

&ListParameters() if ($ENV{rsat_echo} >=2);

#&ReadMatrixFromFile() ;

## demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"variation-scan.cgi");

################# Matrix input
 &Panel1();


################# Select reference organism, taxon and query genes
 &Panel2();

################# Scanning Parameters
 &Panel3();

################# Drawing parameters
##&Panel4();



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

################################################################
### data for the demo

my $descr1 = "<H4>Comment on the demonstration example 1 :</H4>\n";
$descr1 .= "<blockquote class ='demo'>";

$descr1 .= "<p>In this demonstration, we use <i>variation-scan<\i> to assess the effect that a genetic variants have on transcription factor binding.</p>\n

<p> The genetic variants used in this example have been previously reported in the literature. Used motifs correspond to TFs which biniding is expected to be affected</p>\n";

$descr1 .= "</blockquote>";

print $query->start_multipart_form(-action=>"variation-scan_form.cgi");
print $query->hidden(-name=>'queries',-default=>$demo_queries);
print $query->hidden(-name=>'organism',-default=>"Homo_sapiens_GRCh37");


$demo_matrix=`cat demo_files/do798+do735_mmus_hnf6_liver.transfac`;
$demo_var_seq="./demo_files/variation_demo_set_rsat_var.seq";

print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'variants_seq_file', -default=>$demo_var_seq);

print $query->hidden(-name=>'bg_method',-default=>'background"');
print $query->hidden(-name=>'background',-default=>'upstream-noorf');
print $query->hidden(-name=>'markov_order',-default=>'2');

print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


##print "<td><b><a href='tutorials/tut_peak_motif.html'>[TUTORIAL]</a></B></TD>\n";
print "<td><b><a href='help.variation-scan.html'>[MANUAL]</a></B></TD>\n";
#print "<td><b><a href='tutorials/tut_peak-motifs.html'>[TUTORIAL]</a></B></TD>\n";
print "<TD><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



#################### Internal functions #######################

sub Panel1 {

  print "<fieldset>\n<legend><b><a href='help.formats.html'>Matrix </a></b></legend>\n";

  print "<table>
  <tr><td colspan='2' style='text-align:center;'> ";
  &GetMatrix();


  print "</td>
  </tr>";

  print "</td></tr></table>";
  print '<b><font style="font-size:80%"><a href="tutorials/tut_peak-motifs.html#seq" target="_blank"></a></font></b></br>';
  print "</fieldset><p/>";
}

##########################################
sub Panel2 {
  print "<p class=\"clear\"></p>\n";
  print "<div class=\"menu_heading_open\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>Variant Sequence</b> </div>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible_display\">\n";

  print "<p/><fieldset>\n";
  print "<legend><b><a href='help.peak-motifs.html#tasks'>Select reference organism and input variations</a></b></legend>\n";


### Input variant-seqs
  print "<br>";
  print "<A HREF='help.variation-scan.html#input'><B>\n";
  print "Variant sequences</p>";
  print "</B></A>\n";
  if ($variants_seq_file = $query->param('variants_seq_file')) {
      ## Variants_Seq file is already on the server machine
      ## (piped from a previous script)
      $variants_seq_url = $variants_seq_file;
      $variants_seq_url =~ s|$ENV{RSAT}/public_html|$ENV{rsat_www}|;
      $variants_seqChoiceString .=  "<a href=$variants_seq_url>";
      $variants_seqChoiceString .=  " transferred from previous query<BR>\n";
      $variants_seqChoiceString .=  "</a>";   
      print $variants_seqChoiceString ;

  } else {
      print $query->textarea(-name=>'variants_seqs',
			     -default=>"",
			     -rows=>6,
			     -columns=>65);
      
      
      
      print "<br/>";
      print "<BR>Upload variant sequences<BR>\n";
      print $query->filefield(-name=>'uploaded_file',
			      -default=>'',
			      -size=>45,
			      -maxlength=>200);
      print "</UL>\n";
      print "<BR>\n";
  }
      print "</fieldset><p/>";
  print '</div>
</div>
<p class="clear"></p>';
}

##########################################
sub Panel3 {
    print '
<div class="menu_heading_closed" onclick="toggleMenu(\'102\')" id="heading102"><b>Scanning Parameters</b> </div>
<div id="menu102" class="menu_collapsible">';
    
    print "<p/><fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Scanning Parameters </a></b></legend>";
## Lenght of the sequences surranding the variant
    print "<B>Length of sequence around the variant</B>&nbsp;\n";
    print $query->textfield(-name=>'mml',
			    -default=>$default{mml},
			    -size=>5);
    
    ## Threshold table

  my $thresh_matches =
    $query->table({-border=>0,-cellpadding=>1,-cellspacing=>0},
		  $query->Tr({-align=>center,-valign=>MIDDLE},
			     [
			      $query->th([" <A HREF='help.matrix-scan.html#return_fields'>Field</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Lower<BR>Threshold</A> ",
					  " <A HREF='help.matrix-scan.html#thresholds'>Upper<BR>Threshold</A>"]),

			      ### Threshold on score
			      $query->td(['Weight<br>score',
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

   
    print "<br/> <b>Threshold </b> <br/>";
    
    print "<td bgcolor='#F6E6CA'>$thresh_matches</td>";
    
################################################################
## Background model
    print "<hr>";
    
    my %bg_params =("markov" => 1,
		    #"bg_input" => 1,
		    "no_bg_pseudo" => 1,
		    "markov_message" => 1
	);
    &GetBackgroundModel(%bg_params);
    
    
    
    
    print "</p>";
    
    
    print "</fieldset><p/>";
    
    print '
</div></div>';
    
}

