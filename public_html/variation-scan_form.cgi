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
require "matrix_web_forms.lib.pl";
use RSAT::MatrixReader;
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
## Default values for filling the form

## matrix-scan
$default{demo_descr1} = "";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "transfac";
$default{pseudo_distribution} = "pseudo_prior";
$checked{$default{pseudo_distribution}} = "CHECKED";
$default{mml}=30;

#$default{compare_motif_database}="jaspar_core_nonredundant_vertebrates";
#$default{db_choice_set}="jaspar_core_nonredundant_vertebrates";
#$default{custom_motif_db_name}="custom_motif_collection";

## Threshold values for site detection
## Suported uth 
$default{uth_pval} = "1e-3";

## Suported lth
$default{lth_score} = "1";
$default{lth_w_diff} = "1";
$default{lth_pval_ratio} = "10";


## Background model
$default{markov_order} = "1";

$default{leaders} = 'checked';
$default{bg_method}="bgfile";
$checked{$default{bg_method}} = "CHECKED";
$default{organism}="Homo_sapiens_GRCh37";
#$default{organism}="Homo_sapiens";
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
print "Scan sequences bearing variants with a list of motifs to predict which motifs may be affected by the variant.<br/>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>:</br> ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.epernicus.com/am27'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", <a target='_blank' href=''>Jeremy Delerce</a><sup>ci</sup>\n";
print "</CENTER>";
print "</BLOCKQUOTE>\n";

################################################################
## Display the form only if this RSAT instance supports variation
## analysis.
&check_variation_tools();

################################################################
### formheader

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

## Demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"variation-scan.cgi");

################# Matrix input
 &Panel1();


################# Input variant sequences
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

my $descr1 = "<H4>Comment on the demonstration :</H4>\n";
$descr1 .= "<blockquote class ='demo'>";
$descr1 .= "<p>In this demonstration, we use <i>variation-scan<\i> to assess the effect that a genetic variants have on transcription factor binding.</p>\n
<p> The genetic variants used in this example were collected by Weirauch, et al (Cell, 2014), these variants were reported in previous publications as affecting transcription factor binding. Motifs correspond to the transcription factores which biniding was reported to be affected by Weirauch, et al.</p>\n";
$descr1 .= "</blockquote>";

print $query->start_multipart_form(-action=>"variation-scan_form.cgi");
print $query->hidden(-name=>'queries',-default=>$demo_queries);
print $query->hidden(-name=>'organism',-default=>"Homo_sapiens_GRCh37");

#$demo_matrix_file=$ENV{rsat_www}."/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_TFs.tf";
$demo_matrix=`cat demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_TFs.tf`;
$demo_var_seq=`cat demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.varseq`;

print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'transfac');
print $query->hidden(-name=>'variants_seqs', -default=>$demo_var_seq);

#print $query->hidden(-name=>'db_choice', -default=>'custom_motif_db_url' );
#print $query->hidden(-name=>'compare_motif_database', -default=>'custom_motif_db' );
#print $query->hidden(-name=>'custom_motif_db_url', -default=>'on' );
#print $query->hidden(-name=>'',-default=>'on');
#print $query->hidden(-name=>'db_choice', -default=>'jaspar_pbm_mouse' );
#print $query->hidden(-name=>'db_choice' , -default=>'custom_motif_db' ." checked" );
#print $query->hidden(-name=>'jaspar_core_fungi' , -default=>'on'  );

print $query->hidden(-name=>'db_choice', -default=>'custom_motif_text');
print $query->hidden(-name=>'matrix', -default=>$demo_matrix);
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

# sub Panel1 {

#   print "<fieldset>\n<legend><b><a href='help.formats.html'>Matrix </a></b></legend>\n";

#   &GetMatrix();

#   print '<b><font style="font-size:80%"><a href="tutorials/tut_peak-motifs.html#seq" target="_blank"></a></font></b></br>';
#   print "</fieldset><p/>";
# }

################################################################
## Select motif data base to be used to scan variants

sub Panel1 {
  print '
<br/>';
#<div>' #  print "
# <div class=\"menu_heading_closed\" onclick=\"toggleMenu('103')\" id=\"heading103\">

#  <span title=\"A set conformed by different motifs representing binding profiles for several transcription factors will be used to scan the variants. Select here the known motifs collection or provide your own.\"></span> </div>\n
# <div id=\"menu103\" class=\"menu_collapsible\">";
  
  ## Tasks
  print "<fieldset><legend><a href='help.variation-scan.html'><b>Matrix collections</b></a></legend>";
  

  print "<p/> ";
  print "<b>A set conformed by different motifs representing binding profiles for several transcription factors will be used to scan the variants. Provide your own motifs collection or select one to be used.</b><br/>";
   print "<BR>\n";
  
  
  ## Old code to input the motif database via file or url
  # print "<input type='radio' NAME='db_choice' VALUE='custom_motif_db' $checked{file_upload}>";
  # print "<a href=''><b>Use your own motif database file:</b></a><br/>";
  # print $query->filefield(-name=>'custom_motif_db_txt',
  # 			  -size=>10);
  # print "<p/> ";
  # print "<input type='radio' NAME='db_choice' VALUE='custom_motif_db_url' >";   
  # print "<a href=''><b>Use your own motif database from URL source:</b></a><br/>";
  # print $query->textfield(-name=>'custom_motif_db_url_txt',
  # 		  -default=>$default{'custom_motif_db_url'},
  # 						-size=>62);

  my %matrix_args= (
      'db_choice'=>1,
      'no_pseudo'=>1,
      'status_db_choice'=>"checked"
      );
  
  &GetMatrix(%matrix_args);

  
  ## load the various databases that can be compared against
  print "<p/>";
  print "<b>Select one motif collection</b></p>";
  &DisplayMatrixDBchoice("mode"=>"radio");
  print "<p/> ";

 
  print "</fieldset><p/>";

  print '
</div>
</div>
<p class="clear"></p>';
}

##########################################
sub Panel2 {
  print "<p class=\"clear\"></p>\n";
 # print "<div class=\"menu_heading_open\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>Variant Sequence</b> </div>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible_display\">\n";

  print "<p/><fieldset>\n";
  print "<legend><b><a href='help.variation-scan.html'>Input variation sequences</a></b></legend>\n";


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
      $variants_seqChoiceString .=  "<a href=".$variants_seq_url.">";
      $variants_seqChoiceString .=  " transferred from previous query<BR>\n";
      $variants_seqChoiceString .=  "</a>";   
      $variants_seqChoiceString .=  "<INPUT type='hidden' NAME='variants_seq_file' VALUE='".$variants_seq_file."'>\n";
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
   
################################################################
## Background model
    
    #my %bg_params =("markov" => 1,
#		    #"bg_input" => 1,
#		    "no_bg_pseudo" => 1,
#		    "markov_message" => 1
		    #"ensembl"=>1
#	);
    
print "<fieldset>
<legend><b><a href='help.convert-matrix.html#io_format'>Background </a></b></legend>";

my %bg_params =("markov" => 1,
		"bg_input" => 0,
		"markov_message" => 0,
		"no_bg_pseudo"=>1
	       );
&GetBackgroundModel(%bg_params);

print "</fieldset><p/>";
   
    
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

			      ### Threshold on weight score
			      $query->td(['Weight of predicted sites',
					  $query->textfield(-name=>'lth_score',
							    -default=>$default{lth_score},
							    -size=>5),
					  ""
					 ]),
			      
			      ### Threshold on Sig of the score
			      $query->td(['Weight difference between variants',
					  $query->textfield(-name=>'lth_w_diff',
							    -default=>$default{lth_w_diff},
							    -size=>5),
					  ""
					 ]),

                              ### Threshold on P-value of the score
			      $query->td(['P-value of predicted sites',
					  "",
					  $query->textfield(-name=>'uth_pval',
							    -default=>$default{uth_pval},
							    -size=>5)				       
					 ]),

			      ### Threshold on P-value of the score
			      $query->td(['P-value ratio between variants',
					  $query->textfield(-name=>'lth_pval_ratio',
							    -default=>$default{lth_pval_ratio},
							    -size=>5),
					  ""
					 ]),
					  
		      
			     ]
			    )
		 );

   
    print "<br> <b>Thresholds </b> </br>";
    
    print "<td bgcolor='#F6E6CA'>$thresh_matches</td>";
 
    
    print "</p>";
    
    
    print "</fieldset><p/>";
    
    print '
</div></div>';
    
}

