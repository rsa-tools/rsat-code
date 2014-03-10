#!/usr/bin/perl
#### this cgi script fills the HTML form for the program compare-matrices
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### default values for filling the form
$default{ref_db} = "CHECKED";
$default{demo_descr1} = "";
$default{matrix}="";
$default{matrix_file}="";
#$default{pseudo_counts}=1;
#$default{pseudo_distribution}="pseudo_prior";
#$checked{$default{pseudo_distribution}} = "CHECKED";
$default{matrix_format} = "transfac";
$default{bg_format}="oligo-analysis";
$default{bg_method}="from_matrix";
$checked{$default{bg_method}} = "CHECKED";

$default{'return_w'} = "CHECKED"; $default{'lth_w'} = 5;
$default{'return_cor'} = "CHECKED"; $default{'lth_cor'} = 0.7;
$default{'return_Ncor'} = "CHECKED"; $default{'lth_Ncor'} = 0.4;
$default{'return_match_rank'} = "CHECKED"; $default{'uth_match_rank'} = 50;
$default{'return_logoDP'} = "CHECKED";
$default{'return_NSW'} = "CHECKED";
$default{'return_NsEucl'} = "CHECKED";


## motif database
$default{compare_motif_database}="jaspar_core_vertebrates";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
#   if ($query->param($key) =~ /checked/i) {
#     $checked{$key} = "CHECKED";
#   }
}



&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
### print the form ###

################################################################
### header
&RSA_header("compare-matrices", "form");
print "<CENTER>";
print "Comparison between two collections of position-specific scoring matrices.<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</a><sup>t</sup>\n";
print ", <a target='_blank' href='http://www.ibens.ens.fr/spip.php?article26&lang=en'>Denis Thieffry</a><sup>t</sup>\n";
print "and <a target='_blank' href='http://biologie.univ-mrs.fr/view-data.php?id=202'>Carl Herrmann</a><sup>ct</sup>\n";
print "</CENTER>";

## demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"compare-matrices.cgi");


################################################################
#### Matrix specification
print "<hr>";

################################################################
## Query matrices
&GetMatrix('title'=>'Query matrices', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);
print "<hr>";

################################################################
## Database comparison
&DatabaseChoice();
print "<hr>";

################################################################
## Background model
my %bg_params =("from_matrix" => 1,
		"bg_input"=>0,
		"simple"=>1,
		"no_bg_pseudo"=>1,
	       );
&GetBackgroundModel(%bg_params);
print "<hr>";

################################################################
## Selection of output fields and thresholds
local @matching_scores = qw(w
			 cor
			 Ncor
                         logoDP
			 logocor
			 Nlogocor
			 Icor
			 NIcor
			 cov
			 dEucl
			 NdEucl
			 NsEucl
			 SSD
			 SW
			 NSW
			 match_rank
			 offset
			);

local %score_descriptions = ('w'=>'Width = number of aligned columns',
			     'cor'=>'Pearson correlation (computed on residue occurrences in aligned columns)',
			     'Ncor'=>'Relative width-normalized Pearson correlation',
			     'logoDP'=>'dot product of sequence logos',
			     'logocor'=>'correlation computed on sequence logos',
			     'Nlogocor'=>'Relative width-normalized logocor',
			     'Icor'=>'Pearson correlation computed on Information content',
			     'NIcor'=>'Relative width-normalized Icor',
			     'cov'=>'covariance between residues in aligned columns',
			     'dEucl'=>'Euclidian distance between residue occurrences in aligned columns',
			     'NdEucl'=>'Relative width-normalized dEucl',
			     'NsEucl'=>'similarity derived from Relative width-normalized Euclidian distance',
			     'SSD'=>'Sum of square deviations',
			     'SW'=>'Sandelin-Wasserman',
			     'NSW'=>'Relative width-normalized Sandelin-Wasserman',
			     'match_rank'=>'rank of current match among all sorted matches',
			     'offset'=>'offset between first and second matrices',
			    );



&ScoresAndThresholdsDiv("Matching scores and thresholds",
			"help.compare-matrices.html#return_fields",
			\@matching_scores,
			\%score_descriptions);

################################################################
## Other selectable output fields
my @other_fields = qw(matrix_number
		      matrix_id
		      matrix_name
		      matrix_ac
		      strand
		      direction
		      pos
		      consensus
		      alignments_pairwise
		      alignments_1ton
		     );

################################################################
### send results by email only
print "<p>\n";
#&SelectOutput('email', email_only=>1);
&SelectOutput();
#print "<i>Note: email output is preferred for comparisons with motifs collections</i>\n";

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo single

my $descr1 = "<H2>Comment on the demonstration example 1</H2>\n";
$descr1 .= "<blockquote class ='demo'>";

$descr1 .= "In this demo, we compare a set of motifs discovered with
<i>peak-motifs</i> in a set of 1000 peak regions bound by the mouse
transcription factor Otc4 (Chen et al., 2008).  </p>\n";

#$descr1 .= "Discovered motifs are compared to JASPAR vertebrate
#motifs, and sequences are scanned to predict binding sites.</p>\n";

$descr1 .= "</blockquote>";

#print $query->start_multipart_form(-action=>"compare-matrices_form.cgi");
print $query->start_multipart_form(-action=>"compare-matrices_form.cgi");
$demo_matrices=`cat demo_files/compare-matrices_demo.tf`;
print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'matrix',-default=>$demo_matrices);
#print $query->hidden(-name=>'user_email',-default=>'nobody@nowhere');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<td><b><a href='help.compare-matrices.html'>[MANUAL]</a></B></TD>\n";
##print "<td><b><a href='tutorials/tut_compare-matrices.html'>[TUTORIAL]</a></B></TD>\n";
print "<TD><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################


################################################################
## Comparisons with motif databases
sub DatabaseChoice {
  print '<br/>';

  ## Tasks
  print "<fieldset><legend><b><a href='help.compare-matrices.html#tasks'>Reference matrices (database or custom motif collection)</a></b></legend>";
  print "<p/> ";

  ## load the various databases that can be compared against
  &MatrixDBcheckBox("choice_mode"=>"radiobox");

  print ("<INPUT TYPE='radio' NAME='db_choice' VALUE='custom'>");
  print "Custom motif collection (in TRANSFAC format *)\n";
#  print "<ul>\n";
  print $query->filefield(-name=>'upload_custom_motif_file',
			-size=>10);
#  print "&nbsp;"x6, "Matrices should be in <b>Transfac format</b>";
  print "<br>", "&nbsp;"x6, "* Other formats can be converted with <a href='convert-matrix_form.cgi'><i>convert-matrix</i></a>.";
#  print"</ul>\n";
  print "</p>\n";

  print "</fieldset><p/>";

  print '<p class="clear"></p>';
}


################################################################
## Display a collapsable div with selectable scores and thresholds
sub ScoresAndThresholdsDiv {
  my ($title, $help_file, $field_ref, $field_descr_ref) = @_;
#  my ($title, $help_file, @fields) = @_;

  print "<p class=\"clear\"></p>\n";
#  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>$title</b></div>\n";
#  print "<div id=\"menu101\" class=\"menu_collapsible\">\n";
#  print "<p/><fieldset>\n";

  &FieldsThresholdsTable($help_file, $field_ref, $field_descr_ref);
#  &FieldsThresholdsTable($help_file, @fields);

#  print "</fieldset><p/>";
#  print '</div>';
#  print '</div>';
  print "<p class=\"clear\"></p>\n";
  print "<hr>";
}

################################################################
## Display a table with checkboxes and thresholds for a set of
## specified fields
sub FieldsThresholdsTable {
  my ($help_file, $field_ref, $field_descr_ref) = @_;
  my @fields = @{$field_ref};
  my %field_descr = %{$field_descr_ref};
  print "<table align='center'>\n";
  print $query->th([" <A HREF='".$help_file."'>Output<br>fields</A> ",
		    " <A HREF='".$help_file."'>Lower<BR>Threshold</A> ",
		    " <A HREF='".$help_file."'>Upper<BR>Threshold</A> "]);
  foreach my $field (@fields) {
    my $lth = $default{'lth_'.$field} || "none";
    my $uth = $default{'uth_'.$field} || "none";

    print "<tr valign='middle'>";
#    print "<td>", $field, "</td>\n";
    print "<td>", $query->checkbox(-name=>'return_'.$field,
			   -checked=>$default{'return_'.$field},
			   -label=>''), "&nbsp;", $field, "</td>\n";
    print "<td>", $query->textfield(-name=>'lth_'.$field,
				    -default=>$lth,
				    -size=>5), "</td>\n";
    print "<td>", $query->textfield(-name=>'uth_'.$field,
				    -default=>$uth,
				    -size=>5), "</td>\n";
    print "<td>", $field_descr{$field}, "</td>\n";
    print "</tr>\n";
  }
  print "</table>\n";
}
