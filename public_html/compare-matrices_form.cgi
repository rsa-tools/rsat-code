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
use RSAT::matrix;
use RSAT::MatrixReader;
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
#$default{compare_motif_database}="footprintDB-plants"; #"jaspar_core_nonredundant_vertebrates";
$default{compare_motif_database}="jaspar_core_nonredundant_vertebrates"; ## I (JvH) SHOULD ADAPT THIS IN ORDER TO PRESENT DIFFERENT DATABASES DEPENDING ON TAXON SPECIFICITY OF THE SERVER (2015-03-13)

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
print "<a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</a><sup>t</sup>\n";
print ", <a target='_blank' href='http://www.ibens.ens.fr/spip.php?article26&lang=en'>Denis Thieffry</a><sup>t</sup>\n";
print "and <a target='_blank' href='http://biologie.univ-mrs.fr/view-data.php?id=202'>Carl Herrmann</a><sup>ct</sup>\n";
print "</CENTER>";

## demo description
#print $default{demo_descr1};
print "<div id='demo_descr1'>$default{demo_descr}</div>";

print $query->start_multipart_form(-action=>"compare-matrices.cgi");


################################################################
#### Matrix specification
print "<hr>";

################################################################
## Query matrices
&GetMatrix('title'=>'Query matrices', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);
print "<hr>";

################################################################
## Supported databases

&DatabaseChoice_select2();

#&DatabaseChoice();
print "<hr>";

################################################################
## Background model
#my %bg_params =("from_matrix" => 1,
#		"bg_input"=>0,
#		"simple"=>1,
#		"no_bg_pseudo"=>1,
#	       );
#&GetBackgroundModel(%bg_params);
#print "<hr>";

################################################################
## Selection of output fields and thresholds
&PrintMatrixMatchingScores();
print "<hr>";

################################################################
### send results by email only
#print "<p>\n";
#&SelectOutput('email', email_only=>1);
&SelectOutput();
#print "<i>Note: email output is preferred for comparisons with motifs collections</i>\n";

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
### data for the demo single

$demo_matrices = "";
open(my $fh, "demo_files/compare-matrices_demo.tf");
while(my $row = <$fh>){
    chomp $row;
    $demo_matrices .= $row;
    $demo_matrices .= "\\n";
}

print "<TD><b>";
print '<script>
function setDemo(demo_matrices){
    $("#reset").trigger("click");
    descr1 = "<H2>Comment on the demonstration example 1</H2>\n<blockquote class =\"demo\">In this demo, we compare a set of motifs discovered with <i>peak-motifs</i> in a set of 1000 peak regions bound by the mouse transcription factor Otc4 (Chen et al., 2008).  </p>\n</blockquote>";
    demo_descr1.innerHTML = descr1;
    matrix.value = demo_matrices;
    $("#db_choice").val("' . $default{compare_motif_database} . '").trigger("change");
}
</script>';

print '<button type="button" onclick="setDemo('. "'$demo_matrices'" .')">DEMO</button>';
print "</B></TD>\n";


print "<td><b><a class='iframe' href='help.compare-matrices.html'>[MANUAL]</a></B></TD>\n";
##print "<td><b><a href='tutorials/tut_compare-matrices.html'>[TUTORIAL]</a></B></TD>\n";
##print "<TD><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print "<div style='display:none'><div id='matrix_descr'>";
&DisplayMatrixDBchoice_select2("mode" => "list");
print "</div></div>";
print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################


################################################################
## Comparisons with motif databases

sub DatabaseChoice_select2 {
  print '<br/>';

  ## Tasks
  print "<fieldset><legend><b><a class='iframe' href='help.compare-matrices.html#tasks'>Reference matrices (database or custom motif collection)</a></b></legend>";
  print "<p> ";


  &MotifSelection();
  print "<br/><b>Or</b> Custom motif(s) in TRANSFAC format *\n";
  print $query->filefield(-name=>'upload_custom_motif_file',
			  -size=>10, -onchange=>'fileupload(event)');
  print "<br>", "&nbsp;"x6, "* Other formats can be converted with <a class='iframe' href='convert-matrix_form.cgi'><i>convert-matrix</i></a>.";
  print "</p>\n";

  print "</fieldset>";

  print '<p class="clear"></p>';
}

