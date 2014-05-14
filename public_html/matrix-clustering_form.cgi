#!/usr/bin/perl
#### this cgi script fills the HTML form for the program matrix-clustering
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
$default{'return_cor'} = "CHECKED"; $default{'lth_cor'} = "none";
$default{'return_Ncor'} = "CHECKED"; $default{'lth_Ncor'} = "none";
$default{'return_match_rank'} = "CHECKED"; $default{'uth_match_rank'} = 100;
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
&RSA_header("matrix-clustering", "form");
print "<CENTER>";
print "Identify clusters of similar matrices and build consensus motifs by merging the matrices that belong to the same cluster.<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</a><sup>t</sup>\n";
print ", <a target='_blank' href='http://www.ibens.ens.fr/spip.php?article26&lang=en'>Denis Thieffry</a><sup>t</sup>\n";
print "and <a target='_blank' href='http://biologie.univ-mrs.fr/view-data.php?id=202'>Carl Herrmann</a><sup>ct</sup>\n";
print "</CENTER>";

## demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"matrix-clustering.cgi");


################################################################
#### Matrix specification
print "<hr>";

################################################################
## Query matrices
&GetMatrix('title'=>'Query matrices', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);
print "<hr>";

################################################################
## Specific options for matrix-clustering

################################################################
## Selection of output fields and thresholds
&PrintMatrixMatchingScores();

################################################################
## Send results by email only
print "<p>\n";
&SelectOutput();

################################################################
## Action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
## Demo data
my $descr1 = "<H2>Comment on the demonstration example 1</H2>\n";
$descr1 .= "<blockquote class ='demo'>";

$descr1 .= "In this demo, we will apply <i>matrix-clustering</i> to a
set of motifs discovered with <i>peak-motifs</i> in ChIP-seq binding
peaks for the mouse transcription factor Otc4 (data from Chen et al.,
2008).  </p>\n";

$descr1 .= "</blockquote>";

print $query->start_multipart_form(-action=>"matrix-clustering_form.cgi");
$demo_file = "demo_files/peak-motifs_result_Chen_Oct4_matrices.tf";
$demo_matrices=`cat ${demo_file}`;
print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'matrix',-default=>$demo_matrices);
#print $query->hidden(-name=>'user_email',-default=>'nobody@nowhere');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<td><b><a href='help.compare-matrices.html'>[MANUAL]</a></B></TD>\n";
print "<TD><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

