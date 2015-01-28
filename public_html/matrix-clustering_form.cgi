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
$default{matrix_format} = "transfac";
$default{hclust_method}="average";
$default{metric} = "Ncor";
$default{newick} = "";
$checked{$default{bg_method}} = "CHECKED";
$default{heatmap} = "CHECKED";
$default{labels} = "name";
$default{'return_w'} = "CHECKED"; $default{'lth_w'} = 5;
$default{'return_cor'} = "CHECKED"; $default{'lth_cor'} = "0.6";
$default{'return_Ncor'} = "CHECKED"; $default{'lth_Ncor'} = "0.4";
$default{'return_logoDP'} = "CHECKED";
$default{'return_NSW'} = "CHECKED";
$default{'return_NsEucl'} = "CHECKED";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
  if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
}



&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
### print the form ###

################################################################
### header
&RSA_header("matrix-clustering", "form");
print "<CENTER>";
print "Identify groups (clusters) of similarities between a set of motifs and align them.<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</a><sup>t</sup>\n";
print ", <a target='_blank' href='http://www.ibens.ens.fr/spip.php?article26&lang=en'>Denis Thieffry</a><sup>t</sup>\n";
print "and <a target='_blank' href='http://biologie.univ-mrs.fr/view-data.php?id=202'>Carl Herrmann</a><sup>ct</sup>\n";
print "and <a target='_blank'>Jaime Castro</a><sup>cit</sup>\n";
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

##############################################################
## Specific options for motif comparison
print "<h2>", "Motif comparison options", ,"</h2>";

## Allow run compare-matrices-quick
print $query->checkbox(-name=>'quick',
  		       -checked=>$default{quick},
  		       -label=>'');
print "&nbsp;<A'><B>Motif comparison with <i>compare-matrices-quick</i> (100 times faster). Only for <strong>Ncor</strong> and <strong>Cor</strong>.</B></A>";
print "<br><br>\n";
#print "<HR width=550 align=left>\n";

## Selection of output fields and thresholds
&PrintMatrixClusteringMatchingScores();


###############################################################
## Specific options for matrix-clustering
print "<h2>", "Clustering options", ,"</h2>";

## Metric selected to build the hierarchical tree
## print "<b>Metric to build the tree.</b>";
print "<B><A HREF='help.matrix-clustering.html#hclust_method'> Metric to build the trees </A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'hclust_method',
 			 -Values=>["cor", "Ncor"],
 			 -default=>$default{metric});
print "<br><br>\n";

## Hierarchical clusterting agglomeration rule
## print "<b>Agglomeration rule</b>";
print "<B><A HREF='help.matrix-clustering.html#hclust_method'> Aglomeration rule </A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'hclust_method',
 			 -Values=>["complete", "average", "single"],
 			 -default=>$default{hclust_method});
print "<br><br>\n";

################################################################
## Specific options for output files
print "<h2>", "Output file options", ,"</h2>";

## Draw heatmap
print $query->checkbox(-name=>'heatmap',
  		       -checked=>$default{heatmap},
  		       -label=>'');
print "&nbsp;<A'><B>Draw a heatmap showing the distances between the motifs.</B></A>";
print "<br><br>\n";
#print "<HR width=550 align=left>\n";

## Export the trees in Newick format
## By default trees are exported in JSON
print $query->checkbox(-name=>'Newick',
  		       -checked=>$default{newick},
  		       -label=>'');
print "&nbsp;<A'><B>Export the trees in Newick format.</B></A>";
print "<br><br>\n";
#print "<HR width=550 align=left>\n";


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
$demo_file = "demo_files/peak-motifs_Oct4_matrices.tf";
$demo_matrices=`cat ${demo_file}`;
print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'matrix',-default=>$demo_matrices);
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<td><b><a href='help.compare-matrices.html'>[MANUAL]</a></B></TD>\n";
print "<td><b><a href='help.matrix-clustering.html'>[MANUAL]</a></B></TD>\n";
print "<TD><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></B></TD>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

################################################################
## Print the many scores supported by compare-matrices, which are also
## available for matrix-clustering.
sub PrintMatrixClusteringMatchingScores {
  my @matching_scores = qw(w
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

  my %score_descriptions = ('w'=>'Width = number of aligned columns',
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

  &ThresholdsDiv("Matrix comparison scores and thresholds",
			  "help.compare-matrices.html#return_fields",
			  \@matching_scores,
			  \%score_descriptions);
}


################################################################
## Display a collapsable div with selectable scores and thresholds
sub ThresholdsDiv {
  my ($title, $help_file, $field_ref, $field_descr_ref) = @_;

  print "<p class=\"clear\"></p>\n";
  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>",$title,"</b>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible\">\n";
  print "<p/><fieldset>\n";

  &FieldsThresholdsTableMC($help_file, $field_ref, $field_descr_ref);

  print "</fieldset><p/>";
  print '</div></div><p class="clear"></p>';
  print "<hr>";
}

################################################################
## Display a table with checkboxes and thresholds for a set of
## specified fields
sub FieldsThresholdsTableMC {
  my ($help_file, $field_ref, $field_descr_ref) = @_;
  my @fields = @{$field_ref};
  my %field_descr = %{$field_descr_ref};
  print "<table align='center'>\n";
  print $query->th([" <A HREF='".$help_file."'>Metrics</A> ",
		    " <A HREF='".$help_file."'>Lower<BR>Threshold</A> ",
		    " <A HREF='".$help_file."'>Upper<BR>Threshold</A> "]);
  foreach my $field (@fields) {
    my $lth = $default{'lth_'.$field} || "none";
    my $uth = $default{'uth_'.$field} || "none";

    print "<tr valign='middle'>";
    print "<td>".$field."</td>\n";
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
