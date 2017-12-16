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
use RSAT::matrix;
use RSAT::MatrixReader;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### default values for filling the form
$default{demo_1_descr} = "";
$default{demo_2_descr} = "";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "transfac";
$default{hclust_method}= "average";
$default{merge_stat}= "sum";
$default{metric} = "Ncor";
$default{newick} = "";
$default{random} = "";
$default{quick} = "";
$default{heatmap} = "CHECKED";
$default{radial} = "";
$default{consensus} = "";
$default{label_id} = "";
$default{label_name} = "CHECKED";
$default{label_consensus} = "";
$default{html_title} = "";
$default{collection_label} = "";
$default{'return_w'} = "CHECKED"; $default{'lth_w'} = 5;
$default{'return_cor'} = "CHECKED"; $default{'lth_cor'} = "0.6";
$default{'return_Ncor'} = "CHECKED"; $default{'lth_Ncor'} = "0.4";
$default{'return_logoDP'} = "CHECKED";
$default{'return_NSW'} = "CHECKED";
$default{'return_NsEucl'} = "CHECKED";
my $demo_html_title = "";
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
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>&nbsp: ";
print "<a target='_blank' href='http://pedagogix-tagc.univ-mrs.fr/rsat/data/published_data/Castro_2016_matrix-clustering/Application_4/Ceevee10/demo_jaime.html'>Jaime Castro-Mondragon</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://morgane.bardiaux.fr/'>Morgane Thomas-Chollier</a><sup>t</sup>\n";
print ", <a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>\n";
print "</CENTER>";

## demo 1 description
#print $default{demo_1_descr};

print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

## demo 2 description
#print $default{demo_2_descr};


print $query->start_multipart_form(-action=>"matrix-clustering.cgi");

################################################################
#### Analysis title
print "<hr>";
print "<h2 style='margin-left: 50px;'> Analysis Title ";

print $query->textfield(-id=>'html_title',
-name=>'html_title',
			 -default=>$default{html_title},
			 -size=>30) ."</h2>";



################################################################
#### Matrix specification
print "<hr>";

################################################################
## Query matrices collection 1
print "<h2 style='margin-left: 50px;'> Add up to three collections of PSSMs</h2>";
print "<hr>";

################################################################
#### Set Motif collection label
print "<h2 style='margin-left: 1px;'> Motif Collection\nName";

print $query->textfield(-id=>'collection_label', -name=>'collection_label',
			 -default=>$default{collection_label},
			 -size=>30) ."</h2>";

&GetMatrix('title'=>'Input matrices', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);

print "<hr>";


################################################################
## Query matrices collection 2

################################################################
#### Set Motif collection label
print "<h2 style='margin-left: 1px;'> Motif Collection 2\nName";

print $query->textfield(-id=>'collection_2_label', -name=>'collection_2_label',
			 -default=>$default{collection_2_label},
			 -size=>30) ."</h2>";

&GetSecondMatrix('title'=>'Input matrices 2', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);

print "<hr>";


################################################################
## Query matrices collection 3

################################################################
#### Set Motif collection label
print "<h2 style='margin-left: 1px;'> Motif Collection 3\nName";

print $query->textfield(-id=>'collection_3_label', -name=>'collection_3_label',
			 -default=>$default{collection_3_label},
			 -size=>30) ."</h2>";

&GetThirdMatrix('title'=>'Input matrices 3', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);

print "<hr>";

##############################################################
## Specific options for motif comparison
print "<h2>", "Motif comparison options", ,"</h2>";

## Allow run compare-matrices-quick
print $query->checkbox(-id=>'quick', -name=>'quick',
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
print "<B><A class='iframe' HREF='help.matrix-clustering.html#hclust_method'> Metric to build the trees </A>&nbsp;</B>\n";
print $query->popup_menu(-id=>'metric', -name=>'metric',
 			 -Values=>["cor", "Ncor", "dEucl", "NdEucl", "logocor", "Nlogocor", "logoDP", "Icor", "NIcor", "SSD", "mean_zscore", "rank_mean"],
 			 -default=>$default{metric});
print "<br><br>\n";

## Hierarchical clusterting agglomeration rule
## print "<b>Agglomeration rule</b>";
print "<B><A class='iframe' HREF='help.matrix-clustering.html#hclust_method'> Aglomeration rule </A>&nbsp;</B>\n";
print $query->popup_menu(-id=>'hclust_method', -name=>'hclust_method',
 			 -Values=>["complete", "average", "single", "median", "centroid"],
 			 -default=>$default{hclust_method});
print "<br><br>\n";

## Merge matrix operator
print "<B><A class='iframe' HREF='help.matrix-clustering.html#merge_operator'> Merge matrices </A>&nbsp;</B>\n";
print $query->popup_menu(-id=>'merge_stat', -name=>'merge_stat',
 			 -Values=>["sum", "mean"],
 			 -default=>$default{merge_stat});
print "<br><br>\n";

################################################################
## Specific options for output files
print "<h2>", "Output file options", ,"</h2>";

## Draw heatmap
print $query->checkbox(-name=>'heatmap',
  		       -checked=>$default{heatmap},
  		       -label=>'');
print "&nbsp;<A'><B>Heatmap</B></A>";
print "<br><br>\n";


## Draw the tree with aligned consensuses
# print $query->checkbox(-name=>'alignment_consensuses',
#   		       -checked=>$default{consensus},
#   		       -label=>'');
# print "&nbsp;<A'><B>Export a hierarchical tree with the consensuses aligment.</B></A>";
# print "<br><br>\n";

## Export the trees in Newick format
## By default trees are exported in JSON
print $query->checkbox(-name=>'newick',
  		       -checked=>$default{newick},
  		       -label=>'');
print "&nbsp;<A'><B>Export the trees in Newick format.</B></A>";
print "<br><br>\n";

## Export radial tree
print $query->checkbox(-name=>'radial',
  		       -checked=>$default{radial},
  		       -label=>'');
print "&nbsp;<A'><B>Radial Tree (motif browser)</B></A>";
print "<br><br>\n";

## Negative control: Permute the columns of the input motifs
print $query->checkbox(-name=>'random',
  		       -checked=>$default{random},
  		       -label=>'');
print "&nbsp;<A'><B>Negative control: the input motifs columns are randomly permuted.</B></A>";
print "<br><br>\n";
print "<HR width=550 align=left>\n";



####################################
## Labels displayed in logo trees
#print "<h2>", "Labels displayed in the logo tree", ,"</h2>";

## Label: id
#print $query->checkbox(-name=>'label_id',
#  		       -checked=>$default{label_id},
#  		       -label=>'');
#print "&nbsp;<A'><B>Motif ID</B></A>";
#print "<br><br>\n";

## Label: name
#print $query->checkbox(-name=>'label_name',
#  		       -checked=>$default{label_name},
#  		       -label=>'');
#print "&nbsp;<A'><B>Motif name</B></A>";
#print "<br><br>\n";

## Label: consensus
#print $query->checkbox(-name=>'label_consensus',
#  		       -checked=>$default{label_consensus},
#  		       -label=>'');
#print "&nbsp;<A'><B>Consensus</B></A>";
#print "<br><br>\n";


################################################################
## Send results by email only
print "<p>\n";
&SelectOutput();

################################################################
## Action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
## Demo 1 data

$demo_1_file = "demo_files/RSAT_peak-motifs_Oct4_matrices.tf";
#$demo_1_matrices=`cat ${demo_1_file}`;
$demo_1_matrices = "";
open(my $fh, $demo_1_file);
while(my $row = <$fh>){
    chomp $row;
    $demo_1_matrices .= $row;
    $demo_1_matrices .= "\\n";
}
print "<TD><b>";

print '<script>
function setDemo1(demo_1_matrix){
    $("#reset").trigger("click");
    descr_1 = "<H2>Comment on the demonstration example 1</H2>\n \
    <blockquote class =\'demo_1\'>\
    In this demo, we will analyze with <i>matrix-clustering</i> a \
    set of motifs discovered with <a \
    href=\"peak-motifs_form.cgi\"><i>peak-motifs</i></a> in ChIP-seq binding \
    peaks for the mouse transcription factor Oct4 (data from Chen et al., 2008).  </p>\n</blockquote>";
    
    demo_descr.innerHTML = descr_1;
    html_title.value = "\'Oct4 motifs found in Chen 2008 peak sets\'";
    collection_label.value = "\'Oct4_peak_motifs\'";
    matrix.value = demo_1_matrix;
    demo.value = descr_1;
}
</script>';

print '<button type="button" onclick="setDemo1('. "'$demo_1_matrices'" .')">DEMO (one collection)</button>';
print "</B></TD>\n";



################################################################
## Demo 3 data
$demo_2_file_1 = "demo_files/RSAT_peak-motifs_Oct4_matrices.tf";
$demo_2_file_2 = "demo_files/MEME_ChIP_Oct4_matrices.tf";
$demo_2_file_3 = "demo_files/Homer_l13_mis3_hyper_Oct4_matrices.tf";

$demo_2_matrices_1 = "";
open(my $fh, $demo_2_file_1);
while(my $row = <$fh>){
    chomp $row;
    $demo_2_matrices_1 .= $row;
    $demo_2_matrices_1 .= "\\n";
}
close($fh);
$demo_2_matrices_2 = "";
open($fh, $demo_2_file_2);
while(my $row = <$fh>){
    chomp $row;
    $demo_2_matrices_2 .= $row;
    $demo_2_matrices_2 .= "\\n";
}
close($fh);
$demo_2_matrices_3 = "";
open($fh, $demo_2_file_3);
while(my $row = <$fh>){
    chomp $row;
    $demo_2_matrices_3 .= $row;
    $demo_2_matrices_3 .= "\\n";
}
close($fh);

print "<TD><b>";

print '<script>
function setDemo3(demo_3_matrix, demo_3_matrix_2, demo_3_matrix_3){
    $("#reset").trigger("click");
    descr_3 = "<H2>Comment on the demonstration example 3</H2>\n \
    <blockquote class =\'demo_3\'>\
    In this demo, we will cluster two set of motifs discovered with <a\
    href=\'peak-motifs_form.cgi\'><i>peak-motifs</i></a> and <i>Meme-ChIP</i> in ChIP-seq binding\
    peaks for the transcription factor Oct4 (data from Chen et al.,\
    2008).  </p>\n</blockquote>";
    
    demo_descr.innerHTML = descr_3;
    html_title.value = "\'Oct4 motifs found in Chen 2008 peak sets discovered by peak-motifs and meme-chip\'";
    collection_label.value = "\'Oct4_peak_motifs\'";
    collection_2_label.value = "\'Oct4_Meme-chip\'";
    collection_3_label.value = "\'Oct4_Homer\'";

    matrix.value = demo_3_matrix;
    matrix_2.value = demo_3_matrix_2;
    matrix_3.value = demo_3_matrix_3;
    demo.value = descr_3;
}
</script>';
print '<button type="button" onclick="setDemo3('. "'$demo_2_matrices_1','$demo_2_matrices_2', '$demo_2_matrices_3'" .')">DEMO (three collections)</button>';
print "</B></TD>\n";
#print $query->end_form;


################################################################
## Demo negative control
$demo_2_file = "demo_files/peak-motifs_result_Chen_Oct4_permuted_matrices.tf";

$demo_2_matrices = "";
open($fh, $demo_2_file);
while(my $row = <$fh>){
    chomp $row;
    $demo_2_matrices .= $row;
    $demo_2_matrices .= "\\n";
}
close($fh);

print "<TD><b>";
print '<script>
function setDemo2(demo_3_matrix){
    $("#reset").trigger("click");
    descr_2 = "<H2>Comment on the demonstration example 2</H2>\n \
    <blockquote class =\'demo_2\'>\
    Negative control: we will analyze with <i>matrix-clustering</i> a\
    set of motifs discovered with <a\
    href=\'peak-motifs_form.cgi\'><i>peak-motifs</i></a> in ChIP-seq binding\
    peaks for the mouse transcription factor Oct4 (data from Chen et al.,\
    2008). The columns of the motifs are randomly permuted, conserving thus their\
    information content. Note that poor-complexity motifs (i.e. A-rich) are grouped\
    together.  </p>\n</blockquote>";
    
    demo_descr.innerHTML = descr_2;
    html_title.value = "\'Clustering column-permuted matrices discovered in Oct4 ChIP-seq\'";
    collection_label.value = "\'Oct4_peak_motifs_permuted\'";
    matrix.value = demo_3_matrix;
    demo.value = descr_2;
}
</script>';
print '<button type="button" onclick="setDemo2('. "'$demo_2_matrices'" .')">DEMO (negative control)</button>';
print "</B></TD>\n";


print "<td><b><a class='iframe' href='help.matrix-clustering.html'>[MANUAL]</a></B></TD>\n";
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
                         NcorS
                         logoDP
			 logocor
			 Nlogocor
			 Icor
			 NIcor
			 dEucl
			 NdEucl
			 NsEucl
			 SSD
			);

  my %score_descriptions = ('w'=>'Width = number of aligned columns',
			    'cor'=>'Pearson correlation (computed on residue occurrences in aligned columns)',
			    'Ncor'=>'Relative width-normalized Pearson correlation',
			    'NcorS'=>'Relative width-normalized Pearson correlation of the smallest alignment',
			    'logoDP'=>'dot product of sequence logos',
			    'logocor'=>'correlation computed on sequence logos',
			    'Nlogocor'=>'Relative width-normalized logocor',
			    'Icor'=>'Pearson correlation computed on Information content',
			    'NIcor'=>'Relative width-normalized Icor',
			    'dEucl'=>'Euclidian distance between residue occurrences in aligned columns',
			    'NdEucl'=>'Relative width-normalized dEucl',
			    'NsEucl'=>'similarity derived from Relative width-normalized Euclidian distance',
			    'SSD'=>'Sum of square deviations',
      );

  &ThresholdsDiv(" Thresholds to define the clusters",
			  "help.matrix-clustering.html#return_fields",
			  \@matching_scores,
			  \%score_descriptions);
}


################################################################
## Display a collapsable div with selectable scores and thresholds
sub ThresholdsDiv {
  my ($title, $help_file, $field_ref, $field_descr_ref) = @_;

  print "<p class=\"clear\"></p>\n";
  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>",$title,"</b></div>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible\">\n";
  print "<p><fieldset>\n";

  &FieldsThresholdsTableMC($help_file, $field_ref, $field_descr_ref);

  print "</fieldset><p/>";
  print '</div><p class="clear"></p>';
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
  print $query->th([" <A HREF='".$help_file."'>Metrics<br></A> ",
		    " <A HREF='".$help_file."'>Lower<br>Threshold</A> ",
		    " <A HREF='".$help_file."'>Upper<br>Threshold</A> "]);
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
