#!/usr/bin/env perl
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
$default{matrix1}="";
$default{matrix_file}="";
$default{matrix_format1} = "transfac";
$default{matrix_format2} = "transfac";
$default{matrix_format3} = "transfac";
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
&RSA_header_bootstrap("matrix-clustering", "form");

print $query->start_multipart_form(-action=>"matrix-clustering.cgi");

print '
 <!-- Form with bootstrap -->
<div class="container">
	<div class="row">
        <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
            <div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
              <div class="list-group">
                <a href="#" class="list-group-item active text-center">
                  <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>Matrix Clustering
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Mandatory inputs
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tags fa-2x"></i></h4><br/>Optional inputs
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tasks fa-2x"></i></h4><br/>Advanced options
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-play-circle fa-2x"></i></h4><br/>Run analysis
                </a>
              </div>
            </div>
            <div class="col-lg-9 col-md-9 col-sm-9 col-xs-9 bhoechie-tab">


 <!-- ################################################################ -->
 <!-- ### info ### -->

                <div class="bhoechie-tab-content active">

                     <h2> <img src="images/RSAT_logo.jpg" style="max-width:150px;max-height:60px;padding-bottom:10px" alt="RSAT server" border="0"></img>
                      matrix-clustering</h2>
                    <span class="fa-stack fa-lg">
  							<i class="fa fa-info-circle fa-stack-1x"></i>
					</span>
					Identify groups (clusters) of similarities between one (or several) collections of motifs, and align them. <br>
                    <span class="fa-stack fa-lg">
 							 <i class="fa fa-user fa-stack-1x"></i>
					</span>

					<a target="_blank" href="http://folk.uio.no/jamondra/">Jaime A Castro-Mondragon</a>, <a target="_blank" href="http://morgane.bardiaux.fr/">Morgane Thomas-Chollier</a>, <a target="_blank" href="http://jacques.van-helden.perso.luminy.univ-amu.fr/ ">Jacques van Helden</a>.<br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-folder-open fa-stack-1x"></i>
					</span>
					Sample output<br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-book fa-stack-1x"></i>
					</span>
					<a class="iframe" href="help.matrix-clustering.html">User Manual</a><br>

					<!--span class="fa-stack fa-lg">
  							<i class="fa fa-graduation-cap fa-stack-1x"></i>
					</span>
					<a class="iframe" href="help.matrix-clustering.html">Tutorial</a><br-->

					<span class="fa-stack fa-lg">
  							<i class="fa fa-twitter fa-stack-1x"></i>
					</span>
					<a href="https://twitter.com/rsatools" target="_blank">Ask a question to the RSAT team</a><br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-pencil fa-stack-1x"></i>
					</span>
					Cite the publication: <a href="https://twitter.com/rsatools" target="_blank"></a><br>
					<div class="panel panel-default">
  						<div class="panel-body">
    						Castro-Mondragon JA, Jaeger S, Thieffry D, Thomas-Chollier M#, van Helden J#. <i>"RSAT matrix-clustering: dynamic exploration and redundancy reduction of transcription factor binding motif collections."</i>, Nucleic Acid Research, 45:13 e119 (2017) <a href="https://www.ncbi.nlm.nih.gov/pubmed/28591841" target="_blank">[Pubmed]</a><a href="https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkx314" target="_blank">[Full text]</a>

  						</div>
					</div>

                </div>

 <!-- ################################################################ -->
 <!-- ### mandatory inputs ### -->
                <div class="bhoechie-tab-content">
 <!-- title -->
<div class="panel panel-danger">
 <div class="panel-heading">Analysis Title <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Title that will be displayed at the top of the report page." data-original-title=""></i></div>

<div class="panel-body">
<div class="form-group">
 ';
print $query->textfield(-id=>'html_title',
-name=>'html_title', -class=>'form-control',-placeholder=>'Provide a Title for this analysis ', -required=>'true',
			 -default=>$default{html_title}) .'

		</div>
			 </div>
			 </div>

<!-- Motifs -->
<div class="panel panel-danger">
 <div class="panel-heading">Motif Collection
  <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Input here the motif collection of interest to be clustered, you can either paste the motif in the text box or upload it from your computer" data-original-title=""></i>
</div>
			<div class="panel-body">
 <div class="form-group">';
print $query->textfield(-id=>'collection_1_label', -name=>'collection_1_label', -class=>'form-control', -placeholder=>'Provide a name for this Motif Collection',
			 -default=>$default{collection_1_label});
print '</div>';
#&GetMatrix_bootstrap('title'=>'Matrix Format', 'nowhere'=>1,'no_pseudo'=>1, consensus=>1);
&MultiGetMatrix_bootstrap('title'=>'Matrix Format','mat_num'=>1);

print '</div></div></div>

 <!-- ################################################################ -->
 <!-- ### optional inputs ### -->
                <div class="bhoechie-tab-content">

<div id="accordion" role="tablist">
  <div class="card">
    <div class="card-header" role="tab" id="headingOne">
      <h5> <i class="fa fa-tags"></i>
        <a data-toggle="collapse" href="#collapseOne" aria-expanded="true" aria-controls="collapseOne">
          Add a second Motif collection
        </a>
      </h5>
    </div>

    <div id="collapseOne" class="collapse" role="tabpanel" aria-labelledby="headingOne" data-parent="#accordion">
      <div class="card-body">
 <div class="panel panel-warning">
 <div class="panel-heading">Motif Collection 2
  <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Input here the motif collection of interest to be clustered, you can either paste the motif in the text box or upload it from your computer" data-original-title=""></i>
</div>
			<div class="panel-body">
 <div class="form-group">';
print $query->textfield(-id=>'collection_2_label', -name=>'collection_2_label', -class=>'form-control', -placeholder=>'Provide a name for this Motif Collection',
			 -default=>$default{collection_2_label});
print '</div>';
&MultiGetMatrix_bootstrap('title'=>'Matrix Format','mat_num'=>2);

print '</div></div>
      </div>
    </div>
  </div>
  <div class="card">
    <div class="card-header" role="tab" id="headingTwo">
      <h5 class="mb-0"> <i class="fa fa-tags"></i>
        <a class="collapsed" data-toggle="collapse" href="#collapseTwo" aria-expanded="false" aria-controls="collapseTwo">
          Add a third Motif collection
        </a>
      </h5>
    </div>
    <div id="collapseTwo" class="collapse" role="tabpanel" aria-labelledby="headingTwo" data-parent="#accordion">
      <div class="card-body">
 <div class="panel panel-warning">
 <div class="panel-heading">Motif Collection 3
  <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Input here the motif collection of interest to be clustered, you can either paste the motif in the text box or upload it from your computer" data-original-title=""></i>
</div>
			<div class="panel-body">
 <div class="form-group">';
print $query->textfield(-id=>'collection_3_label', -name=>'collection_3_label', -class=>'form-control', -placeholder=>'Provide a name for this Motif Collection',
			 -default=>$default{collection_3_label});
print '</div>';
&MultiGetMatrix_bootstrap('title'=>'Matrix Format','mat_num'=>3);

print '
      </div>
    </div>
  </div>
</div>
                </div></div>
                </div>

 <!-- ################################################################-->
 <!-- ### advanced options ###-->

                   <!-- ADVANCED OPTIONS -->

                <div class="bhoechie-tab-content">

  <div id="accordion" role="tablist">


 <!-- Matrix clustering-->
  <div class="card">
    <div class="card-header" role="tab" id="headingFour">
      <h5 class="mb-0">  <i class="fa fa-tasks"></i>
        <a class="collapsed" data-toggle="collapse" href="#collapseFour" aria-expanded="false" aria-controls="collapseFour">
          Options for the clustering step
        </a>
      </h5>
    </div>
    <div id="collapseFour" class="collapse" role="tabpanel" aria-labelledby="headingFour" data-parent="#accordion">
      <div class="card-body">

<div class="panel panel-warning">
 <div class="panel-heading">Clustering options</div>
			<div class="panel-body">';

#Metric selected to build the hierarchical tree
 print "<div class='form-row'>
 <label for='metric' class='col-sm-9 control-label'>Metric for the motif-to-motif similarity matrix <A class='badge badge-primary iframe' HREF='help.matrix-clustering.html#metric_build_tree-metric'>Info</a></label>\n";
  print "<div class='col-sm-3'>";

print $query->popup_menu(-id=>'metric', -name=>'metric',
 			 -Values=>["cor", "Ncor", "dEucl", "NdEucl", "logocor", "Nlogocor", "logoDP", "Icor", "NIcor", "SSD", "mean_zscore", "rank_mean"],
 			 -class=>'form-control',
 			 -default=>$default{metric});
print "</div></div>\n";

# Hierarchical clusterting agglomeration rule
print "<div class='form-row'>
 <label for='hclust_method' class='col-sm-9 control-label'>Agglomeration (linkage) rule to build the hierachical tree<A class='badge badge-primary iframe' HREF='help.matrix-clustering.html#hclust_method'>Info</a></label>\n";
  print "<div class='col-sm-3'>";

 print $query->popup_menu(-id=>'hclust_method', -name=>'hclust_method',
  			 -Values=>["complete", "average", "single", "median", "centroid"],
  			  -class=>'form-control',
 			 -default=>$default{hclust_method});
print "</div></div>\n";

# Merge matrix operator
print "<div class='form-row'>
 <label for='merge_stat' class='col-sm-9 control-label'>Merge matrices <A class='badge badge-primary iframe' HREF='help.matrix-clustering.html#merge_operator'>Info</a></label>\n";
  print "<div class='col-sm-3'>";

print $query->popup_menu(-id=>'merge_stat', -name=>'merge_stat',
  			 -Values=>["sum", "mean"],
  			  -class=>'form-control',
 			 -default=>$default{merge_stat});
print '
      </div>
    </div>
  </div>
</div>';

print'
<!-- Compare matrices-->
<div class="card">
  <div class="card-header" role="tab" id="headingThree">
    <h5> <i class="fa fa-tasks"></i>
      <a data-toggle="collapse" href="#collapseThree" aria-expanded="true" aria-controls="collapseThree">
       Options for the motif comparison step (program: compare-matrices)
      </a>
    </h5>
  </div>

  <div id="collapseThree" class="collapse" role="tabpanel" aria-labelledby="headingThree" data-parent="#accordion">
    <div class="card-body">
    <div class="panel panel-warning">
    <div class="panel-heading">Motif comparison options</div>
    <div class="panel-body">';


# Allow run compare-matrices-quick
#print $query->checkbox(-id=>'quick', -name=>'quick',
#   		       -checked=>$default{quick},
#   		       -label=>'');
# print "Motif comparison with <i>compare-matrices-quick</i> (100 times faster). Only for <strong>Ncor</strong> and <strong>Cor</strong>.";
#print "<hr>";
# Selection of output fields and thresholds
&PrintMatrixClusteringMatchingScores();

print '</div></div>
    </div>
  </div>
  </div>

                </div></div>


<!--close panel-->
</div></div>
 </div>

 <!--################################################################-->
 <!--### output & run ###-->

                <div class="bhoechie-tab-content">

 <!-- ## Specific options for output files-->
<div class="panel panel-info">
 <div class="panel-heading">Output options</div>

<div class="panel-body">
<div class="form-group">

<!--A class="badge badge-primary iframe" HREF="help.matrix-clustering.html#merge_operator">Info</a></label></h5-->


<!--## Draw heatmap-->
 <div class="form-check">
    <label class="form-check-label">';
    print $query->checkbox(-name=>'heatmap',
  		       -checked=>$default{heatmap},
  		       -class=>form-check-input,
  		        -label=>'');
print   'Heatmap</label></div>';


## Export the trees in Newick format
## By default trees are exported in JSON
#print '<div class="form-check">
#    <label class="form-check-label">';
#    print $query->checkbox(-name=>'newick',
#  		       -checked=>$default{newick},
#  		       -class=>form-check-input,
#  		        -label=>'');
#print   'Export the trees in Newick format</label></div>';

## Export radial tree
print '<div class="form-check">
    <label class="form-check-label">';
    print $query->checkbox(-name=>'radial',
  		       -checked=>$default{radial},
  		       -class=>form-check-input,
  		        -label=>'');
print   'Export Radial Tree (motif browser)</label></div>';

## Negative control: Permute the columns of the input motifs
print '<div class="form-check">
    <label class="form-check-label">';
    print $query->checkbox(-name=>'random',
  		       -checked=>$default{random},
  		       -class=>form-check-input,
  		        -label=>'');
print   'Negative control: the input motifs columns are randomly permuted.</label></div>';

print "<HR>\n";

################################################################
## Send results by email only
print "<p>\n";
&SelectOutput();
print " </div>
  </div> </div>";


################################################################
## Action buttons
#print "<TABLE class='formbutton'>\n";
#print "<TR VALIGN=MIDDLE>\n";
#print "<TD>", $query->submit(-label=>"GO", -class=>"btn btn-success"), "</TD>\n";
print $query->submit(-label=>"GO", -class=>"btn btn-success", -type=>"button");
print " ";
print $query->reset(-id=>"reset",-class=>"btn btn-warning", -type=>"button");
print $query->end_form;

 print ' </div>
        </div>
  </div>
</div>
';

################################################################
## Demo area
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr' class='col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container'></div>";



#
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

print '<script>
function setDemo1(demo_1_matrix){
    $("#reset").trigger("click");
    descr_1 = "<H4>Demonstration: input = one collection</H4>\n \
    <blockquote class =\'blockquote text-justify small\'>\
    In this demo, we will cluster with <i>matrix-clustering</i> one \
    set of motifs discovered with the program <a \
    href=\"peak-motifs_form.cgi\"><i>peak-motifs</i></a> in ChIP-seq \
    peaks for the mouse transcription factor Oct4 (data from Chen et al., 2008). Check the panel <b>Mandatory inputs</b> and then <b>Run analysis</b></blockquote>";

    demo_descr.innerHTML = descr_1;
    html_title.value = "\'Oct4 motifs found in Chen 2008 peak sets\'";
    collection_1_label.value = "\'Oct4_peak_motifs\'";
    matrix1.value = demo_1_matrix;
    demo.value = descr_1;
}
</script>';

print ' <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container">

<button type="button" class="btn btn-info" onclick="setDemo1('. "'$demo_1_matrices'" .')">DEMO (one collection)</button> ';



################################################################
## Demo 2 data
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
function setDemo2(demo_2_matrix, demo_2_matrix_2, demo_2_matrix_3){
    $("#reset").trigger("click");
    descr_2 = "<H4>Demonstration: input = three collections</H4>\n \
    <blockquote class =\'demo_2 blockquote text-justify small\'>\
    In this demo, we will cluster three sets of motifs discovered with <a\
    href=\'peak-motifs_form.cgi\'><i>peak-motifs</i></a>, <i>Meme-ChIP</i> and <i>HOMER</i> in ChIP-seq\
    peaks for the mouse transcription factor Oct4 (data from Chen et al.,\
    2008). Check the panel <b>Mandatory inputs</b>, then <b>Optional inputs</b> and finally <b>Run analysis</b></p>\n</blockquote>";

    demo_descr.innerHTML = descr_2;
    html_title.value = "\'Oct4 motifs found in Chen 2008 peak sets discovered by peak-motifs, meme-chip and Homer\'";
    collection_1_label.value = "\'Oct4_peak_motifs\'";
    collection_2_label.value = "\'Oct4_Meme-chip\'";
    collection_3_label.value = "\'Oct4_Homer\'";

    matrix1.value = demo_2_matrix;
    matrix2.value = demo_2_matrix_2;
    matrix3.value = demo_2_matrix_3;
    demo.value = descr_2;
}
</script>';
print '<button type="button"class="btn btn-info" onclick="setDemo2('. "'$demo_2_matrices_1','$demo_2_matrices_2', '$demo_2_matrices_3'" .')">DEMO (three collections)</button> ';

#print $query->end_form;


################################################################
## Demo negative control
$demo_3_file = "demo_files/peak-motifs_result_Chen_Oct4_permuted_matrices.tf";

$demo_3_matrices = "";
open($fh, $demo_3_file);
while(my $row = <$fh>){
    chomp $row;
    $demo_3_matrices .= $row;
    $demo_3_matrices .= "\\n";
}
close($fh);

print "<TD><b>";
print '<script>
function setDemo3(demo_3_matrix){
    $("#reset").trigger("click");
    descr_3 = "<H4>Demonstration: input motifs permuted</H4>\n \
    <blockquote class =\'demo_3 blockquote text-justify small\'>\
    Negative control: we will analyze with <i>matrix-clustering</i> a\
    set of motifs discovered with <a\
    href=\'peak-motifs_form.cgi\'><i>peak-motifs</i></a> in ChIP-seq\
    peaks for the mouse transcription factor Oct4 (data from Chen et al.,\
    2008). The columns of the motifs are randomly permuted, thus conserving  their\
    information content. Note that poor-complexity motifs (i.e. A-rich) will be grouped\
    together. Check the panel <b>Mandatory inputs</b> and then <b>Run analysis</b></p>\n</blockquote>";

    demo_descr.innerHTML = descr_3;
    html_title.value = "\'Clustering column-permuted matrices discovered in Oct4 ChIP-seq\'";
    collection_1_label.value = "\'Oct4_peak_motifs_permuted\'";
    matrix1.value = demo_3_matrix;
    demo.value = descr_3;
}
</script>';
print '<button type="button"  class="btn btn-info" onclick="setDemo3('. "'$demo_3_matrices'" .')">DEMO (negative control)</button>';
print "</B></TD>\n";



print "</div> ";

print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

################################################################
## Print the many scores supported by compare-matrices, which are also
## available for matrix-clustering.
sub PrintMatrixClusteringMatchingScores {
  my @matching_scores_main = qw(w
			 cor
			 Ncor
			);
my @matching_scores = qw(logoDP
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

  &ThresholdsDiv("Main Thresholds to define the clusters",
    			  "help.matrix-clustering.html#return_fields",
    			  \@matching_scores_main,
    			  \%score_descriptions);

  &ThresholdsDiv("Supplementary thresholds to define the clusters (slower version of the program)",
			  "help.matrix-clustering.html#return_fields",
			  \@matching_scores,
			  \%score_descriptions);
}


################################################################
## Display a collapsable div with selectable scores and thresholds
sub ThresholdsDiv {
  my ($title, $help_file, $field_ref, $field_descr_ref) = @_;
  print "<p><fieldset class='form-group'>\n<b>".$title."</b>";
  &FieldsThresholdsTableMC($help_file, $field_ref, $field_descr_ref);
print '</fieldset><p/>';
}

################################################################
## Display a table with checkboxes and thresholds for a set of
## specified fields
sub FieldsThresholdsTableMC {
  my ($help_file, $field_ref, $field_descr_ref) = @_;
  my @fields = @{$field_ref};
  my %field_descr = %{$field_descr_ref};
  print "<table class='table table-striped table-sm' style=';font-size:12px'>\n";
  print $query->th(["Metrics",
		    "Lower<br>Threshold",
		    "Upper<br>Threshold",
		    "description"]);

  foreach my $field (@fields) {
    my $lth = $default{'lth_'.$field} || "none";
    my $uth = $default{'uth_'.$field} || "none";

    print "<tr valign='middle'>";
    print "<td>".$field."</td>\n";
    print "<td>", $query->textfield(-name=>'lth_'.$field, -class=>'form-control',
				    -default=>$lth,
				    -size=>5), "</td>\n";
    print "<td>", $query->textfield(-name=>'uth_'.$field,-class=>'form-control',
				    -default=>$uth,
				    -size=>5), "</td>\n";
    print "<td>", $field_descr{$field}, "</td>\n";
    print "</tr>\n";
  }
  print "</table>\n";
}
