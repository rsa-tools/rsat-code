#!/usr/bin/perl
#### updated by Bruno Jan2018

BEGIN{
    if ($0 =~ /([^(\/)]+)$/) {
        push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{query_classes} = "";
$default{upload_query_classes} = "";
$default{ref_classes} = "";
$default{upload_ref_classes} = "";

$default{occ} = "checked";
$default{lth_occ} = 1;
$default{uth_occ} = "none";
$default{freq} = "checked";
$default{sig} = "checked";
$default{lth_sig} = 0;
$default{uth_sig} = "none";
$default{proba} = "checked";
$default{freq} = "checked";
$default{jac} = "checked";
$default{entropy} = "checked";
$default{members} = "";
$default{sort_key} = "sig";
$default{pop_size} = "auto";
$default{entropy} = "checked";
$default{jac} = "checked";

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

### print the form ###
&RSA_header_bootstrap("compare-classes", 'form');

print $query->start_multipart_form(-action=>"compare-classes.cgi");

print '
 <!-- Form with bootstrap -->
<div class="container">
        <div class="row">
        <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
            <div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
              <div class="list-group">
                <a href="#" class="list-group-item active text-center">
                  <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>compare-classes
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

    <h2> <img src="images/RSAT_logo.jpg" style="max-width:150px;max-height:60px;padding-bottom:10px" alt="RSAT server" border="0"></img>compare-classes</h2>
    <span class="fa-stack fa-lg">
        <i class="fa fa-info-circle fa-stack-1x"></i>
    </span>
    Compare two classifications (clustering results, functional classes, etc), and assess the statistical significance of common members between pairs of classes.<br>
    <span class="fa-stack fa-lg">
        <i class="fa fa-user fa-stack-1x"></i>
    </span>
    <a target="_blank" href="http://jacques.van-helden.perso.luminy.univ-amu.fr/ ">Jacques van Helden</a> with help from Joseph Tran.<br>
    <span class="fa-stack fa-lg">
        <i class="fa fa-folder-open fa-stack-1x"></i>
    </span>
    Sample output<br>
    <span class="fa-stack fa-lg">
        <i class="fa fa-book fa-stack-1x"></i>
    </span>
    <a class="iframe" href="help.compare-classes.html">User Manual</a><br>
    <!--span class="fa-stack fa-lg">
        <i class="fa fa-graduation-cap fa-stack-1x"></i>
    </span>
    <a class="iframe" href="help.compare-classes.html">Tutorial</a><br-->
    <span class="fa-stack fa-lg">
        <i class="fa fa-twitter fa-stack-1x"></i>
    </span>
    <a href="https://twitter.com/rsatools" target="_blank">Ask a question to the RSAT team</a><br>
    <span class="fa-stack fa-lg">
        <i class="fa fa-pencil fa-stack-1x"></i>
    </span>
    <!--Cite the publication: <a href="https://twitter.com/rsatools" target="_blank"></a><br>
    <div class="panel panel-default">
        <div class="panel-body">
        Castro-Mondragon JA, Jaeger S, Thieffry D, Thomas-Chollier M#, van Helden J#. <i>"RSAT matrix-clustering: dynamic exploration and redundancy reduction of transcription factor binding motif collections."</i>, Nucleic Acid Research, 45:13 e119 (2017) <a href="https://www.ncbi.nlm.nih.gov/pubmed/28591841" target="_blank">[Pubmed]</a><a href="https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkx314" target="_blank">[Full text]</a>
        </div>
    </div>-->
</div>

<!-- ################################################################ -->
<!-- ### mandatory inputs ### -->
                <div class="bhoechie-tab-content">
<!-- title -->
<div class="panel panel-danger">
    <div class="panel-heading">Analysis Title <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Title that will be displayed at the top of the report page." data-original-title=""></i></div>

    <div class="panel-body">
        <div class="form-group">';

print $query->textfield(-id=>'html_title',
-name=>'html_title', -class=>'form-control',-placeholder=>'Provide a Title for this analysis ', -required=>'true',
                         -default=>$default{html_title}); 

print '
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
#                      -checked=>$default{quick},
#                      -label=>'');
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
    print $query->checkbox(-name=>'Heatmap',
                       -checked=>$default{heatmap},
                       -class=>form-check-input,
                        -label=>'');
print   'Heatmap</label></div>';

## Export the trees in Newick format
### By default trees are exported in JSON
print '<div class="form-check">
    <label class="form-check-label">';
    print $query->checkbox(-name=>'newick',
                       -checked=>$default{newick},
                       -class=>form-check-input,
                        -label=>'');
print   'Export the trees in Newick format</label></div>';

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
### Send results by email only
print "<p>\n";
&SelectOutput();
print " </div>
  </div> </div>";


################################################################
### Action buttons






















#### upload query classifcation file
print "<a href='help.compare-classes.html#upload_query_classes'>Query classification file</a><BR>";
print $query->filefield(-name=>'upload_query_classes',
			-default=>$default{upload_query_classes},
			-size=>30,
			-maxlength=>200);
print "<p>";

#### upload reference classifcation file
print "<a href='help.compare-classes.html#upload_ref_classes'>Reference classification file</a><BR>";
print $query->filefield(-name=>'upload_ref_classes',
			-default=>$default{upload_ref_classes},
			-size=>30,
			-maxlength=>200);
print "<p>";

#### table with all the statistics and thresholds
print "<h4>Return</h4>\n";

print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-classes.html#return_fields'>Fields</A> "]),

				### occurrences
				$query->td([$query->checkbox(-name=>'occ',
							     -checked=>$default{occ},
							     -label=>' Occurrences ')
					    ]),

				### Frequencies
				$query->td([$query->checkbox(-name=>'freq',
							     -checked=>$default{freq},
							     -label=>' Frequencies ')
					    ]),

				### Probabilities
				$query->td([$query->checkbox(-name=>'proba',
							     -checked=>$default{proba},
							     -label=>' Probabilities ')
					    ]),


				### Jaccard index
				$query->td([$query->checkbox(-name=>'jac',
							     -checked=>$default{jac},
							     -label=>' Jaccard index ')
					    ]),
				### Entropy
				$query->td([$query->checkbox(-name=>'entropy',
							     -checked=>$default{entropy},
							     -label=>' Entropy ')
					    ]),

				### Members
				$query->td([$query->checkbox(-name=>'members',
							     -checked=>$default{members},
							     -label=>' Members '),
					    ]),

			 ]
			)
		);

print "<h4>Thresholds</h4>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-classes.html#return_fields'>Fields</A> ",
					    " <A HREF='help.compare-classes.html#thresholds'>Lower<BR>Threshold</A> ",
					    " <A HREF='help.compare-classes.html#thresholds'>Upper<BR>Threshold</A> ",
					    ]),
				
				### Query class size
				$query->td([' Query size ',
					    $query->textfield(-name=>'lth_q',
							      -default=>$default{lth_q},
							      -size=>5),
					    $query->textfield(-name=>'uth_q',
							      -default=>$default{uth_q},
							      -size=>5),
					    ]),
				### Reference class size
				$query->td([' Reference size ',
					    $query->textfield(-name=>'lth_r',
							      -default=>$default{lth_r},
							      -size=>5),
					    $query->textfield(-name=>'uth_r',
							      -default=>$default{uth_r},
							      -size=>5),
					    ]),
				### Intersection size
				$query->td([' Intersection size ',
					    $query->textfield(-name=>'lth_qr',
							      -default=>$default{lth_qr},
							      -size=>5),
					    $query->textfield(-name=>'uth_qr',
							      -default=>$default{uth_qr},
							      -size=>5),
					    ]),
				### Significance 
				$query->td([' Significance ',
					    $query->textfield(-name=>'lth_sig',
							      -default=>$default{lth_sig},
							      -size=>5),
					    $query->textfield(-name=>'uth_sig',
							      -default=>$default{uth_sig},
							      -size=>5),
					    ]),

				### P-value 
				$query->td([' P-value ',
					    $query->textfield(-name=>'lth_pval',
							      -default=>$default{lth_pval},
							      -size=>5),
					    $query->textfield(-name=>'uth_pval',
							      -default=>$default{uth_pval},
							      -size=>5),
					    ]),

				### E-value 
				$query->td([' E-value ',
					    $query->textfield(-name=>'lth_eval',
							      -default=>$default{lth_eval},
							      -size=>5),
					    $query->textfield(-name=>'uth_eval',
							      -default=>$default{uth_eval},
							      -size=>5),
					    ]),
				### Jaccard index
				$query->td([' Jaccard index ',
					    $query->textfield(-name=>'lth_jac',
							      -default=>$default{lth_jac},
							      -size=>5),
					    $query->textfield(-name=>'uth_jac',
							      -default=>$default{uth_jac},
							      -size=>5),
					    ]),
				$query->td([' Mutual information ',
					    $query->textfield(-name=>'lth_mi',
							      -default=>$default{lth_mi},
							      -size=>5),
					    $query->textfield(-name=>'uth_mi',
							      -default=>$default{uth_mi},
							      -size=>5),
					    ]),

			 ]
			)
		);



################################################################
## sort key
print "<b><a href='help.compare-classes.html#sort_key'>Sort key </a></b>";
print  $query->popup_menu(-name=>'sort_key',
			  -Values=>['sig',
				    'E_val', 
				    'P_val',
				    'Jaccard index',
				    'Mutual information',
				    'names'
				    ],
			  -default=>$sequence_format);

################################################################
## population size
print "&nbsp"x8, "<b><a href='help.compare-classes.html#pop_size'>Population size </a></b>";
print $query->textfield(-name=>'pop_size',
			-default=>$default{pop_size},
			-size=>5);

### send results by email or display on the browser
print "<HR width=550 align=left>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
# print $query->start_multipart_form(-action=>"compare-classes_form.cgi");
# print "<TD><B>";
# print $query->hidden(-name=>'sort_key',-default=>"sig");
# print $query->submit(-label=>"DEMO");
# print "</B></TD>\n";
# print $query->end_form;


print "<TD><B><A HREF='help.compare-classes.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_compare-classes.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "<HR>";

print $query->end_html;

exit(0);


