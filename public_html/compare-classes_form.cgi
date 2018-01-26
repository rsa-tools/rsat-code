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
#$default{pipe} = "";

$default{metric} = 'QR';
$default{occ} = 1;
$default{sort} = 1;
$default{proba} = 1;
$default{jac} = 1;

#$default{members} = "";
#$default{sort_key} = "sig";
#$default{pop_size} = "auto";

$default{lth_q} = 1;
$default{uth_q} = "none";
$default{lth_r} = 1;
$default{uth_r} = "none";
$default{lth_qr} = 1;
$default{uth_qr} = "none";
$default{lth_sig} = 0;
$default{uth_sig} = "none";

# TOBEDONE: check which tools might produce output pipeable to this form

### replace defaults by parameters from the cgi call, if defined
#foreach $key (keys %default) {
#    if ($query->param($key)) {
#        $default{$key} = $query->param($key);
#    }
#    if ($query->param($key) =~ /checked/i) {
#        $checked{$key} = "CHECKED";
#    }
#}

&ListParameters() if ($ENV{rsat_echo} >= 2);

### print the form as in matrix-clustering_form.cgi
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
                        <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>Compare classes
                    </a>
                    <a href="#" class="list-group-item text-center">
                        <h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Main input
                    </a>
                    <!--
                    <a href="#" class="list-group-item text-center">
                        <h4 class="glyphicon"><i class="fa fa-tags fa-2x"></i></h4><br/>Optional input
                    </a> 
                    -->
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
    <a target="_blank" href="http://jacques.van-helden.perso.luminy.univ-amu.fr/ ">Jacques van Helden</a> with help from Joseph Tran and Bruno Contreras-Moreira.<br>
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
    <!--<span class="fa-stack fa-lg">
        <i class="fa fa-pencil fa-stack-1x"></i>
    </span>
    Cite the publication: <a href="https://twitter.com/rsatools" target="_blank"></a><br>
    <div class="panel panel-default">
        <div class="panel-body">
        # citation should go here
        </div>
    </div>-->
            </div>

<!-- ################################################################ -->
<!-- ### main input ### -->
    
            <div class="bhoechie-tab-content">

                <!-- query classes -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Query classes <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="A tab-delimited text file containing the description of query classes." data-original-title=""></i></div>
                    <div class="panel-body">
                        <div class="form-group">';

                        print $query->textarea( -id=>'classesQ',-name=>'classesQ',-rows=>5,-cols=>60, -required=>'true',
                            -placeholder=>'Paste here your query classes, or select a file to upload below',
                            -default=>$default{query_classes});
                        print "<br><b>Or</b> select a file to upload<br>\n";
                        print $query->filefield(-name=>'Qclass_file',-default=>'',-size=>40);

print '                 </div>
                    </div>
                </div>

                <!-- reference classes -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Reference classes <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="A tab-delimited text file containing the description of reference classes." data-original-title=""></i></div>
                    <div class="panel-body">
                        <div class="form-group">';

                        print $query->textarea( -id=>'classesR',-name=>'classesR',-rows=>5,-cols=>60, -required=>'true',
                            -placeholder=>'Paste here your reference classes, or select a file to upload below',
                            -default=>$default{ref_classes});
                        print "<br><b>Or</b> select a file to upload<br>\n";
                        print $query->filefield(-name=>'Rclass_file',-default=>'',-size=>40);

print '                 </div>
                    </div>
                </div>

                <!-- output format -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Format</div>
                    <div class="panel-body">
                        <div class="form-group">';

                        my %output_labels = (
                            'classes',' Pairwise class comparison table',
                            'matrix',' Matrix with reference classes as rows and query classes as columns' );

                        print $query->radio_group( -name => 'outformat',-values  => ['classes','matrix'],-default => 'classes',
                            -labels=>\%output_labels)."<br>";

print '                 </div>
                    </div>
                </div>
            </div>

<!-- ################################################################ 
### optional input ### 

            <div class="bhoechie-tab-content">

                <div class="panel panel-danger">
                    <div class="panel-heading">Score column <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Column of the input files containing a score associated to each member. Must be valid for both query and reference classes. It is used for some metrics like the dot product." data-original-title=""></i></div>
                    <div class="panel-body">
                        <div class="form-group">';
                        print $query->textfield(-id=>'score_col',-name=>'score_col',-size=>10,-placeholder=>'optional') .'
                        </div>
                    </div>
                </div>

                <div class="panel panel-danger">
                    <div class="panel-heading">Type of comparison</div>
                    <div class="panel-body">
                        <div class="form-group"> -->';

#                        my %self_compa_labels = ( 
#                            'off',' Cross-compare query classes to reference classes',
#                            'on',' Self-compare query classes to query classes' );
#                        print $query->radio_group( -name => 'self_compa',
#                            -values  => ['off', 'on'],-default => 'off',
#                            -labels=>\%self_compa_labels)."<br>";
#print '                 </div>
#                    </div>
#                </div>
#            </div>

print '
<!-- ################################################################-->
<!-- ### advanced output options  ###-->

            <div class="bhoechie-tab-content">

                <!-- comparison type -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Type of comparison</div>
                    <div class="panel-body">
                        <div class="form-group">';
                        my %self_compa_labels = (
                            'off',' Cross-compare query classes to reference classes',
                            'on',' Self-compare query classes to query classes' );
                        print $query->radio_group( -name => 'self_compa',
                            -values  => ['off', 'on'],-default => 'off',
                            -labels=>\%self_compa_labels)."<br>";

                        # commented out as they seem confusing,Bruno jan2018
                        #print $query->checkbox(-name=>'distinct',-checked=>1,-value=>'on',
                        #   -label=>'Prevent self-comparison of classes')."<br>";
                        #print $query->checkbox(-name=>'triangle',-checked=>1,-value=>'on',
                        #   -label=>'Prevent reciprocal comparison of classes, only applies to self');   

print '                 </div>
                    </div>
                </div>

                <!-- score column -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Score column <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Column of the input files containing a score associated to each member. Must be valid for both query and reference classes. It is used for some metrics like the dot product." data-original-title=""></i></div>
                    <div class="panel-body">
                        <div class="form-group">';
                        print $query->textfield(-id=>'score_col',-name=>'score_col',-size=>10,-placeholder=>'optional') .'
                        </div>
                    </div>
                </div>

                <!-- matrix metric  -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Metric in matrix output</div>
                    <div class="panel-body">
                        <div class="form-group">';
my %metric_labels = (
'qr',' Intersection',
'sig',' Significance',
'jac_sim',' Jaccard similarity',
'dotprod',' Dot product of score column',
'eval',' E-value',
'pval',' P-value',
'mi',' Mutual information ');

                        print $query->popup_menu(-id=>'matrix_metric', -name=>'matrix_metric',
                            -Values=>['qr','sig','jac_sim','dotprod','eval','pval','mi'],
                            -class=>'form-control',
                            -default=>$default{metric},
                            -labels=>\%metric_labels);
print '                 </div>
                    </div>
                </div>

                <!-- classes output fields -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Return fields of pairwise class comparison</div>
                    <div class="panel-body">
                        <div class="form-group">';

print $query->checkbox(-name=>'occ',-checked=>$default{'occ'},-value=>'on',-label=>'Occurrences').'<br>';
print $query->checkbox(-name=>'freq',-checked=>0,-value=>'on',-label=>'Frequencies').'<br>';
print $query->checkbox(-name=>'proba',-checked=>$default{'proba'},-value=>'on',-label=>'Hypergeometric probability').'<br>';
print $query->checkbox(-name=>'sort',-checked=>$default{'sort'},-value=>'on',-label=>'Sorting criterion').'<br>';
print $query->checkbox(-name=>'jac_sim',-checked=>$default{'jac'},-value=>'on',-label=>'Jaccard similarity').'<br>';
#print $query->checkbox(-name=>'sor_sim',-checked=>0,-value=>'on',-label=>'Sorensen similarity').'<br>';
print $query->checkbox(-name=>'dotprod',-checked=>0,-value=>'on',-label=>'Dot product, relevant if a score column is specified').'<br>';
print $query->checkbox(-name=>'entropy',-checked=>0,-value=>'on',-label=>'Entropy').'<br>';
print $query->checkbox(-name=>'members',-checked=>0,-value=>'on',-label=>'Class members, might generate large result files');

print '                 </div>
                    </div>
                </div>

                <!-- thresholds -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Thresholds for the pairwise class comparison table</div>
                    <div class="panel-body">
                        <div class="form-group">';

                        &PrintThresholdTableForm();

print '                 </div>
                    </div>
                </div>
            </div>

<!--################################################################-->
<!--### run & reset ###-->

            <div class="bhoechie-tab-content">

                <!-- results delivery  -->
                <div class="panel panel-danger">
                    <div class="panel-heading">Results delivery options</div>
                    <div class="panel-body">
                        <div class="form-group">';

                        &SelectOutput();

                        print $query->submit(-label=>"GO", -class=>"btn btn-success", -type=>"button");
                        print " ";
                        print $query->reset(-id=>"reset",-class=>"btn btn-warning", -type=>"button");

print "                 </div>
                    </div>
                </div>
            </div>

</div></div></div></div>";

print $query->end_form;

print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr' class='col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container'></div>";


################################################################
### Demo data

my $demo_fileQ = "demo_files/gavin_mcl_clusters_inf2.1.tab";
my $demo_fileR = "demo_files/mips_complexes.tab";
my ($demoQ,$demoR);

open(FILEQ, $demo_fileQ);
while(my $row = <FILEQ>){
    $demoQ .= $row;
}
close(FILEQ);

open(FILER, $demo_fileR);
while(my $row = <FILER>){
    $demoR .= $row;
}
close(FILER);

print '<script>
function setDemo(demoQ, demoR){
    $("#reset").trigger("click");
    descr = "<H4>Demonstration:</H4>\n \
    <blockquote class =\'blockquote text-justify small\'>\
    This demo consists on the comparison between protein clusters \
    obtained after application of the <a href=\'http://micans.org/mcl\' \
    target=\'top\'>MCL</a> clustering algorithm to the <a target=\'_blank\' \
    href=\'https://www.ncbi.nlm.nih.gov/pubmed/16429126\'>Gavin et al \
    (2006)</a> interaction network and the complexes annotated in the \
    <a target=\'_blank\' href=\'http://mips.gsf.de\'>MIPS</a> database. \
    Check the panel <b>Main input</b> and <b>Run analysis</b></blockquote>";

    demo_descr.innerHTML = descr;
    html_title.value = "\'clusters of protein complexes in Gavin et al dataset\'";
    collection_label.value = "\'yeast protein complexes\'";
    matrix1.value = demo_1_matrix;
    demo.value = descr;
}
</script>';

print ' <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container">

<button type="button" class="btn btn-info" onclick="setDemo1('. "'$demo_1_matrices'" .')">DEMO (one collection)</button> ';





#print "</div> </div> </div>";

print $query->end_html;

exit(0);

## Print table of supported return fields of pairwise comparison tables
sub PrintThresholdTableForm {

    my @vars = qw( q r qr sig eval pval jac_sim mi dotprod );

    my %descriptions = (
        'q',' Query occurrences (Q)', # not sure if this is equal to occ?
        'r',' Reference occurrences (R)',
        'qr',' Intersection occurrences (QR)',
        'sig',' Significance',
        'pval',' P-value of the intersection, hypergeometric function',
        'eval',' E-value = P-value * nb_tests',
        'jac_sim',' Jaccard similarity = (Q and R)/(Q or R) ',
        'mi',' Mutual information of class Q and R',
        'dotprod',' Dot product (if score column is set)' );

    &ThresholdsDiv("Thresholds of return fields",
        "help.compare_classes.html#thresholds",
        \@vars,
        \%descriptions);
}

## Display a collapsable div with selectable scores and thresholds
sub ThresholdsDiv {
  my ($title, $help_file, $field_ref, $field_descr_ref) = @_;
  print "<p><fieldset class='form-group'>\n<b>".$title."</b>";
  &FieldsThresholdsTableMC($help_file, $field_ref, $field_descr_ref);
    print '</fieldset><p/>';
}

## Display a table with checkboxes and thresholds for a set of return fields
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
    my ($lth,$uth);
    if(defined($default{'lth_'.$field})){
        $lth = $default{'lth_'.$field};
    }
    else{ $lth = "none" }
    if(defined($default{'uth_'.$field})){
        $uth = $default{'uth_'.$field};
    }
    else{ $uth = "none" }

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





