#!/usr/bin/env perl

############################################ 
## redirect error log to a file
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
### redirect error log to a file
BEGIN {
  $ERR_LOG = "/dev/null";
#    $ERR_LOG = "/tmp/RSA_ERROR_LOG.txt";
  use CGI::Carp qw(carpout);
  open (LOG, ">> $ERR_LOG")
      || die "Unable to redirect log\n";
  carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
require RSAT::util;
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

################
## Restrict the number of input matrices treated on the Web server.
local $max_matrices = 300;

################################################################
## Result page header

## Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("matrix-clustering result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();


################################################################
## Output paths
$command = $ENV{RSAT}."/perl-scripts/matrix-clustering";
$return_fields = "-return json";

$output_prefix = "matrix-clustering";
$output_path = &RSAT::util::make_temp_file("",$output_prefix, 1); $output_dir = &ShortFileName($output_path);

## We need to create the output directory before starting
## matrix-clustering, since it will generate multiple output files.
#&RSAT::util::CheckOutDir($output_dir);
system("rm -f $output_path; mkdir -p $output_path"); ## We have to delete the file created by &make_temp_file() to create the directory with same name

################################################################
## Command line parameters
local $parameters .= " -v 1";
$parameters .= " -max_matrices ".$max_matrices;

############################################
## Collection Motif 1 (Mandatory)
############################################

## Matrix input format
local $query_matrix_format_1 = lc($query->param('matrix_format1'));

#### Query matrix file
local $matrix_file_1 = &GetMatrixFile($output_path."/".$output_prefix."_query_matrices.".$query_matrix_format_1,1);

## Add motif collection label
local $collection_1_label = lc($query->param('collection_1_label'));
if($collection_1_label){
    $collection_1_label =~ s/\s+/_/g;
} else {
    $collection_1_label = "Collection_1";
}

$parameters .= " -matrix $collection_1_label $matrix_file_1 $query_matrix_format_1";

######################
## Add collection 2 ##
######################

## check if second collection has been provided
if ($query->param('uploaded_file_matrix2') || ($query->param('matrix2') =~ /\S/)) {

## Matrix input format
local $query_matrix_2_format = lc($query->param('matrix_format2'));

#### Query matrix file
local $matrix_file_2 = &GetMatrixFile($output_path."/".$output_prefix."_second_matrices.".$query_matrix_format_2,2);

## Add motif collection label
local $collection_2_label = lc($query->param('collection_2_label'));
if($collection_2_label){
    $collection_2_label =~ s/\s+/_/g;
} else {
    $collection_2_label = "Collection_2";
}

    $parameters .= " -matrix $collection_2_label $matrix_file_2 $query_matrix_2_format";
}

######################
## Add collection 3 ##
######################

## Matrix input format
local $query_matrix_3_format = lc($query->param('matrix_format3'));

## check if third collection has been provided
if ($query->param('uploaded_file_matrix3') || ($query->param('matrix3') =~ /\S/)) {

#### Query matrix file
local $matrix_file_3 = &GetMatrixFile($output_path."/".$output_prefix."_third_matrices.".$query_matrix_3_format,3);

## Add motif collection label
local $collection_3_label = lc($query->param('collection_3_label'));
if($collection_3_label){
    $collection_3_label =~ s/\s+/_/g;
} else {
    $collection_3_label = "Collection_3";
}

    $parameters .= " -matrix $collection_3_label $matrix_file_3 $query_matrix_3_format";
}

push @result_files, ("Input file",$matrix_file);
push @result_files, ("Result file",$result_file);


################################################################
## Agglomeration rule for hierarchical clustering
local $hclust_method  = lc($query->param('hclust_method'));
if ($hclust_method) {
  $parameters .= " -hclust_method ".$hclust_method;
}


################################################################
## Merge operator
local $merge_stat  = lc($query->param('merge_stat'));
if ($merge_stat) {
  $parameters .= " -calc ".$merge_stat;
}


###############
## Add title
local $title = lc($query->param('html_title'));
if($title){
    $title =~ s/\s+/_/g;
    $parameters .= " -title '".$title."'";
}

############################
## Add the metric used to 
## cluster the motifs
local $metric_tree = $query->param('metric');
if($metric_tree){
    $parameters .= " -metric_build_tree '".$metric_tree."'";
}


################################################################
## Specify the thresholds on all parameters for compare-matrices
my @main_threshold_fields = qw(w
			 cor
			 Ncor
			);

my @threshold_fields = qw(logoDP
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
my $thresholds = "";

## main thresholds compatible with the quick version
foreach my $field (@main_threshold_fields) {
  ## Selected field
  if ($query->param('return_'.$field)) {
    push @selected_output_fields, $field;
  }
  
    ## Lower threshold
  my $lth = $query->param('lth_'.$field);
  if (&IsReal($lth)) {
    $thresholds .= " -lth ".$field." ".$lth;
  }
  
  ## Upper threshold
  my $uth = $query->param('uth_'.$field);
  if (&IsReal($uth)) {
    $thresholds .= " -uth ".$field." ".$uth;
  }
}

## other thresholds not compatible
my $forbidden = 0;
foreach my $field (@threshold_fields) {
  ## Selected field
  if ($query->param('return_'.$field)) {
    push @selected_output_fields, $field;
    $forbidden = 1;
  }
  
  ## Lower threshold
  my $lth = $query->param('lth_'.$field);
  if (&IsReal($lth)) {
    $thresholds .= " -lth ".$field." ".$lth;
  }
  
  ## Upper threshold
  my $uth = $query->param('uth_'.$field);
  if (&IsReal($uth)) {
    $thresholds .= " -uth ".$field." ".$uth;
  }
}
$parameters .= $thresholds;


## Run the C code compare-matrices-quick if the thresholds allow for it
if ($forbidden == 0) {
    $parameters .= " -quick";
}


##############################
## Add options from toolbox

## Heatmap selection
$heatmap = $query->param('heatmap');
if ($heatmap) {
    $return_fields .= ",heatmap";
}


## Export newick selection
#$newick = $query->param('newick');
#if ($newick) {
#    $return_fields .= ",newick";
#}


## Random permutation
$rand = $query->param('random');
if ($rand) {
    $parameters .= " -rand";
}

## Heatmap selection
$radial = $query->param('radial');
if ($radial) {
    $parameters .= " -radial_tree_only";
}

## Insert labels
#my @labs = ();
my $lab = "name";
#$label_id = $query->param('label_id');
#$label_name = $query->param('label_name');
#$label_consensus = $query->param('label_consensus');


#if($label_name){
#    push(@labs, "name");
#}

#if($label_id){
#    push(@labs, "id");
#}

#if($label_consensus){
#    push(@labs, "consensus");
#}

#$lab = join(",", @labs);


$parameters .= " -label_in_tree ".$lab;

## Insert fields to return 
$parameters .= " ".$return_fields." ";

################################################################
## Output file
$parameters .= " -o ".$output_path."/".$output_prefix;

## Add an error-log file for matrix-clustering
$err_file = $output_path."/".$output_prefix."_err.txt";
$parameters .= " 2> ".$err_file;

## Report the full command before executing
&ReportWebCommand($command." ".$parameters, $err_file);

################################################################
## Display or send result by email
$index_file = $output_path."/".$output_prefix."_SUMMARY.html";
my $mail_title = join (" ", "[RSAT]", "matrix-clustering", &AlphaDate());
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere", "", title=>$mail_title, index=>$index_file, no_email=>1);
} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), "", title=>$mail_title,index=>$index_file);
}

################################################################
## Result page footer
print $query->end_html;

exit(0);

