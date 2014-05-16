#!/usr/bin/perl

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

################################################################
## result page header
### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("matrix-clustering result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

#
## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

################################################################
## Output paths
$command = "$ENV{RSAT}/perl-scripts/matrix-clustering";

$output_prefix = "matrix-clustering";
$output_path = &RSAT::util::make_temp_file("",$output_prefix, 1); $output_dir = &ShortFileName($output_path);

## We need to create the output directory before starting
## matrix-clustering, since it will generate multiple output files.
#&RSAT::util::CheckOutDir($output_dir);
system("rm -f $output_path; mkdir -p $output_path"); ## We have to delete the file created by &make_temp_file() to create the directory with same name

################################################################
## Command line paramters
$parameters .= " -v 1";

## Read parameters
local $parameters = " -v 1";

################################################################
## Matrix input format
local $query_matrix_format = lc($query->param('matrix_format'));
($query_matrix_format) = split (/\s+/, $query_matrix_format);
$parameters .= " -format ".$query_matrix_format;

################################################################
#### Query matrix file
$matrix_file = $output_path."/".$output_prefix."_query_matrices.".$query_matrix_format;
if ($query->param('matrix')) {
  open MAT, "> $matrix_file";
  print MAT $query->param('matrix');
  close MAT;
  &DelayedRemoval($matrix_file);
  $parameters .= " -i $matrix_file";
} else {
  &RSAT::error::FatalError('You did not enter any data in the matrix box');
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
## Specify the thresholds on all parameters for compare-matrices
my @selected_output_fields = ();
my $thresholds = "";
foreach my $field (@supported_output_fields) {
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
$parameters .= $thresholds;

## Add selected output fields
push @selected_output_fields, qw(
				 matrix_id
				 matrix_name
				 width
				 strand
				 offset
				 consensus
				);

my $selected_output_fields = join (",", @selected_output_fields);
$parameters .= " -return ".$selected_output_fields;

################################################################
## Output formats
$parameters .= " -export newick";
$parameters .= " -d3_base file";
$parameters .= " -labels name,consensus";



################################################################
## Output file
$parameters .= " -o ".$output_path."/".$output_prefix;

## Add an error-log file for matrix-clustering
  $err_file = $output_path."/".$output_prefix."_err.txt";
  $parameters .= " >& ".$err_file;

## Report the full command before executing
  &ReportWebCommand($command." ".$parameters);

################################################################
## Display or send result by email
  $index_file = $output_path."/".$output_prefix."_index.html";
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

