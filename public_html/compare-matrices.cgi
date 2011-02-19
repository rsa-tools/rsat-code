#!/usr/bin/perl

############################################ imports
### redirect error log to a file
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

################################################################
## configuration
$command = "$ENV{RSAT}/perl-scripts/compare-matrices";
$output_directory = sprintf "compare-matrices.%s", &AlphaDate();

## We need to create the output directory before starting
## compare-matrices, since it will generate multiple output files.
$output_path = $TMP."/".$output_directory;
$output_path =~ s|\/\/|\/|g;
system("mkdir -p $output_path");
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/perl-scripts" . ":$ENV{RSAT}/python-scripts";
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/bin";

$output_prefix = $output_path."/compare-matrices";

################################################################
## result page header
### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("compare-matrices result", "results");
&ListParameters() if ($ENV{rsat_echo} >=0);

### update log file
&UpdateLogFile();

################################################################
## command line paramters
### read parameters
$parameters = " -v 1";


#### read parameters ####
local $parameters = " -v 1";


################################################################
## Matrix input format
local $query_matrix_format = lc($query->param('matrix_format'));
($query_matrix_format) = split (/\s+/, $query_matrix_format);
$parameters .= " -format1 ".$query_matrix_format;

################################################################
#### Query matrix file
$matrix_file = $output_prefix."_query_matrices.".$query_matrix_format;
if ($query->param('matrix')) {
    open MAT, "> $matrix_file";
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
    $parameters .= " -file1 $matrix_file";
} else {
    &RSAT::error::FatalError('You did not enter any data in the matrix box');
}
push @result_files, ("Input file",$matrix_file);
push @result_files, ("Result file",$result_file);


################################################################
## Pseudo-counts
if (&IsReal($query->param('pseudo_counts'))) {
  $parameters .= " -pseudo ".$query->param('pseudo_counts');
} else {
  &FatalError("Pseudo-count should be a real number");
}
if ($query->param('pseudo_distribution') eq "equi_pseudo") {
  $parameters .= " -equi_pseudo ";
}

################################################################
## Background model method
&SetBackgroundModel();

################################################################
## bg_pseudo
if (&IsReal($query->param('bg_pseudo'))) {
  $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
}

################################################################
## Reference motifs
if ($query->param('db_choice') eq "custom") {
  ## Upload custom reference motif file
  local $ref_motif_file = "${TMP}/${tmp_file_name}_ref_motif_file.txt";
  local $upload_ref_motif_file = $query->param('upload_ref_motif_file');
  if ($upload_ref_motif_file) {
    if ($upload_ref_motif_file =~ /\.gz$/) {
      $ref_motif_file .= ".gz";
    }
    local $type = $query->uploadInfo($upload_ref_motif_file)->{'Content-Type'};
    open REF_MOTIF_FILE, ">$ref_motif_file" ||
      &cgiError("Cannot store reference motif fie file in temp dir.");
    while (<$upload_ref_motif_file>) {
      print REF_MOTIF_FILE;
    }
    close REF_MOTIF_FILE;
    $parameters .= " -file2 $ref_motif_file";
  } else {
    &RSAT::error::FatalError("You did not specify the custom matrix file with the Browse button");
  }
} else {
  my ($mat_db_params, @selected_db) = &GetMatrixDBfromBox();
  if (scalar(@selected_db) > 0) {
    $parameters .= " -file2 ".$selected_db[0];
  }
}

### other default parmaters
$parameters .= " -strand DR";

### output directory
$parameters .= " -o $output_path";

## Report the full command before executing
print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >=1);

################################################################
## display or send result
$index_file = $output_directory."/".$output_prefix."_synthesis.html";
my $mail_title = join (" ", "[RSAT]", "compare-matrices", &AlphaDate());
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere",$index_file, title=>$mail_title ,no_email=>1);
} else {

&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);
# $debug = "$command $parameters 2> $TMP/log.txt";
# print $debug;
# `$debug`;
}

################################################################
## result page footer
print $query->end_html;

exit(0);

