#!/usr/bin/perl
#### redirect error log to a file
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("compare-features result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
@result_files = ();
$command = "$SCRIPTS/compare-features";
#$tmp_file_name = sprintf "compare-features.%s", &AlphaDate();
$prefix = "compare-features";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

#### read parameters ####
$parameters = " -v 1";

#### return a confusion table
# if ($query->param('return') eq "matrix") {
#     $parameters .= " -matrix"; 
# } 

### fields to return
$return_fields = " -return ";

### statistics
if ($query->param('stats')) {
    $return_fields .= "stats,";
} 

### intersection
if ($query->param('inter')) {
    $return_fields .= "inter,";
}

### differences
if ($query->param('diff')) {
	$return_fields .= "diff,";
}

### lower threshold on interection leength (size)
if ($query->param('inter_len') =~ /^\d+$/) {
  $inter_len = $query->param('inter_len');
  if (&IsNatural($inter_len)) {
    $parameters .= " -lth inter_len ".$inter_len;
  } else {
    &FatalError("Lower threshold on inter_len: $inter_len invalid value.");
  }
}

### lower threshold on interection leength (size)
if ($query->param('inter_cov') =~ /\d+/) {
  $inter_cov = $query->param('inter_cov');
  if (&IsReal($inter_cov)) {
    $parameters .= " -lth inter_cov ".$inter_cov;
  } else {
    &FatalError("Lower threshold on inter_cov: $inter_cov invalid value.");
  }
}


#### load the query feature file
$tmp_query_features = $tmp_file_path."_upload_query_features.tab";
push @result_files, "Query features", $tmp_query_features;
$uploaded_file = $query->param('upload_ref_features');
if ($uploaded_file) {
  $upload_query_features = $query->param('upload_query_features');
  if ($upload_query_features) {
    if ($upload_file =~ /\.gz$/) {
      $tmp_query_features .= ".gz";
    }
    $type = $query->uploadInfo($upload_query_features)->{'Content-Type'};
    open FEATURES, ">$tmp_query_features" ||
      &cgiError("Cannot store query feature file in temp dir.");
    while (<$upload_query_features>) {
      print FEATURES;
    }
    close FEATURES;
  }

## Pasted query features
}elsif ($query->param('featQ') =~/\S/) {
#  $tmp_query_features = "${TMP}/${tmp_file_name}_pasted_query_features.tab";
  open FEATURES, "> $tmp_query_features";
  print FEATURES $query->param('featQ');
  close FEATURES;
  &DelayedRemoval($tmp_query_features);
}else {
  &FatalError ("Please select the query feature file on your hard drive with the Browse button or paste features in the text area");
}
$parameters .= " -i $tmp_query_features";

## Load reference feature file
$tmp_ref_features = $tmp_file_path."_uploaded_ref_features.tab";
push @result_files, "Reference features", $tmp_ref_features;
$uploaded_file = $query->param('upload_ref_features');
if ($uploaded_file) {
  $upload_ref_features = $query->param('upload_ref_features');
  if ($upload_ref_features) {
    if ($upload_file =~ /\.gz$/) {
      $tmp_ref_features .= ".gz";
    }
    $type = $query->uploadInfo($upload_ref_features)->{'Content-Type'};
    open FEATURES, ">$tmp_ref_features" ||
      &cgiError("Cannot store expected frequency file in temp dir.");
    while (<$upload_ref_features>) {
      print FEATURES;
    }
    close FEATURES;
  }
} elsif ($query->param('featRef') =~/\S/) {
  #    $tmp_ref_features = "${TMP}/${tmp_file_name}_pasted_ref_features.tab";
  open FEATURES, "> $tmp_ref_features";
  print FEATURES $query->param('featRef');
  close FEATURES;
  &DelayedRemoval($tmp_ref_features);
}else {
  &FatalError ("Please select the reference feature file on your hard drive with the Browse button or paste features in the text area");
}
$parameters .= " -ref $tmp_ref_features";

$command .= " -iformat ".$query->param('feature_format');

$command .= " ".$return_fields." ".$parameters;

&ReportWebCommand($command);

$result_file = $tmp_file_path.".tab";
push @result_files, "Comparison result", $result_file;

if ($query->param('output') =~ /display/i) {

#    &PipingWarning();

    ### execute the command ###
    open RESULT, "$command |";

    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);

    &PrintURLTable(@result_files);
#    &PipingForm();
    print '<HR SIZE=3>';

} else {
    &EmailTheResult("$command ", $query->param('user_email'), $tmp_file_path);
}

print $query->end_html;

exit(0);


