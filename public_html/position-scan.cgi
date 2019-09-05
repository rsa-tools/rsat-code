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

################################################################
## Result page header

## Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("position-scan result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

################################################################
## Output paths
$command = $ENV{RSAT}."/perl-scripts/position-scan";
$tasks = "all";

$output_prefix = "position-scan";
$output_path = &RSAT::util::make_temp_file("",$output_prefix, 1); $output_dir = &ShortFileName($output_path);

## We need to create the output directory before starting
## position-scan, since it will generate multiple output files.
#&RSAT::util::CheckOutDir($output_dir);
system("rm -f $output_path; mkdir -p $output_path"); ## We have to delete the file created by &make_temp_file() to create the directory with same name

################################################################
## Command line paramters
local $parameters .= " -v 1";
$parameters .= " -task ".$tasks;

################################################################
## Matrix input format
local $query_matrix_format = lc($query->param('matrix_format'));
($query_matrix_format) = split (/\s+/, $query_matrix_format);
$parameters .= " -matrix_format ".$query_matrix_format;

################################################################
#### Query matrix file
local $matrix_file = &GetMatrixFile($output_path."/".$output_prefix."_query_matrices.".$query_matrix_format);

$parameters .= " -matrix ".$matrix_file;

push @result_files, ("Input matrix file",$matrix_file);
push @result_files, ("Result file",$result_file);

################################################################
## Agglomeration rule for hierarchical clustering
local $class_interval  = lc($query->param('class_interval'));
if ($class_interval) {
  $parameters .= " -bin ".$class_interval;
}

###############
## Add title
local $title = lc($query->param('html_title'));
if($title){
    $title =~ s/\s+/_/g;
    $parameters .= " -title '".$title."'";
}

################################
## Add motif collection label
local $pval = lc($query->param('thresh_value'));
if($pval){
    $parameters .= " -pval ".$pval;
}

##############################
## Sequence file and format
($sequence_file, $sequence_format) = &GetSequenceFile();
push @result_files, ("Input sequences", $sequence_file);
$parameters .= " -seq ".$sequence_file." -seq_format ".$sequence_format;

################################################################
## Background model method
################################################################
## Background model method
local $bg_method = $query->param('bg_method');
if ($bg_method eq "from_matrix") {
    
}   elsif ($bg_method eq "bginput") {
    $parameters .= " -bginput";
    $parameters .= " -markov ".$markov_order;

}   elsif ($bg_method eq "bgfile") {
    ## Select pre-computed background file in RSAT genome directory
    local $organism_name = $query->param("organism");
    local $noov = "ovlp";
    local $background_model = $query->param("background");
    #local $oligo_length = 1;
    local $oligo_length = $markov_order + 1;
    $bg_file = &ExpectedFreqFile($organism_name,
				 $oligo_length, $background_model,
				 noov=>$noov, str=>"-1str");
    $parameters .= " -bgfile ".$bg_file.".gz";
    $parameters .= " -bg_format ".'oligo-analysis';
    
} elsif ($bg_method =~ /upload/i) {
    ## Upload user-specified background file
    local $bgfile = $tmp_file_name."_bgfile.txt";
    local $upload_bgfile = $query->param('upload_bgfile');
    if ($upload_bgfile) {
	if ($upload_bgfile =~ /\.gz$/) {
	    $bgfile .= ".gz";
	}
	local $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
	open BGFILE, ">$bgfile" ||
	    &cgiError("Cannot store background file in temp dir.");
	while (<$upload_bgfile>) {
	    print BGFILE;
	}
	close BGFILE;
	$parameters .= " -bgfile $bgfile";
	$parameters .= " -bg_format ".$query->param('bg_format');
	
    } 
} elsif ($bg_method =~ /url/i) {
    ## Retrieve user-specified URL for background file
    my $url = $query->param('bgmodel_url');
    
    &RSAT::message::Info("Fetching background from URL ".$url) if ($ENV{rsat_echo} >= 1);
      my $bgmodel = "";
    local $bgfile = $tmp_file_name."_bgfile.txt";
    
    if (open BGM, ">$bgfile") {
	$bg = get($url);
	if ($bg =~ /\S/) {
	    print BGM $bg;
	    close BGM;
	} else {
	    &RSAT::error::FatalError("No background model could be downloaded from the URL ".$url);
	}
	
    }
    
    close BGFILE;
	$parameters .= " -bgfile $bgfile";
	$parameters .= " -bg_format ".$query->param('bg_format');
    
} else {
    &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
}


######################
## Origin selection
$origin = $query->param('origin');
# if ($heatmap) {
#     $return_fields .= ",heatmap";
# }


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
$index_file = $output_path."/".$output_prefix."_report.html";
my $mail_title = join (" ", "[RSAT]", "position-scan", &AlphaDate());
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere", "", title=>$mail_title, index=>$index_file, no_email=>1);
} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), "", title=>$mail_title,index=>$index_file);
}

################################################################
## Result page footer
print $query->end_html;

exit(0);

