#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("convert-features result", 'results');


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/convert-features -v 1";
$prefix = "convert-features";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
$tmp_file_name = sprintf "convert-features.%s", &AlphaDate();
@result_files = ();

#### read parameters ####
my $parameters;

################################################################
## feature input format
my $input_format = lc($query->param('feature_format'));
$parameters .= " -from ".$input_format;


################################################################
## feature output format
my $output_format = lc($query->param('output_format'));
$parameters .= " -to ".$output_format;


################################################################
#### Feature from input box
my $input_file = $tmp_file_path.".".$input_format;
push @result_files, ("Input features ($input_format)",$input_file);

if ($query->param('feature')){
  open FEAT, "> $input_file";
  print FEAT $query->param('feature');
  close FEAT;
} else  {
  ## Upload user-specified  file
  my $upload_file = $query->param('uploaded_file');
  if ($upload_file) {
    if ($upload_file =~ /\.gz$/) {
      $input_file .= ".gz";
    }
    my $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open FEAT, ">$input_file" ||
      &cgiError("Cannot store feature file in temp dir.");
    while (<$upload_file>) {
      print FEAT;
    }
    close FEAT;
  } else {
    &FatalError ("If you want to upload a file, you should specify the location of this file on your hard drive with the Browse button");
  }
}
&DelayedRemoval($input_file);
$parameters .= " -i $input_file";

###########################
  ## coord   file
my $coord_file = $tmp_file_path."_coord.bed";


if ($query->param('bed_coord')){
  open FEAT, "> $coord_file";
  print FEAT $query->param('bed_coord');
  close FEAT;
  push @result_files, ("Input Bed file with genomic coordinates",$coord_file);
  &DelayedRemoval($coord_file);  
  $parameters .= " -coord ".$coord_file;
} else  {
  ## Upload user-specified  file
  my $upload_file = $query->param('bed_file');
  if ($upload_file) {
    if ($upload_file =~ /\.gz$/) {
      $input_file .= ".gz";
    }
    my $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open FEAT, ">$coord_file" ||
      &cgiError("Cannot store coordinate file in temp dir.");
    while (<$upload_file>) {
      print FEAT;
    }
    close FEAT;
    push @result_files, ("Input Bed file with genomic coordinates",$coord_file);
    &DelayedRemoval($coord_file);  
 	$parameters .= " -coord ".$coord_file;
  }
}


my $output_file = $tmp_file_path.".".$output_format;
push @result_files, ("Output features ($output_format)",$output_file);


&ReportWebCommand($command." ".$parameters);

### execute the command ###
if ($query->param('output') eq "display") {
#    &PipingWarning();

 ## prepare figures
    ### prepare data for piping
    open RESULT, "$command $parameters |";
    print '<H4>Result</H4>';
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $output_file, 1);
    close(RESULT);
#    print '<PRE>';
#    while (<RESULT>) {
#      print $_;
#    }
#    print '</PRE>';
#    close(RESULT);

    &PrintURLTable(@result_files);
#    &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $output_file);
}
print $query->end_html;

exit(0);

