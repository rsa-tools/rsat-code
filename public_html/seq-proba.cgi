#!/usr/bin/perl
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/seq-proba -v 1";
$tmp_file_name = sprintf "seq-proba.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("seq-proba result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();

#### read parameters ####
$parameters = "";

################################################################
## sequence file
($in_sequence_file,$sequence_format) = &GetSequenceFile();

#### parameters
$parameters .= " -i $in_sequence_file -seq_format $sequence_format";
&DelayedRemoval("$in_sequence_file");


################################################################
## Background model method
my $bg_method = $query->param('bg_method');
if ($bg_method eq "bgfile") {
  ## Select pre-computed background file in RSAT genome directory
  my $organism_name = $query->param("organism");
  my $noov = "ovlp";
  my $background_model = $query->param("background");
  my $markov_order = $query->param('markov_order');
  my $oligo_length = $markov_order + 1;
  $bg_file = &ExpectedFreqFile($organism_name,
			       $oligo_length, $background_model,
			       noov=>$noov, str=>"-1str");
  $parameters .= " -bgfile ".$bg_file;

  } elsif ($bg_method =~ /upload/i) {
    ## Upload user-specified background file
    my $bgfile = "${TMP}/${tmp_file_name}_bgfile.txt";
    my $upload_bgfile = $query->param('upload_bgfile');
    if ($upload_bgfile) {
      if ($upload_bgfile =~ /\.gz$/) {
	$bgfile .= ".gz";
      }
      my $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
      open BGFILE, ">$bgfile" ||
	&cgiError("Cannot store background file in temp dir.");
      while (<$upload_bgfile>) {
	print BGFILE;
      }
      close BGFILE;
      $parameters .= " -bgfile $bgfile";
      $parameters .= " -bg_format ".$query->param('bg_format');
    } else {
      &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
    }
		
  } else {
    &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
  }

## Return fields
@return_fields = qw(id proba_b log_proba len seq detail);
foreach my $field (@return_fields) {
  if ($query->param($field)) {
    $parameters .= " -return ".$field;
  }
}


print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >= 1);

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

#    $ENV{RSA_OUTPUT_CONTEXT} = "text";


    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$command $parameters |";

    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

#} elsif 
#    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

exit(0);


