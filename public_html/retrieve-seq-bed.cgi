#!/usr/bin/env perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
#BEGIN {
#    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
#    use CGI::Carp qw(carpout);
#    open (LOG, ">> $ERR_LOG")
#	|| die "Unable to redirect log\n";
#    carpout(*LOG);
#}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = "$SCRIPTS/retrieve-seq-bed";
$prefix="retrieve-seq-bed";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);

@result_files = ();


## Read the CGI query
$query = new CGI;

## Print the header
&RSA_header("retrieve-seq-bed result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

############################################################
## Read parameters
$parameters = " -v 1 ";

################################################################
## Organism
my $organism = $query->param("organism");
&RSAT::error::FatalError("No organism was specified") unless ($organism);
$parameters .= " -org ".$organism;

############################################################
## Genomic fragments

## Coordinate file
my $coordinate_format = "bed"; ## Later I will add an option to the form
my $coordinate_file = &MultiGetInputFile(1, $tmp_file_path."_coordinates.txt", 1);
$parameters .= " -i ".$coordinate_file;
&RSAT::message::Info("coordinate file", $coordinate_file) if ($main::echo >= 2);

## Check that a coordinate file has been given
unless ($coordinate_file) {
    &RSAT::error::FatalError("Input coordinates must be specified");
}

push @result_files, ("Coordinate file ($coordinate_format)",$coordinate_file);

## Check whether there is a list of common chr names for this org
my $common_names_file=$ENV{RSAT}."/data/genomes/$organism/genome/common_to_ensemble_chromosome_names.tab";
#print "$common_names_file";die;
if(-s $common_names_file) {
    $parameters .= " -common_chr $common_names_file ";
}
else{ $parameters .= " -check_chr "; }

## Compute genome fragment lengths from the coordinates
## this only works for BED, Bruno Mar2018
#my $length_file = $tmp_file_path."_lengths.tab";
#push @result_files, ("Genome fragment lengths",$length_file);
#my $seqlength_cmd = $SCRIPTS."/sequence-lengths -v 1 -i ".$coordinate_file;
#$seqlength_cmd .= " -in_format ".$coordinate_format;
#$seqlength_cmd .= " -o ".$length_file;
#system($seqlength_cmd);

## Repeats
if ($query->param('rm') =~ /on/) {
  $parameters .= " -rm ";
}


############################################################
## Output
$output_format = "fasta";

## Output file
$result_file = $tmp_file_path.".".$output_format;
$log_file = $tmp_file_path."_log.txt";
$parameters .= " -o ".$result_file;
&RSAT::message::Info("result_file", $result_file) if ($echo >= 0);
push @result_files, ("Result sequences (".$output_format.")", $result_file);
push @result_files, ("Command log (text)", $log_file); 

## Error log
$err_file = $tmp_file_path.".".$output_format."_err.txt";
$parameters .= " 2> ".$err_file;
push @result_files, ("Error log (text)",$err_file);

############################################################
## Report the command
&ReportWebCommand($command." ".$parameters);

################################################################
## Run the command
system($command." ".$parameters);

if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
  &PipingWarning();

  ## Print any errors to help user
  #if(-s $err_file) {
  #  print '<H4>Check the warnings/errors:</H4>';
  #
  #  print "<PRE>";
  #  open RESULT, $err_file;
  #  while (<RESULT>) {
  #      print $_ if(/^;INFO/ || /WARNING/); #unless ($query->param('output') =~ /server/i);
  #  }
  #  print "</PRE>";
  #  close RESULT;
  #}

  ## Print table with links to the result files
  &PrintURLTable(@result_files);

  ## Prepare data for piping
  &PipingFormForSequence($result_file, "fasta", "getfasta");

  print "<HR SIZE = 3>";

} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}

print $query->end_html;

exit(0);
