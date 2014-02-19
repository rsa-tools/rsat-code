#!/usr/bin/perl
############################################################
#
# implant-sites
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
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

### print the result page
&RSA_header("implant-sites result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$command = "$ENV{RSAT}/python-scripts/implant-sites";
#$convert_matrix_command = "$SCRIPTS/convert-matrix -return counts";
$prefix = "implant-sites";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "implant-sites.%s", &AlphaDate();

################################################################
#
# Read sequences
#
($sequence_file, $in_format) = &GetSequenceFile("fasta", no_format=>1, add_rc=>0);
push (@result_files, "Input sequence ($in_format)", $sequence_file);

################################################################
#
# Read sites
#
$sites_file = $tmp_file_path.".sites";
push (@result_files, "Input sites", $sites_file);
if ($query->param('sites')) {
    open SIT, "> $sites_file";
    print SIT $query->param('sites');
    close SIT;
    &DelayedRemoval($sites_file);
} else {
    &RSAT::error::FatalError('You did not enter any data in the sites box');
}

################################################################
#
# read parameters
#
$parameters = '';

## sequences
$parameters .= " -i " . $sequence_file;

## sites
$parameters .= " -s " . $sites_file;

# ## expected number of sites
# local $sites_nb = $query->param('sites_nb');
# if (&RSAT::util::IsReal($query->param('sites_espp'))) {
#   $parameters .= " --espp=". $query->param('sites_espp');
# } else {
#   &RSAT::error::FatalError($query->param('sites_espp'), "Is not a valid value for sites_espp. Should be a Real number.");
# }

## Concatenate parameters to the command
$command .= " ".$parameters;
$result_file = $TMP."/".$tmp_file_name.".fasta";
push @result_files, "Result (fasta)", $result_file;
$out_format = "fasta"; ## implant-site always exports fasta, but this variable is required for the piping form
$command  .= " -o ".$result_file;

&ReportWebCommand($command);

if ($query->param('output') eq "display") {
    &PipingWarning();

    ### execute the command ###
    system "$command";

    ### Print result on the web page
    print '<h4>Result</h4>';
    print "<pre>";
    print `cat $result_file`;
    print "</pre>";


    ## Add link to the result file
    $result_file = $tmp_file_name.".". "fasta";
    &PrintURLTable(@result_files);

    ## Form for sending results to other programs
    &PipingFormForSequence();

    #&DelayedRemoval($tab_result_file);
    &DelayedRemoval($result_file);

    print "<hr size=\"3\">";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

print $query->end_html;
exit(0);


