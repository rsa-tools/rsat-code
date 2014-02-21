#!/usr/bin/perl
############################################################
#
# random-sites
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
&RSA_header("random-sites result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();

$command = "$ENV{RSAT}/python-scripts/random-sites";
$convert_matrix_command = "$SCRIPTS/convert-matrix -return counts";
$prefix = "random-sites";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

################################################################
#
# Read matrix
#
$convert_matrix_parameters = '';

$matrix_file = $tmp_file_path.".input";
$tmp_matrix_file = $tmp_file_path.".matrix";

if ($query->param('matrix')) {
    open MAT, "> $matrix_file";
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
    $convert_matrix_parameters .= " -i $matrix_file";
} else {
    &RSAT::error::FatalError('You did not enter any data in the matrix box');
}

## Pseudo-counts
if (&IsReal($query->param('pseudo_counts'))) {
    $convert_matrix_parameters .= " -pseudo ".$query->param('pseudo_counts');
} else {
    &FatalError("Pseudo-count should be a real number");
}
if ($query->param('pseudo_distribution') eq "equi_pseudo") {
    $convert_matrix_parameters .= " -equi_pseudo ";
}

## input format
$matrix_input_format = lc($query->param('matrix_format'));
($matrix_input_format) = split (/\s+/, $matrix_input_format);
$convert_matrix_parameters .= " -from ".$matrix_input_format;

## output format
$convert_matrix_parameters .= " -to tab";
$convert_matrix_parameters .= " -o " .$tmp_matrix_file;

push (@result_files, "Input matrix ($matrix_input_format)", $matrix_file);
push (@result_files, "Converted matrix file (tab)", $tmp_matrix_file);

################################################################
#
# read parameters
#
$parameters = '';

## Number of sites
local $sites_nb = $query->param('sites_nb');
if (&RSAT::util::IsNatural($sites_nb)) {
  $parameters .= " -n ".$sites_nb;
} else {
  &RSAT::error::FatalError($sites_nb, "Is not a valid value for sites_nb. Should be a Natural number.");
}

## Concatenate parameters to the command
$command .= " ".$parameters;
$result_file = $tmp_file_path."_sites.fasta";
push (@result_files, "Sites (fasta)", $result_file);

#$result_file = $tab_result_file;
$command  .= " -o ".$result_file;
#$command  .= " -m ".$tmp_file_name_matrix;

## Convert the matrices
#if ($matrix_input_format ne "tab") {
  $command = $convert_matrix_command . $convert_matrix_parameters  . "; " . $command . " -m ". $tmp_matrix_file;
#} else {
  #$command = $command . " -m " . $matrix_file;
#}

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
    &PrintURLTable(@result_files);

    ## Form for sending results to other programs
    &PipingForm();

    #&DelayedRemoval($tab_result_file);
    &DelayedRemoval($result_file);

    print "<hr size=\"3\">";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

print $query->end_html;
exit(0);



### prepare data for piping
sub PipingForm {
  my $sites_data = `cat $result_file`;
#  $title = $query->param('title');
#  $title =~ s/\"/\'/g;
    print <<End_of_form;
<hr size="3">
<table class="Nextstep">
<tr>
<td colspan="3">
<h3>Next step</h3>
</td>
</tr>

<tr>

<td valign="bottom" align="center">
<form method="POST" action="implant-sites_form.cgi">
<input type="hidden" name="test" value="boum">
<input type="hidden" name="sites_file" value="$result_file">
<input type="hidden" name="sites" value="$sites_data">
<input type="submit" value="implant sites">
</form>
</td>
</tr>

</table>
End_of_form

}
