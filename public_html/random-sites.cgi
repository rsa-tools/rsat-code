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
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = "$ENV{RSAT}/python-scripts/random-sites";
$convert_matrix_command = "$SCRIPTS/convert-matrix -return counts";
$tmp_file_name = sprintf "random-sites.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("random-sites result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();

################################################################
#
# Read matrix
#
$convert_matrix_parameters = '';

$matrix_file = "$TMP/$tmp_file_name.input";
$tmp_matrix_file = "$TMP/$tmp_file_name.matrix";

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
$convert_matrix_parameters .= " -o " . $tmp_matrix_file;

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
$short_result_file = $tmp_file_name."_sites.fasta";
$result_file = $TMP."/".$short_result_file;
#$result_file = $tab_result_file;
$command  .= " -o ".$result_file;
#$command  .= " -m ".$tmp_file_name_matrix;

## Convert the matrices
#if ($matrix_input_format ne "tab") {
  $command = $convert_matrix_command . $convert_matrix_parameters  . "; " . $command . " -m ". $tmp_matrix_file;
#} else {
  #$command = $command . " -m " . $matrix_file;
#}

print "<pre>$command\n</pre>" if ($ENV{rsat_echo} >=1);

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
    &PrintURLTable(fasta=>$short_result_file);

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
