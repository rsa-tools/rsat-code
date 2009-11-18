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
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = "$ENV{RSAT}/python-scripts/implant-sites";
#$convert_matrix_command = "$SCRIPTS/convert-matrix -return counts";
$tmp_file_name = sprintf "implant-sites.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("implant-sites result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();


################################################################
#
# Read sequences
#
($sequence_file, $sequence_format) = &GetSequenceFile("fasta", no_format=>1, add_rc=>0);

################################################################
#
# Read sites
#
$sites_file = "$TMP/$tmp_file_name.sites";

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

## expected number of sites
local $sites_nb = $query->param('sites_nb');
if (&RSAT::util::IsReal($query->param('sites_espp'))) {
  $parameters .= " --espp=". $query->param('sites_espp');
} else {
  &RSAT::error::FatalError($query->param('sites_espp'), "Is not a valid value for sites_espp. Should be a Real number.");
}

## Concatenate parameters to the command
$command .= " ".$parameters;
$result_file = $TMP."/".$tmp_file_name.".fasta";
$command  .= " -o ".$result_file;

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

    ################################################################
    ## Table with links to the raw result files in various formats
    #$tab_result_URL = $ENV{rsat_www}."/tmp/".$tmp_file_name.".tab";
    print "<center><table class=\"nextstep\">\n";
    print "<tr><td colspan='3'><h3>Result file(s)</h3> </td></tr>";
    # print ("<tr>",
    #      "<td>tab</td>",
    #      "<td>","<a href='".$tab_result_URL."'>".$tab_result_URL."</a>","</td>",
    #      "</tr>");
      $result_URL = $ENV{rsat_www}."/tmp/".$tmp_file_name.".". "fasta";
      print ("<tr>",
	     "<td>".$output_format."</td>",
	     "<td>","<a href='".$result_URL."'>".$result_URL."</a>","</td>",
	     "</tr>");
    print "</table></center>";

    ## Form for sending results to other programs
    &PipingForm();

    #&DelayedRemoval($tab_result_file);
    &DelayedRemoval($result_file);

    print "<hr size=\"3\">";
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

print $query->end_html;
exit(0);



### prepare data for piping
sub PipingForm {
  $title = $query->param('title');
  $title =~ s/\"/\'/g;
    print <<End_of_form;
<hr size="3">
<table class="Nextstep">
<tr>
<td colspan="3">
<h3>Next step</h3>
</td>
</tr>

<tr>

<!--
<td valign="bottom" align="center">
<b><font color=red>new</a></b>
<form method="POST" action="matrix-scan_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="submit" value="pattern matching (matrix-scan)">
</form>
</td>

-->
</tr>
</table>
End_of_form

}
