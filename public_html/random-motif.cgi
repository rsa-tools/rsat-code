#!/usr/bin/perl
############################################################
#
# random-motif
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
&RSA_header("random-motif result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$command = "$ENV{RSAT}/python-scripts/random-motif";
$convert_matrix_command = "$SCRIPTS/convert-matrix -from gibbs -return counts";
$prefix = "random-motif";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

################################################################
#
# read parameters
#
$parameters = '';

## Motif width
local $width = $query->param('width');
if (&RSAT::util::IsNatural($width)) {
  $parameters .= " -l ".$width;
} else {
  &RSAT::error::FatalError($width, "Is not a valid value for motif width. Should be a Natural number.");
}


## Conservation
local $conservation = $query->param('conservation');
if ((&RSAT::util::IsReal($conservation))
    && ($conservation >= 0) 
    && ($conservation <= 1)) {
  $parameters .= " -c ".$conservation;
} else {
  &RSAT::error::FatalError($conservation, "Is not a valid value for motif conservation. Should be a Real number between 0 and 1.");
}


## Multiplier
local $multiplier = $query->param('multiplier');
if ((&RSAT::util::IsNatural($multiplier)) 
    && ($multiplier >= 1)) {
  $parameters .= " -m ".$multiplier;
} else {
  &RSAT::error::FatalError($multiplier, "Is not a valid value for motif multiplier. Should be a Natural number >= 1.");
}

## Round numbers
if (lc($query->param("round")) eq "on") {
  $parameters .= " --round";
}


## Motif number
local $motif_nb = $query->param('motif_nb');
if (&RSAT::util::IsNatural($motif_nb)) {
  $parameters .= " -n ".$motif_nb;
} else {
  &RSAT::error::FatalError($motif_nb, "Is not a valid value for motif motif_nb. Should be a Natural number.");
}

## Concatenate parameters to the command
$command .= " ".$parameters;

local $tab_result_file = $tmp_file_path.".tab";
push @result_files, 'tab', $tab_result_file;

$command  .= " -o ".$tab_result_file;


## Convert the matrices
$output_format = $query->param('output_format');
$result_file = $tmp_file_path.".".$output_format;
if ($output_format ne "tab") {
  push @result_files, "$output_format", $result_file;
}
if ($output_format ne "tab") {
  $command .= "; ".$convert_matrix_command;
  $command  .= " -i ".$tab_result_file;
  $command .= " -prefix rand_";
  $command .= " -from tab -to ".$output_format;
  $command  .= " -o ".$result_file;
}

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

    ################################################################
    ## Table with links to the raw result files in various formats
    &PrintURLTable(@result_files);


    ## Form for sending results to other programs
    &PipingForm();

    &DelayedRemoval($tab_result_file);
    &DelayedRemoval($result_file);

    print "<hr size=\"3\">";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}

print $query->end_html;
exit(0);



### prepare data for piping
sub PipingForm {
  my $command = "$ENV{RSAT}/perl-scripts/convert-matrix -i $tab_result_file -from tab -to tab -top 1 -return counts";
  my $matrix_content = `$command`;
  $matrix_content =~ s|//\n||gm;
  $matrix_content =~ s|;.*\n||gm;
#  print "<pre>".$command."</pre>";
#  print "<pre>".$matrix_content."</pre>";

  $title = $query->param('title');
  $title =~ s/\"/\'/g;
    print <<End_of_form;
<hr size="3">
<table class="Nextstep">
<tr>
<td colspan="4">
<h3>Next step</h3>
</td>
</tr>

<tr>


<td valign="bottom" align="center">
<form method="POST" action="random-sites_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="submit" value="random sites">
</form>
</td>

<td valign="bottom" align="center">
<b><font color=red>new</a></b>
<form method="POST" action="matrix-scan_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="submit" value="pattern matching (matrix-scan)">
</form>
</td>

<td valign=bottom align=center>
<form method="POST" action="convert-matrix_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="hidden" name="logo" value="on" checked="checked">
<input type="submit" value="convert-matrix">
</form>
</td>

<td valign=bottom align=center>
<form method="post" target='_blank' action="http://meme.nbcr.net/meme4_3_0/cgi-bin/tomtom.cgi">
<input type="hidden" name="query" value="$matrix_content">
<input type="hidden" name="DIST" value="sandelin">
<input type="submit" value="TOMTOM">
</form>
Compare a single matrix to a motif database.
</td>
</tr>

</table>
End_of_form

}
