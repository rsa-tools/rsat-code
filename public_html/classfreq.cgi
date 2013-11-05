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

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("classfreq result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/classfreq";
$prefix = "classfreq";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();

#### read parameters ####
$parameters = " -v 1";

################################################################
## Input data table

## Define data file
$data_file = $tmp_file_path."_input.tab";


if ($query->param('transferred_file') =~ /\S/) {
  ## Use local data file (transferred from previous query)
  $data_file = $query->param('transferred_file');
} elsif ($query->param('uploaded_file')) {
  ## Upload file from the client computer
  open DATA, ">$data_file";
  $upload_data_file = $query->param('uploaded_file');
  $type = $query->uploadInfo($upload_data_file)->{'Content-Type'};
  while (<$upload_data_file>) {
    print DATA;
  }
  close DATA;
  &DelayedRemoval($data_file);

} elsif ( $query->param('data') =~ /\S/) {
  ## Fill file with data from the text area
  open DATA, ">".$data_file;
  print DATA $query->param('data');
  close DATA;
  &DelayedRemoval($data_file);

} else {
    &cgiError("The data box cannot be empty\n");
}

$parameters .= " -i ".$data_file;
push @result_files, ("data",$data_file);

## Parameters

## Class interval
if (&IsReal($query->param('ci'))) {
  $parameters .= " -ci ".$query->param('ci');
}

## Data column
if (&IsReal($query->param('col'))) {
  $parameters .= " -col ".$query->param('col');
}

## Min
if (&IsReal($query->param('min'))) {
  $parameters .= " -min ".$query->param('min');
}

## Max
if (&IsReal($query->param('max'))) {
  $parameters .= " -max ".$query->param('max');
}

## From
if (&IsReal($query->param('from'))) {
  $parameters .= " -from ".$query->param('from');
}

## To
if (&IsReal($query->param('to'))) {
  $parameters .= " -to ".$query->param('to');
}

$result_file = $tmp_file_path."_classfreq.tab";
push @result_files, ("class frequencies",$result_file);

&ReportWebCommand($command." ".$parameters);

################################################################
## Run the command
if ($query->param('output') eq "display") {
    &PipingWarning();

    open RESULT, "$command $parameters |";

    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);

    &PrintURLTable(@result_files);
    &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}
print $query->end_html();

exit(0);


################################################################
#
# Pipe the result to other commands
#
sub PipingForm {
    my $data = `cat $result_file`;
    ### prepare data for piping
    print <<End_of_form;
<hr>
<table class='nextstep'>

<tr><td colspan=2>
<h3>Next step</h3>
</td></tr>

<tr>

<td><form method="post" action="XYgraph_form.cgi">
<input type="hidden" name="data" VALUE="$data">
<input type="hidden" name="xcol" VALUE="3">
<input type="hidden" name="ycol" VALUE="4,5,6">
<input type="hidden" name="title" VALUE="Frequency distribution">
<input type="hidden" name="xsize" VALUE="600">
<input type="hidden" name="xleg1" VALUE="Values">
<input type="hidden" name="ysize" VALUE="400">
<input type="hidden" name="yleg1" VALUE="Frequencies">
<input type="hidden" name="lines" VALUE="on">
<input type="submit" value="XYgraph">
</form></td>

</tr>


</TABLE>
End_of_form

}

