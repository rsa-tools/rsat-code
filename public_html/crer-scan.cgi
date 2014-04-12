#!/usr/bin/perl
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

### Print the header
&RSA_header("crer-scan result", "results");

## Check security issues
&CheckWebInput($query);

$command = $ENV{RSAT}."/python-scripts/crer_scan.py";
$prefix = "crer-scan";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
@result_files = ();

#### read parameters ####
$parameters = " -v 1";

################################################################
## Sites
if ($query->param('sites') =~ /\S/) {
  $query_file = $tmp_file_path."_query.txt";
  push @result_files, ("query",$query_file);
  open QUERY, ">".$query_file;
  print QUERY $query->param('sites');
  close QUERY;
  &DelayedRemoval($query_file);
  $parameters .= " -i ".$query_file;
} else {
    &cgiError("Please enter sites.\n");
}

################################################################
## site format
if ($query->param('in_format')) {
    my ($in_format) = split " ", $query->param('in_format'); ### take the first word
    $parameters .= " -in_format ".$in_format;
}


$result_file = $tmp_file_path.".tab";
push @result_files, ("gene info",$result_file);

&ReportWebCommand($command." ".$parameters);


## Update log file
&UpdateLogFile();

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


##############f##################################################
#
# Pipe the result to other commands
#
sub PipingForm {
    my $genes = `cat $result_file`;
    ### prepare data for piping
    print <<End_of_form;
<hr>
<table class='nextstep'>

<tr><td colspan=2>
<h3>Next step</h3>
</td></tr>

<tr>

<td><form method="post" action="retrieve-seq_form.cgi">
<input type="hidden" name="organism" VALUE="$organism">
<input type="hidden" name="genes" VALUE="selection">
<input type="hidden" name="gene_selection" VALUE="$genes">
<INPUT type="submit" value="retrieve sequences">
</form></td>

<td><form method="post" action="get-orthologs_form.cgi">
<input type="hidden" name="organism" VALUE="$organism">
<input type="hidden" name="sites" VALUE="$genes">
<INPUT type="submit" value="get orthologs">
</form></td>

</tr>


</TABLE>
End_of_form

}

