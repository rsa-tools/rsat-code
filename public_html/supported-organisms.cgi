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

### print the header
&RSA_header("Supported organisms", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);


##$tmp_file_name = sprintf "supported-organisms.%s", &AlphaDate();
$prefix = "supported-organisms";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

$font{variable} = 1;
$command = "$SCRIPTS/supported-organisms -v 1";

$parameters = " -return ID,nb,source,last_update,taxonomy,up_from,up_to";

################################################################
## treat taxon specificity of the server if required
if ($ENV{GROUP_SPECIFICITY}) {
  $command .= " -group ".$ENV{GROUP_SPECIFICITY};
}
$command .= " ".$parameters;
&ReportWebCommand($command) if ($ENV{rsat_echo} >= 1);

print "<pre>", `$command`, "</pre>";

### execute the command ###
$command .= " | awk '{print \$0\"\t<a href=$ENV{rsat_www}/data/genomes/\"\$1\"/>data</a>\"}'";
$command .= " | perl -pe 's|<a href=.+/data/genomes/;/>data</a>||'";

open RESULT, "$command|";

### Print result on the web page
print "<CENTER>\n";
&PrintHtmlTable(RESULT, $tmp_file_path, 0, 10000);
print "</CENTER>\n";

close(RESULT);

print '<HR SIZE=3>';


print $query->end_html;

exit(0);

