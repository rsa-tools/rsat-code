#!/usr/bin/perl
#### redirect error log to a file
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

$print_env_command = "$SCRIPTS/print-env";

### Print the header
$query = new CGI;

### print the result page
&RSA_header("print-env result");
&ListParameters if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

################################################################
## Execute the command
print "<PRE>\n";
$result = `$print_env_command`;
print $result;
print "</PRE>\n";

print $query->end_html;

exit(0);

