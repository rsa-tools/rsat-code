#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
$print_env_command = "$SCRIPTS/print-env";
$tmp_file_name = sprintf "print-env.%s", &AlphaDate;

### Print the header
$query = new CGI;
print $query->header;
print $query->start_html;


##### update log file ####
#&UpdateLogFile;

#### execute the command #####

print "<PRE>\n";
$result = `$print_env_command`;
print $result;
print "Hello\n";
print "</PRE>\n";

print $query->end_html;

exit(0);

