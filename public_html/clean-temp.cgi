#!/usr/bin/perl
#### redirect error log to a file
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = $ENV{RSAT}."/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log to $ERR_LOG";
    carpout(*LOG);
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
use RSAT::util;
use RSAT::server;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";


### Print the header
$query = new CGI;

### print the result page
&RSA_header("clean-temp result");
&ListParameters if ($ENV{rsat_echo} >=2);


#
## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

#### execute the command #####

## Parameters
my $clean_limit = 3;
$clean_temp_command = "echo '<br>started'; date; \n";
$clean_temp_command .= "find ".$ENV{RSAT}."/public_html/tmp/ -mtime +".${clean_limit}." -type f -exec rm -f {} \\; ; \n";
$clean_temp_command .= "find ".$ENV{RSAT}."/public_html/tmp/ -mtime +".${clean_limit}." -type d -exec rm -rf {} \\; ; \n";
$clean_temp_command .= "rm -f ".$ENV{RSAT}."/public_html/tmp/serialized_genomes/*.serial ; \n";
$clean_temp_command .= "echo '<br>done\n'; date; \n";

print "<h2>Disk free before cleaning</h2>";
print "<pre>"; system("df -h $ENV{RSAT}"); print "</pre>";

&RSAT::message::TimeWarn("Cleaning temporary directory from files older than ".$clean_limit." days + all serialized files");

print "<PRE>\n";
print $clean_temp_command if ($ENV{rsat_echo} >= 2);
print "</PRE>\n";

$err = system($clean_temp_command);

print "<hr><h2>Disk free after cleaning</h2>";
print "<pre>"; system("df -h $ENV{RSAT}"); print "</pre>";

&RSAT::message::TimeWarn("Cleaning finished");


print $query->end_html;

exit(0);

