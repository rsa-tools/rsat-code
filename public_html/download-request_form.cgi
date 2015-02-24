#!/usr/bin/perl
#### this cgi script fills the HTML form for the program convert-matrix
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;


################################################################
## Print the form


################################################################
## Header
&RSA_header("Download request", "form");
print "<center>";
print "Request a stand-alone version of the RSAT/NeAT software suites.<P>\n";
print "</center>";
print "<blockquote>\n";

print $query->start_multipart_form(-action=>"download-request.cgi");

################################################################
## User-specified parameters

print '<table border=0>';

print "<tr><td>\n";
print "<b>First name</b> <font color='red'>*</font>\n";
print "</td>\n<td>";
print $query->textfield(-name=>'first_name', -size=>50);
print "</td></tr>\n";

print "<tr><td>\n";
print "<br><b>Last name</b> <font color='red'>*</font>\n";
print "</td>\n<td>";
print $query->textfield(-name=>'last_name', -size=>50);
print "</td></tr>\n";

print "<tr><td>\n";
print "<br><b>Email</b> <font color='red'>*</font>\n";
print "</td>\n<td>";
print $query->textfield(-name=>'email_address', -size=>50);
print "</td></tr>\n";

print "<tr><td>\n";
print "<br><b>Institution</b> <font color='red'>*</font>\n";
print "</td>\n<td>";
print $query->textfield(-name=>'institution', -size=>50);
print "</td></tr>\n";

print "<tr><td>\n";
print "<br><b>Country</b> <font color='red'>*</font>\n";
print "</td>\n<td>";
print $query->textfield(-name=>'country', -size=>50);
print "</td></tr>\n";

print "</table>";


################################################################
## Action buttons
print "<ul><ul><table class='formbutton'>\n";
print "<tr valign=middle>\n";
print "<td>", $query->submit(-label=>"GO"), "</td>\n";
print "<td>", $query->reset, "</td>\n";
print $query->end_form;

################################################################
### data for the demo 

print $query->end_html;

exit(0);

