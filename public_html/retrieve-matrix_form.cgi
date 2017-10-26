#!/usr/bin/perl
#### this cgi script fills the HTML form for the program retrieve-matrix
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
### default values for filling the form
$default{output}="display";
$default{input}="";
$default{id}="";
$default{id_file} = "tab";
$default{table} = 1;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
    if ($query->param($key)) {
        $default{$key} = $query->param($key);
    }
}

################################################################
### header
&RSA_header("retrieve-matrix", "form");
print "<style>
	table.result {
		border-collapse: collapse;		
	}
	table.result th, table.result td {
		border: 1px solid #cbcbb4;
		padding: 15px;
	}
	table.resultlink td, table.resultlink th{
		font-size: 100%;
	}
</style>";

print "<CENTER>";
print "Retrieve the matrices with identifiers.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

################################################################
#### collections
print "<hr>";
&MotifSelection("retrieve" => 1);

print "<div id='result' style='display:none'></div>";

print "<div id='outputurl'></div>";

################################################################
### action buttons
#print "<UL><UL><TABLE class='formbutton'>\n";
#print "<TR VALIGN=MIDDLE>\n";
#print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
#print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
#print $query->end_form;

################################################################
### data for the demo

#print "<TD><B>";

print $query->end_html;

exit(0);

