#!/usr/bin/perl
#### this cgi script fills the HTML form for the program random-genes
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{organism} = "Saccharomyces cerevisiae";
$default{gene_nb} = 25;
$default{replacement} = "";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("random gene selection");

### head
print "<CENTER>";
print "Select random genes for a given organism.<P>\n";
print "</CENTER>";


print $query->start_multipart_form(-action=>"random-genes.cgi");

print "<FONT FACE='Helvetica'>";


#### number of genes
print " <A HREF='help.random-genes.html#gene_nb'>Number of genes</A> ";
print $query->textfield(-name=>'gene_nb',
			-default=>$default{gene_nb},
			-size=>5);
#### organism
print "&nbsp;&nbsp;&nbsp;\n";
print $query->checkbox(-name=>'replacement',
		       -checked=>$default{'replacement'},
		       -label=>''
		       );
print "&nbsp;<a href=help.random-genes.html#replacement>With replacement</a>";



 


#### organism
print "<P>\n";
&OrganismPopUp();



### send results by email or display on the browser
print "<P>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.random-genes.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_random-genes.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

