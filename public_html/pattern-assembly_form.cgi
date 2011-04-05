#!/usr/bin/perl
#### this cgi script fills the HTML form for the program pattern-assembly
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{patterns} = "";
$default{pattern_file} = "";
$default{local_pattern_file} = "";
$default{maxfl} = 1;
$default{maxpat} = 200;
$default{subst} = 1;
$default{strand} = "insensitive";
$default{sc} = "none";

### print the form ###
&RSA_header("pattern-assembly", "form");
print "<CENTER>";
print "Assembly of patterns (oligos or dyads).<P>\n";
print "</CENTER>";
print "<HR>";
print "<blockquote>";

&ListParameters if ($ENV{rsat_echo} >=2);

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 


### if a pattern file is specified in the query,
### read patterns from this file
if (($pattern_file = $query->param("local_pattern_file")) &&
    (-e $pattern_file)) {
    open PAT, $pattern_file;
    while (<PAT>) {
	$default{patterns} .= $_;
    }
    close PAT;
}

print $query->start_multipart_form(-action=>"pattern-assembly.cgi");


### text area to enter the patterns
print "<A HREF='help.pattern-assembly.html#patterns'><B>\n";
print "Query pattern(s)</B></A><BR>\n";
print $query->textarea(-name=>'patterns',
		       -default=>$default{patterns},
		       -rows=>5,
		       -columns=>60);
print "<BR>\n";

#### upload patterns from a file on the client machine
print "<a href='help.pattern-assembly.html#pattern_file'>Upload pattern file</a><BR>";

print $query->filefield(-name=>'pattern_file',
			-default=>'starting value',
			-size=>30,
			-maxlength=>200);
print "<p>";
### maximum flanking size
print "<B><A HREF='help.pattern-assembly.html#maxfl'>Maximum flanking residues</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'maxfl',
			 -Values=>[0,1,2,3,4,5,6,7,8],
			 -default=>$default{maxfl});
print "<br>";

### maximum substitutions
print "<B><A HREF='help.pattern-assembly.html#subst'>Maximum substitutions</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'subst',
			 -Values=>[0,1,2,3],
			 -default=>$default{maxfl});
print "<br>";

### maximal number of patterns 
print "<B><A HREF='help.pattern-assembly.html#maxpat'>Maximum number of patterns</A>&nbsp;</B>\n";
print $query->textfield(-name=>'maxpat',
			-size=>4,
			-default=>$default{maxpat});
print "<br>";

### score column
print "<B><A HREF='help.pattern-assembly.html#sc'>Score column</A>&nbsp;</B>\n";
print $query->textfield(-name=>'sc',
			-size=>4,
			-default=>$default{sc});
print "<br>";

### strand ###
print "<B><A HREF='help.pattern-assembly.html#count_strands'>Strands</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['sensitive',
				  'insensitive'],
			 -default=>$default{strand});
print "<br>";


print "<HR width=550 align=left>\n";



### send results by email or display on the browser
&SelectOutput();


### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"pattern-assembly_form.cgi");
$demo_patterns = ";seq	score
acgtgc	2.82
tgccaa	2.52
ctgcac	1.30
acgtgg	1.10
cgcacg	1.00
cacgtg	0.94
aaacgt	0.48
aacgtg	0.11
cccacg	0.04
";
print "<TD><B>";
print $query->hidden(-name=>'patterns',-default=>$demo_patterns);
print $query->hidden(-name=>'maxfl',-default=>1);
print $query->hidden(-name=>'subst',-default=>1);
print $query->hidden(-name=>'sc',-default=>2);
print $query->hidden(-name=>'maxpat',-default=>100);
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.pattern-assembly.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.pattern-assembly.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_pattern-assembly.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);

