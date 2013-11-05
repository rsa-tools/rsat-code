#!/usr/bin/perl
#### this cgi script fills the HTML form for the program classfreq
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

################################################################
### default values for filling the form
$default{ci} = 'auto';
$default{col} = '1';
$default{min} = 'auto';
$default{max} = 'auto';
$default{from} = 'auto';
$default{to} = 'auto';


## Replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

################################################################
## Print the form


################################################################
## Header
&RSA_header("Frequency distribution (<i>classfreq</i>)");
print "<center>";
print "<p>Compute frequency distribution of numerical values provided in a given column of a tab-delimited text file.</p>\n";
print "</center>";
print "<blockquote>\n";

print $query->start_multipart_form(-action=>"classfreq.cgi");

print "<font face='Helvetica'>";


################################################################
## Input data table


#### data from pipe (compare-graphs)
if ($query->param('transferred_file')) {
  my $transferred_file = $query->param('transferred_file');
  my $file_url = $transferred_file;
  $file_url =~ s|$ENV{RSAT}/public_html|$ENV{rsat_www}|;
  $file =~ s|$ENV{rsat_www}|$ENV{RSAT}/public_html|;
  print "<ul><a href=$file_url>";
  print " transferred from previous query<BR>\n";
  print "</a></ul>";
  print "<input type='hidden' NAME='transferred_file' VALUE='$transferred_file'>\n";

} else {
  $default{data} = $query->param('data');
  $default{data} =~ s/\"//g; #### remove quotes for security reasons (avoid imbedded command)
  $default{data} =~ s/\r//g; #### remove quotes for security reasons (avoid imbedded command)
  print $query->textarea(-name=>'data',
			 -default=>$default{data},
			 -rows=>10,
			 -columns=>65);

  ### option to upload a file with the data from the client machine
  print "<BR>Upload data from file<BR>\n";
  print $query->filefield(-name=>'uploaded_file',
			  -default=>$default{uploaded_file},
			  -size=>45,
#			  -maxlength=>200,
			 );
}

# ## Fill data in text area
# print "<B><A HREF='help.classfreq.html#data'>Input data table</A></B>&nbsp;";
# print "<BR>\n";
# print $query->textarea(-name=>'data',
# 		       -default=>$default{data},
# 		       -rows=>10,
# 		       -columns=>60);
#
# ### option to upload a file
# print "<BR>Upload data file<BR>\n";
# print $query->filefield(-name=>'uploaded_file',
# 			-default=>'',
# 			-size=>45,
# 			-maxlength=>200);
# print "<BR>\n";

################################################################
## Parameters

print "<h2>Parameters</h2>\n";


print "Class interval",
  $query->textfield(-name=>'ci',
		    -default=>$default{ci},
		    -size=>5);
print "Data column",
  $query->textfield(-name=>'col',
		    -default=>$default{col},
		    -size=>5);

## Values to take into account
print "<br><b>Values to take into account</b>";
print "&nbsp"x5, "min", 
  $query->textfield(-name=>'min',
		    -default=>$default{min},
		    -size=>5);
print "&nbsp"x5, "max", 
  $query->textfield(-name=>'max',
		    -default=>$default{max},
		    -size=>5);

## Values to report in the class frequency table
print "<br><b>Values to report in the class frequency table</b>";
print "&nbsp"x5, "from", 
  $query->textfield(-name=>'from',
		    -default=>$default{from},
		    -size=>5);
print "&nbsp"x5, "to", 
  $query->textfield(-name=>'to',
		    -default=>$default{to},
		    -size=>5);


$default{from} = 'auto';
$default{to} = 'auto';


################################################################
### send results by email or display on the browser
print "<hr>";
&SelectOutput();

################################################################
### action buttons
print "<ul><ul><table>\n";
print "<tr valign='middle'>\n";
print "<td>", $query->submit(-label=>"GO"), "</td>\n";
print "<td>", $query->reset, "</td>\n";
print $query->end_form;

################################################################
## Data for the demo
print $query->start_multipart_form(-action=>"classfreq_form.cgi");
$demo_data = `cat demo_files/allup500_Saccharomyces_cerevisiae_some_pattern_counts.tab`;
print "<TD><B>";
print $query->hidden(-name=>'data',-default=>$demo_data);
print $query->hidden(-name=>'ci',-default=>'1');
print $query->hidden(-name=>'col',-default=>'4');
print $query->submit(-label=>"DEMO");
print "</b></td>\n";
print $query->end_form;


print "<TD><B><A HREF='help.classfreq.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_classfreq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE class='formbutton'></UL></UL>\n";

print "</BLOCKQUOTE>\n";
print "</FONT>\n";

print $query->end_html;

exit(0);

