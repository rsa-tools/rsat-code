#!/usr/bin/perl
#### this cgi script fills the HTML form for the program convert-matrix
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### default values for filling the form
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "consensus";
$default{alignment}="checked";
$default{consensus}="checked";
$default{frequencies}="";
$default{information}="";
$default{parameters}="";
$default{profile}="";
$default{weights}="";
$default{pseudo_weight}=1;
$default{margins}="checked";
$default{max_profile}=10;
$default{decimals}=2;


&ReadMatrixFromFile();

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 


################################################################
### print the form ###


################################################################
### header
&RSA_header("convert-matrix");
print "<CENTER>";
print "Convert different types of position-specific scoring matrices (PSSM), and calculate statistical parameters.<P>\n";
print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"convert-matrix.cgi");

#print "<FONT FACE='Helvetica'>";

################################################################
### input matrix
print "<B><A HREF='help.convert-matrix.html#matrix'>Input matrix</A></B>&nbsp;";

### Input matrix format 
print "&nbsp;"x10;
print "<B><A HREF='help.convert-matrix.html#matrix_format'>format</A></B>&nbsp;";
print $query->popup_menu(-name=>'matrix_format',
			 -Values=>['consensus', 
				   'gibbs',
				   'MEME',
				   'clustal'],
			 -default=>$default{matrix_format});
print "<BR>\n";

#### The matrix
print $query->textarea(-name=>'matrix',
		       -default=>$default{matrix},
		       -rows=>10,
		       -columns=>60);

################################################################
#### Return fields
print "<p><B><A HREF='help.convert-matrix.html#return'>Return fields</A></B>&nbsp;<br>\n";
my $i = 0;
foreach my $stat qw(alignment frequencies weights information margins consensus parameters profile) {
    print $query->checkbox(-name=>$stat,
			   -checked=>$default{$stat},
			   -label=>'');
    print "&nbsp;<A HREF='help.convert-matrix.html#",$stat,"'><B>", $stat, "</B></A>\n";
    print "<br>\n";
}



################################################################
### send results by email or display on the browser
print "<p>\n";
&SelectOutput("display");

################################################################
### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"convert-matrix_form.cgi");


#$demo_matrix = "A   \|   1   3   2   0   8   0   0   0   0   0   1   2
#C   \|   2   2   3   8   0   8   0   0   0   2   0   2
#G   \|   1   2   3   0   0   0   8   0   5   4   5   2
#T   \|   4   1   0   0   0   0   0   8   3   2   2   2";

$demo_matrix=`cat convert-matrix_demo_data.txt`;
print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.patser.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_patser.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

