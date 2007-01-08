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
$default{matrix_format} = "tab";
$default{counts}="checked";
$default{consensus}="checked";
$default{frequencies}="";
$default{information}="";
$default{parameters}="";
$default{profile}="";
$default{weights}="";
$default{pseudo_weight}=1;
$default{margins}="checked";
$default{max_profile}=10;
$default{decimals}=4;


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
			 -Values=>['tab',
				   'Cluster-Buster',
				   'consensus (warning: requires the whole consensus file)', 
				   'gibbs',
				   'MEME',
				   'clustal',
				   'assembly',
				   ],
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
foreach my $stat qw(counts frequencies weights information margins consensus parameters profile) {
    print $query->checkbox(-name=>$stat,
			   -checked=>$default{$stat},
			   -label=>'');
    print "&nbsp;<A HREF='help.convert-matrix.html#",$stat,"'><B>", $stat, "</B></A>\n";
    print "<br>\n";
}


## Pseudo weight
print "<p><A HREF='help.convert-matrix.html#item_weight'><B>Pseudo-weight</B></A>&nbsp; ";
print $query->textfield(-name=>'pseudo_weight',
			-default=>$default{pseudo_weight},
			-size=>2);

## Decimals
print "&nbsp;"x5;
print "<A HREF='help.convert-matrix.html#decimals'><B>Decimals</B></A>&nbsp; ";
print $query->textfield(-name=>'decimals',
			-default=>$default{decimals},
			-size=>2);

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
$demo_matrix=`cat convert-matrix_demo_data.txt`;
print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'matrix_format',-default=>'tab');
print $query->hidden(-name=>'information',-default=>"on");
print $query->hidden(-name=>'weights',-default=>"on");
print $query->hidden(-name=>'parameters',-default=>"on");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.convert-matrix.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_PSSM.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

