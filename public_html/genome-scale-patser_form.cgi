#!/usr/bin/perl
#### this cgi script fills the HTML form for the program genome-scale-patser
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "$RSA/public_html/genome-scale.lib.pl";

### Read the CGI query
$query = new CGI;

### default values for filling the form
$default{matrix_format} = "consensus";
$default{matrix} = "";
$default{strands} = "both";
$default{return} = "all matching positions";
$default{lthreshold} = "0";
$default{uthreshold} = "none";
$default{alphabet} = "a:t 0.3 c:g 0.2";
$default{pseudo_counts} = 1;


$default{sequence_format} = "wc"; ### Important ! the format supported by patser

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### if a matrix file is specified in the query,
### read matrix from this file
if (($matrix_file = $query->param("matrix_file")) &&
     (-e $matrix_file)) {
  open PAT, $matrix_file;
  while (<PAT>) {
    $default{matrix} .= $_;
  }
  close PAT;
}

### print the form ###
&RSA_header("genome-scale patser");

### head
print "<CENTER>";
print "Scan all upstream or downstream regions with a profile matrix<BR>\n";
print "Program developed by <A HREF='mailto:hertz\@colorado.edu (Jerry Hertz)'>Jerry Hertz</A><BR>";
print "Web interface by <A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>Jacques van Helden</A><P>";
print "</CENTER>";

#&ListParameters;


print $query->start_multipart_form(-action=>"genome-scale-patser.cgi");

################################################################
#
# retrieve-seq options
#

&DisplayRetrieveSeqOptions();

################################################################
#
# patser options
#
### text area to enter the matrix
print "<A HREF='help.patser.html#matrix'><B>\n";
print "Matrix</B></A>\n";
print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
print "<B>Format</B>&nbsp;";
print $query->popup_menu(-name=>'matrix_format',
			 -Values=>['consensus',
				   'gibbs',
				   'transfac'
				   ],
			 -default=>$matrix_format);
print "<BR>\n";
print $query->textarea(-name=>'matrix',
		       -default=>$default{matrix},
		       -rows=>4,
		       -columns=>60);
print "<BR>\n";

### strands
print "<A HREF='help.patser.html#strands'><B>Search strands</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'strands',
			 -Values=>['single',
				   'both'],
			 -default=>$default{strands});

### return
print "<A HREF='help.patser.html#return'><B>Return</B></A>&nbsp;\n";
print $query->popup_menu(-name=>'return',
			 -Values=>['top value for each sequence',
				   'all matching positions'],
			 -default=>$default{return});

### thresholds
print "<BR>\n";
print CGI::table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       CGI::Tr({-align=>left,-valign=>MIDDLE},
		       [
		      CGI::td({-align=>left,-valign=>MIDDLE},
			      [
			       "<B><A HREF='help.patser.html#pseudo_counts'>Pseudo-counts</A>\n",
			       $query->textfield(-name=>'pseudo_counts',
						       -default=>$default{pseudo_counts},
						       -size=>2),
			       "&nbsp;&nbsp;<b>Thresholds</b>",
			       "<A HREF='help.patser.html#lthreshold'><B> lower</B></A>",
			       $query->textfield(-name=>'lthreshold',
						 -default=>$default{lthreshold},
						 -size=>2),
			       "<A HREF='help.patser.html#uthreshold'><B> upper</B></A>",
			       $query->textfield(-name=>'uthreshold',
						 -default=>$default{uthreshold},
						 -size=>2)
			       
			       ]),
			])
		 );
print "<BR>\n";

### alphabet
print "<B><A HREF='help.patser.html#alphabet'>\n";
print "Alphabet</A>\n";
print $query->textfield(-name=>'alphabet',
			-default=>$default{alphabet},
			-size=>50);
print "<BR>\n";

### send results by e-mail or display on the browser
&SelectOutput;

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"genome-scale-patser_form.cgi");

$demo_matrix = "A    |    1    3    2    0    8    0    0    0    0    0    1    2
C    |    2    2    3    8    0    8    0    0    0    2    0    2
G    |    1    2    3    0    0    0    8    0    5    4    5    2
T    |    4    1    0    0    0    0    0    8    3    2    2    2";

print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'sequence',-default=>$demo_sequence);
print $query->hidden(-name=>'lthreshold',-default=>9);
print $query->hidden(-name=>'alphabet',-default=>"a:t 0.325 c:g 0.175");
print $query->hidden(-name=>'sequence_format',-default=>$default{sequence_format});
print $query->hidden(-name=>'set_name',-default=>'upstream sequences from the yeast PHO genes');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.patser.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_genome-scale-patser.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);





