#!/usr/bin/perl

################################################################
#### this cgi script fills the HTML form for the program roc-stats
################################################################

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
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
## Output fields
my @output_fields = qw(
		       occ
		       TP
		       FP
		       FN
		       Sn
		       PPV
		       FPR
		       Acc_g
		       Acc_a
		       graphs
		      ); # add AUC when possible
my %field_description = ();
$field_description{occ} = "Inverse cumulative occurrences";
$field_description{TP} = "TP (True Positives)";
$field_description{FP} = "FP (False Positives)";
$field_description{FN} = "FN (False Negatives)";
$field_description{Sn} = "Sn (Sensitivity)";
$field_description{PPV} = "PPV (Positive Predictive Value)";
$field_description{FPR} = "FPR (False Positive Rate)";
$field_description{Acc_g} = "Geometic Accuracy";
$field_description{Acc_a} = "Arithmetic Accuracy";
$field_description{AUC} = "AUC (Area Under the Curve)";
$field_description{graphs} = "Graphs (distribution of statistics, ROC and PR curves)";

################################################################
## Output graphs
my @output_graphs = qw(
		       ROC
		       PR
		       Stats
		       Stats_xlog
		      );
$field_description{ROC} = "ROC curves";
$field_description{PR} = "Precision-Recall";
$field_description{Stats} = "Distribution of evaluation statistics";
$field_description{Stats_xlog} = "Distribution of evaluation statistics (logarithm scale on X axis)";

################################################################
### default values for filling the form
$default{occ} = "checked";
$default{TP} = "checked";
$default{FP} = "checked";
$default{FN} = "checked";
$default{Sn} = "checked";
$default{PPV} = "checked";
$default{FPR} = "checked";
$default{Acc_g} = "checked";
$default{Acc_a} = "checked";
$default{AUC} = "";

$default{scores} = $query->param("scores");
$default{scores} =~ s/\"//g; #### remove quotes for security reasons (avoid imbedded command)
$default{scores} =~ s/\r//g; #### remove quotes for security reasons (avoid imbedded command)

$default{sc_col} = "1";
$default{status_col} = "2";
$default{pos} = 'pos';
$default{neg} = 'neg';
$default{total} = "";
$default{graphs}="checked";
$default{img_format}="png";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}



################################################################
### Print the form
################################################################

################################################################
### Header
# print "<head>
#    <title>NeAT - roc-stats</title>
#    <link rel=\"stylesheet\" type=\"text/css\" href = \"main_grat.css\" media=\"screen\">
# </head>";

&RSA_header("roc-stats", "form");
print "<CENTER>";
print "This program takes as input a set of scored results associated with validation labels (pos for positive, neg for negative) and computes, for each score value, the derived statistics (Sn, PPV, FPR), which can be further used to draw a ROC curve.<P>\n";
print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "<br>This program was developed by <a target=_blank href='http://www.bigre.ulb.ac.be/people/Members/rekins'>Rekins Janky</a> and <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.</center>";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"roc-stats.cgi");

print "<hr>";

################################################################
### Input scores
print "<B><A HREF='help.roc-stats.html#scores'>Input scores</A></B><br>";
print $query->textarea(-name=>'scores',
		       -default=>$default{scores},
		       -rows=>6,
		       -columns=>65);

### option to upload a file with the scores from the client machine
print "<BR>Upload scores from file<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);

################################################################
#### Input parameters
print "<HR><b><a href = 'help.roc-stats.html#params'>Input parameters</a></b><br>";
print "<table><tr><td><B><a href = 'help.roc-stats.html#scores'>Scores column</a></B></td><td><input type = 'text' name='sc_col' value = '".$default{sc_col}."' size = 1></input></td>";#</tr>";
print "<td><B><a href = 'help.roc-stats.html#pos'>Positive labels</a></B></td><td><input type = 'text' name='pos' value = '".$default{pos}."' size = 1></input></td>";
print "<td><B><a href = 'help.roc-stats.html#total'>Total Number of elements</a></B></td><td><input type = 'text' name='total' value = '".$default{total}."' size = 1></input></td></tr>";
print "<tr><td><B><a href = 'help.roc-stats.html#status'>Status column</a></B></td><td><input type = 'text' name='status_col' value = '".$default{status_col}."' size = 1></input></td>";
print "<td><B><a href = 'help.roc-stats.html#neg'>Negative labels</a></B></td><td><input type = 'text' name='neg' value = '".$default{neg}."' size = 1></input></td></tr></table>";

################################################################
#### Return fields
print "<HR><BR>\n";
print "<B><A HREF='help.roc-stats.html#return'>Return fields</A></B>&nbsp;<br><UL>\n";
my $i = 0;
foreach my $field (@output_fields) {
  $i++;
  print $query->checkbox(-name=>$field,
			 -checked=>$default{$field},
			 -label=>'');
  print "&nbsp;<A HREF='help.roc-stats.html#",$field,"'><B>", $field_description{$field}, "</B></A>\n";
  if ($field eq "graphs"){
    print "&nbsp;&nbsp;<A HREF='help.roc-stats.html#img'>Image format (if graphs is checked)</A>&nbsp;\n";
    print $query->popup_menu(-name=>'img_format',
			     -Values=>['png',
				       'jpg',
				       'gif',
				       'eps',
				       'pdf'],
			     -default=>$default{img_format});
  }
print "<BR>\n";
}
print "</ul><p>\n";

################################################################
#### Return graphs
#print "<p><B><A HREF='help.roc-stats.html#graphs'>Graphs parameters (if graphs output is checked)</A></B>&nbsp;<br>\n";

################################################################
### send results by email or display on the browser
&SelectOutput("display");

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"roc-stats_form.cgi");
my $demo_scores="0.95	pos
0.85	neg
0.75	pos
0.65	pos
0.55	pos
0.45	neg
0.35	pos
0.35	neg
0.25	neg
0.15	pos
0.05	neg
";
#print $demo_scores;
print "<TD><B>";
print $query->hidden(-name=>'scores',-default=>$demo_scores);
print $query->hidden(-name=>'sc_col',-default=>'1');
print $query->hidden(-name=>'status_col',-default=>'2');
print $query->hidden(-name=>'pos',-default=>'pos');
print $query->hidden(-name=>'neg',-default=>'neg');
print $query->hidden(-name=>'AUC',-default=>"checked");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.roc-stats.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_roc.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

