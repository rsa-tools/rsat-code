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
my @output_fields = qw(stats
		       graphs
		      ); # add AUC when possible
#my @output_fields = qw(
#		       occ
#		       TP
#		       FP
#		       FN
#		       Sn
#		       PPV
#		       FPR
#		       Acc_g
#		       Acc_a
#		       graphs
#		      ); # add AUC when possible
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
$field_description{stats} = "Table of statistics (Inverse cumulative occurrences, Sn, PPV, Accuracy)";
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
$default{stats} = "checked";
$default{AUC} = "";
$default{uploaded_file}="";

$default{demo_comment} = $query->param('demo_comment');
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
&NeAT_header("roc-stats", "form");
print "<CENTER>";
print "This program takes as input a set of scored results associated with validation labels (pos for positive, neg for negative) and computes, for each score value, the derived statistics (Sn, PPV, FPR), which can be further used to draw a ROC curve.<P>\n";
#print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "<br>This program was developed by <a target=_blank href='http://www.bigre.ulb.ac.be/people/Members/rekins'>Rekins Janky</a> and <a target=_blank href=http://www.bigre.ulb.ac.be/Users/jvanheld/>Jacques van Helden</a>.</center>";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"roc-stats.cgi");

print "<hr>";

if ($query->param('demo_comment')){
  print "<H4>Comment on the demonstration example : </H4><blockquote class ='demo'>In this demonstration, we use the results of the comparison of the String database network on <i>E.coli</i> genes and the gene-gene network inferred from RegulonDB database (performed in November 2007).</blockquote><br>\n";
}

################################################################
### Input data

print "<B><A HREF='help.roc-stats.html#input_format'>Input data</A></B><br>";

#### data from pipe (compare-graphs)
if ($query->param('roc-stats_graph_file')) {
  my $file_url = $query->param('roc-stats_graph_file');
  $file =~ s|$ENV{rsat_www}|$ENV{RSAT}/public_html|;
  print "<ul><a href=$file_url>";
  print " transferred from previous query<BR>\n";
  print "</a></ul>";
#  $default{uploaded_file} = $file;
} else {
  $default{data} = $query->param('data');
  $default{data} =~ s/\"//g; #### remove quotes for security reasons (avoid imbedded command)
  $default{data} =~ s/\r//g; #### remove quotes for security reasons (avoid imbedded command)
  print $query->textarea(-name=>'data',
			 -default=>$default{data},
			 -rows=>6,
			 -columns=>65);

  ### option to upload a file with the data from the client machine
  print "<BR>Upload data from file<BR>\n";
  print $query->filefield(-name=>'uploaded_file',
			  -default=>$default{uploaded_file},
			  -size=>45,
			  -maxlength=>200);
}

################################################################
#### Input parameters
print "<HR><b><a href = 'help.roc-stats.html#input_format'>Input parameters</a></b><UL>";
print "<table><tr><td><B><a href = 'help.roc-stats.html#scores'>Scores column</a></B></td><td><input type = 'text' name='sc_col' value = '".$default{sc_col}."' size=1></input></td>";#</tr>";
print "<td><B><a href = 'help.roc-stats.html#pos'>Positive labels</a></B></td><td>";
print $query->textarea(-name=>'pos',
		       -default=>$default{pos},
		       -rows=>1,
		       -columns=>5);
#"<input type = 'text' name='pos' value = '".$default{pos}."' size=1></input>
print "</td>";
print "<td><B><a href = 'help.roc-stats.html#total'>Total Number of elements</a></B></td><td><input type = 'text' name='total' value = '".$default{total}."' size=1></input></td></tr>";
print "<tr><td><B><a href = 'help.roc-stats.html#status'>Status column</a></B></td><td><input type = 'text' name='status_col' value = '".$default{status_col}."' size=1></input></td>";
print "<td><B><a href = 'help.roc-stats.html#neg'>Negative labels</a></B></td><td>";
print $query->textarea(-name=>'neg',
		       -default=>$default{neg},
		       -rows=>1,
		       -columns=>5);
#"<input type = 'text' name='neg' value = '".$default{neg}."' size=1></input>
print "</td></tr></table>";
print "</UL>";

################################################################
#### Return fields
print "<HR>\n";
print "<B><A HREF='help.roc-stats.html#output_format'>Return</A></B>&nbsp;<UL>\n";
my $i = 0;
foreach my $field (@output_fields) {
  $i++;
  print $query->checkbox(-name=>$field,
			 -checked=>$default{$field},
			 -label=>'');
  print "&nbsp;<A HREF='help.roc-stats.html#",$field,"'><B>", $field_description{$field}, "</B></A>\n";
#  if ($field eq "graphs"){
#    print "<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Image format (if graphs is checked)&nbsp;\n";
#    print $query->popup_menu(-name=>'img_format',
#			     -Values=>['png',
#				       'jpg',
#				       'gif',
#				       'eps',
#				       'pdf'],
#			     -default=>$default{img_format});
#  }
print "<BR>\n";
}
print "</ul><p>\n";

################################################################
#### Return graphs
#print "<p><B><A HREF='help.roc-stats.html#graphs'>Graphs parameters (if graphs output is checked)</A></B>&nbsp;<br>\n";

################################################################
### send results by email or display on the browser
print "<HR>\n";
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
my $demo_data=`grep -v "^;" $ENV{RSAT}/public_html/demo_files/roc-stats_demo_regulonDB_cg_stringcoex_inter-Q.tab`;
# "0.95	pos
# 0.85	neg
# 0.75	pos
# 0.65	pos
# 0.55	pos
# 0.45	neg
# 0.35	pos
# 0.35	neg
# 0.25	neg
# 0.15	pos
# 0.05	neg
# ";
#print "$ENV{RSAT}/public_html/data/demo_files/roc-stats_demo_regulonDB_cg_stringcoex_inter-Q.tab\n<BR>".$demo_data;
print "<TD><B>";
print $query->hidden(-name=>'data',-default=>$demo_data);
print $query->hidden(-name=>'sc_col',-default=>'3');
print $query->hidden(-name=>'status_col',-default=>'6');
print $query->hidden(-name=>'pos',-default=>"Q.and.R\nR.not.Q");
print $query->hidden(-name=>'neg',-default=>"Q.not.R");
print $query->hidden(-name=>'AUC',-default=>"checked");
print $query->hidden(-name=>'demo_comment',-default=>1);
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B><A HREF='help.roc-stats.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_roc.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:rekins\@bigre.ulb.ac.be,jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

################################################################
#### print Header

sub NeAT_header {
  my $css_body_class = "form";
  my ($title) = shift;
  $title =~ s/\"//g;
  $title =~ s/\'//g;
  if (scalar @_ > 0) {
    $css_body_class = shift;
  }


#  print &html_header();
  print $query->header();
  print sorttable_script();
  ### print the header of the result page
  print $query->start_html(-title=>"Network Analysis Tools : $title",
			   -class => "$css_body_class",
			   -author=>'jacques.van.helden@ulb.ac.be',
			   -style => { 	-src => "$ENV{rsat_www}/main.css",
                             	       	-type => 'text/css',
                             		-media => 'screen' });
  print "<H3 ALIGN='center'><A HREF='$ENV{rsat_www}/NeAT_home.html'>Network Analysis Tools</A> - $title</H3>";
}
