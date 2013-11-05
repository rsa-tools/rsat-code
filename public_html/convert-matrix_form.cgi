#!/usr/bin/perl
#### this cgi script fills the HTML form for the program convert-matrix
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
use RSAT::matrix;
use RSAT::MatrixReader;

### Read the CGI query
$query = new CGI;

local @supported_output_formats = sort(keys( %RSAT::matrix::supported_output_format));

################################################################
### default values for filling the form
$default{output}="display";
$default{matrix}="";
$default{matrix_file}="";
$default{matrix_format} = "tab";
$default{output_format} = "tab";
$default{counts}="checked";
$default{comments}="checked";
$default{consensus}="checked";
$default{frequencies}="";
$default{info}="";
$default{parameters}="";
$default{profile}="";
$default{weights}="";
$default{pseudo_counts}=1;
$default{header}="checked";
$default{margins}="";
#$default{links}="checked";
$default{max_profile}=10;
$default{decimals}=1;
$default{rc} = "";
$default{multiply} = 1;
$default{perm} = 0;
$default{pseudo_prior} = "pseudo_prior";
$checked{$default{pseudo_prior}} = "CHECKED";
$default{bg_pseudo} = "0.01";
$default{bg_format}="oligo-analysis";
$default{bg_method}="from_matrix";
$checked{$default{bg_method}} = "CHECKED";
$default{logo}="checked";
$default{error_bar}="checked";
$default{small_correc}="checked";
$default{stretch}="";

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
&RSA_header("convert-matrix", "form");
print "<CENTER>";
print "Convert different types of position-specific scoring matrices (PSSM), and calculate statistical parameters.<P>\n";
#print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"convert-matrix.cgi");

#print "<FONT FACE='Helvetica'>";

################################################################
#### Matrix specification
print "<hr>";

## Input matrix
&GetMatrix();
print "<hr>";

## Background model
my %bg_params =("from_matrix" => 1);
&GetBackgroundModel(%bg_params);

print "<br/><b>Note:</b> Only Bernoulli models are supported. Higher-order Markov models are converted into Markov 0 (Bernoulli).";
print "<hr>";

### Output matrix format
print "<br>";
print "<b><a href='help.convert-matrix.html#output_format'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
			 -Values=>[@supported_output_formats],
			 -default=>$default{output_format});
print "<BR>\n";

################################################################
## Output fields
print "<p><b><a href='help.convert-matrix.html#return'>Output fields</a></B>&nbsp;<br>\n";
my $i = 0;
foreach my $stat qw(counts frequencies weights info header margins consensus parameters profile comments) {
  print $query->checkbox(-name=>$stat,
			 -checked=>$default{$stat},
			 -label=>'');
  print "&nbsp;<A HREF='help.convert-matrix.html#",$stat,"'><B>", $stat, "</B></A>\n";
  print "<br>\n";
}

print $query->checkbox(-name=>'logo',
		       -checked=>$default{logo},
		       -label=>'');
print "&nbsp;<A HREF='help.convert-matrix.html#logo'><B>", "logo", "</B></A>&nbsp;(using <a target='_blank' href='http://weblogo.berkeley.edu/'>Weblogo</a>)\n";
print "&nbsp;&nbsp; (<b>options</b>:";
print $query->checkbox(-name=>'error_bar',
		       -checked=>$default{error_bar},
		       -label=>'Error bar');

print '&nbsp;'x3, $query->checkbox(-name=>'small_correc',
		       -checked=>$default{small_correc},
		       -label=>'Small sample correction');

print '&nbsp;'x3, $query->checkbox(-name=>'stretch',
		       -checked=>$default{stretch},
		       -label=>'Stretching of logos to entire length'); 
print ")<br>\n";

print "<br/>";
print "<A HREF='help.convert-matrix.html#decimals'><B>score decimals</B></A>\n";
print $query->popup_menu(-name=>'decimals',
			 -Values=>['0',
				   '1','2'],
			 -default=>$default{decimals});

#### Compute reverse complement
print "<BR>\n";
print $query->checkbox(-name=>"rc",
		       -checked=>$default{rc},
		       -label=>'');
print "<B><A HREF='help.convert-matrix.html#rc'>Compute reverse complement</A></b>\n";

#### Multiply counts
print "<BR>\n";
print "<B><A HREF='help.convert-matrix.html#multiply'>Multiply counts</A></b>\n";
print $query->textfield(-name=>'multiply',
			-default=>$default{multiply},
			-size=>2);
print " (convert frequency matrices into count matrices)\n";

#### permutations
print "<BR>\n";
print "<B><A HREF='help.convert-matrix.html#permutations'>Number of permutations</A></b>\n";
print $query->textfield(-name=>'perm',
			-default=>$default{perm},
			-size=>2);
print "<B>(returns 'counts' field only)</b>\n";

################################################################
### send results by email or display on the browser
print "<p>\n";
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
print $query->start_multipart_form(-action=>"convert-matrix_form.cgi");
my $demo_matrix=`cat demo_files/convert-matrix_demo_data.txt`;
print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->hidden(-name=>'input_format',-default=>'tab');
print $query->hidden(-name=>'output_format',-default=>'tab');
print $query->hidden(-name=>'header',-default=>"off");
print $query->hidden(-name=>'margins',-default=>"off");
#print $query->hidden(-name=>'info',-default=>"on");
#print $query->hidden(-name=>'weights',-default=>"on");
print $query->hidden(-name=>'parameters',-default=>"on");
#print $query->hidden(-name=>'links',-default=>"on");
print $query->hidden(-name=>'logo',-default=>"on");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.convert-matrix.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_PSSM.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

