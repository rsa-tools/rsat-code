#!/usr/bin/env perl
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

#local @supported_input_formats = sort(keys(%RSAT::feature::supported_input_format));
local @supported_input_formats = qw(bed bed3col dnapat ft gft gff gff3bed swembl galaxy_seq  ucsc_seq);

#local @supported_output_formats = sort(keys(%RSAT::feature::supported_output_format));
local @supported_output_formats = qw(bed bed3col dnapat ft gft gff gff3 great);

##my $input_formats = join (",",@supported_input_formats);

################################################################
### default values for filling the form
$default{output}="display";
$default{feature_format} = "dnapat";
$default{feature} = "";
$default{bed_coord} = "";
$default{input_format}="ft";
$default{output_format}="gff3";

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
&RSA_header("convert-features", "form");
print "<CENTER>";
print "Interconversions between formats of feature descriptions.<P>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

print $query->start_multipart_form(-action=>"convert-features.cgi", -id=>"form");


################################################################
#### Features
print "<hr>";

print "<B>Feature</B>\n";
print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";

#### feature format (pop-up menu)
print "<A class='iframe' HREF='help.convert-features.html'><B>Format</B></a>&nbsp;";
print  $query->popup_menu(-name=>'feature_format',-id=>'feature_format',
			 -Values=>[@supported_input_formats],
			 -default=>$default{input_format});
print "<br/>";

### text area to copy-paste the feature
print  "Paste your feature in the box below<BR>\n";
print $query->textarea(-name=>'feature',-id=>'feature',
		       -default=>$default{feature},
			   -rows=>4,
			-columns=>55);
						  
print  "<BR>\n";

### option to upload the feature file from the client machine 
print "Or select a file to upload<BR>\n";
print  $query->filefield(-name=>'uploaded_file',
			 -default=>'',
			 -size=>45,
			 -maxlength=>200);

print "<HR/>";

### change coordinates

print "<B>(Optional) Conversion from relative to genomic coordinates</B>\n<p/>";

  print "The file to convert (from the box above) must contain features which coordinates are <i>relative</i> to larger fragments. To transform these relative coordinates into <i>genomic</i> coordinates, enter below a BED file (zero-based) containing the genomic coordinates of these larger fragments." ;
  print "&nbsp;"x3, "<br>The 4th column of this BED file (feature name) must correspond to the name of the feature in the file to convert.<br/>";

  print $query->textarea(-name=>'bed_coord',-id=>'bed_coord',
		       -default=>$default{bed_coord},
			   -rows=>4,
			-columns=>55);
						  
print  "<BR>\n";
print "Or select a file to upload<BR>\n";
  print $query->filefield(-name=>'bed_file',
				      -size=>10);

print "<HR/>";

### Output bg format
print "<BR>";

print "<B><A class='iframe' HREF='help.convert-features.html'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',-id=>'output_format',
			 -Values=>[@supported_output_formats],
			 -default=>$default{output_format});
print "<BR/>\n";



################################################################
### send results by email or display on the browser
print "<p>\n";
&SelectOutput("display");

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 

################################################################
### data for the demo 
$demo="";
open($fh, "demo_files/convert-features_demo.pat");
while($row = <$fh>){
    chomp $row;
    $demo .= $row;
    $demo .= "\\n";
}
print '<script>
function setDemo(demo){
    $("#reset").trigger("click");
    feature.value = demo;
    bed_coord.value = "";
    $("#feature_format").val("dnapat");
    $("#output_format").val("ft");
}
</script>';

print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo'" .')">DEMO</button>';
print "</B></TD>\n";

#### demo 2
my $demo_file1= $ENV{RSAT}."/public_html/demo_files/seq_mm9_galaxy_matrix-scan.ft";
my $demo_file1_content="";
open($fh, $demo_file1);
while($row = <$fh>){
    chomp $row;
    $demo_file1_content .= $row;
    $demo_file1_content .= "\\n";
}
my $demo_file2= $ENV{RSAT}."/public_html/demo_files/seq_mm9_galaxy.bed";
my $demo_file2_content="";
open($fh, $demo_file2);
while($row = <$fh>){
    chomp $row;
    $demo_file2_content .= $row;
    $demo_file2_content .= "\\n";
}

print '<script>
function setDemo2(demo1, demo2){
    $("#reset").trigger("click");
    feature.value = demo1;
    bed_coord.value = demo2;
    $("#feature_format").val("ft");
    $("#output_format").val("bed");
}
</script>';
print "<TD><b>";
print '<button type="button" onclick="setDemo2('. "'$demo_file1_content'" .','."'$demo_file2_content'".')">DEMO genomic coordinates conversion</button>';
print "</B></TD>\n";


print "<TD><B><A class='iframe' HREF='help.convert-features.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);
