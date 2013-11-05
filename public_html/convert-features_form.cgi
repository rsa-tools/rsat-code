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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

#local @supported_input_formats = sort(keys(%RSAT::feature::supported_input_format));
local @supported_input_formats = qw(bed dnapat ft galaxy_seq gft gff gff3bed swembl ucsc_seq);

#local @supported_output_formats = sort(keys(%RSAT::feature::supported_output_format));
local @supported_output_formats = qw(bed dnapat ft fasta gft gff gff3 great);

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

print $query->start_multipart_form(-action=>"convert-features.cgi");


################################################################
#### Features
print "<hr>";

print "<B>Feature</B>\n";
print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
print "<A HREF='help.convert-features.html'><B>Format</B></a>&nbsp;";


#### feature format (pop-up menu)
print  $query->popup_menu(-name=>'feature_format',
			 -Values=>[@supported_input_formats],
			 -default=>$default{input_format});
print "<br/>";
### text area to copy-paste the feature
print  "Paste your feature in the box below<BR>\n";
print $query->textarea(-name=>'feature',
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

  print $query->textarea(-name=>'bed_coord',
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

print "<B><A HREF='help.convert-features.html'>Output format</A></B>&nbsp;";
print $query->popup_menu(-name=>'output_format',
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
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo 

################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"convert-features_form.cgi");
$demo="; dna-pattern  -v -pl tmp/dna-pattern.2008_01_31.171040.pat -i tmp/dna-pattern.2008_01_31.171040.seq -format fasta -return sites -origin -0 -N 4 -noov -2str -subst 0
; Citation: van Helden et al. (2000). Yeast 16(2), 177-187.
; Input file           	tmp/dna-pattern.2008_01_31.171040.seq
; Input format         	fasta
; Pattern file         	tmp/dna-pattern.2008_01_31.171040.pat
; Search method        	regexp
; Threshold            	0
; Allowed substitutions	0
; Return fields
;                     	sites
; Patterns
; 	seq	id	score
; 	CACGTG	CACGTG	1
; 	CACGTT	CACGTT	1
; 
; Matching positions
; PatID	Strand	Pattern	SeqID	Start	End	matching_seq	Score
CACGTG	DR	CACGTG	PHO5	-253	-248	ctcaCACGTGggac	1.00
CACGTT	D	CACGTT	PHO5	-362	-357	ttagCACGTTttcg	1.00
CACGTT	R	CACGTT	PHO5	-724	-719	gggtCACGTTtctc	1.00
CACGTG	DR	CACGTG	PHO8	-534	-529	gggcCACGTGcagc	1.00
CACGTT	R	CACGTT	PHO8	-380	-375	atctCACGTTtctc	1.00
CACGTG	DR	CACGTG	PHO11	-283	-278	ttcaCACGTGggtt	1.00
CACGTT	D	CACGTT	PHO11	-416	-411	ttacCACGTTttcg	1.00
CACGTG	DR	CACGTG	PHO81	-344	-339	atggCACGTGcgaa	1.00
CACGTT	R	CACGTT	PHO81	-8	-3	tgCACGTTtatc	1.00
CACGTG	DR	CACGTG	PHO84	-436	-431	gttcCACGTGgacg	1.00
CACGTG	DR	CACGTG	PHO84	-414	-409	ccagCACGTGgggc	1.00
CACGTT	D	CACGTT	PHO84	-587	-582	tacgCACGTTggtg	1.00
CACGTT	R	CACGTT	PHO84	-262	-257	tacgCACGTTttta	1.00
; Job started	2008_01_31.171041
; Job done   	2008_01_31.171041
";
print "<TD><B>";
print $query->hidden(-name=>'feature',-default=>$demo);
print $query->hidden(-name=>'feature_format',-default=>'dnapat');
print $query->hidden(-name=>'output_format',-default=>"ft");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;

print "<TD><B>";
print $query->hidden(-name=>'matrix',-default=>$demo_matrix);
print $query->start_multipart_form(-action=>"convert-features_form.cgi");
my $demo_file1= $ENV{RSAT}."/public_html/demo_files/seq_mm9_galaxy_matrix-scan.ft";
my $demo_file1_content=`cat $demo_file1`;
my $demo_file2= $ENV{RSAT}."/public_html/demo_files/seq_mm9_galaxy.bed";
my $demo_file2_content=`cat $demo_file2`;
print "<TD><b>";
print $query->hidden(-name=>'feature',-default=>$demo_file1_content);
print $query->hidden(-name=>'bed_coord',-default=>$demo_file2_content);
print $query->hidden(-name=>'feature_format',-default=>'ft');
print $query->hidden(-name=>'output_format',-default=>"bed");
print $query->submit(-label=>"DEMO genomic coordinates conversion");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.convert-features.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);
