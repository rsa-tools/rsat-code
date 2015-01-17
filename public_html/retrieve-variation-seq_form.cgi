#!/usr/bin/perl
#### this cgi script fills the HTML form for the program retrieve-variation-seq
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
$default{organism} = "Homo_sapiens_GRCh37";
$default{input_type}="bed";
$default{mml}=30 ; ## Length of the sequence sorounding the variant, 
                   ## has to be consistent with the longest matrix to be used
### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("retrieve variation sequence", 'form');

### head
print "<CENTER>";
print "Given a set of IDs for polymorphic variations, retrieve the corresponding variants and their flanking sequences, in order to scan them wiht the tool <a href='variation-scan_form.cgi'>variation-scan</a> .<P>\n";
print "</CENTER>";


print $query->start_multipart_form(-action=>"retrieve-variation-seq.cgi");


#print "<FONT FACE='Helvetica'>";

#### Select organims to retrieve variants sequences from

print "&nbsp;"x0, &OrganismPopUpString();
print "<p>\n";


### Query variants
### Variants can be input as a list of rs numbers, rsa variation file or bed regiones
### from where variants annotated in ensembl variations will be extracted.

print "<p>";
print "<B>Variants or regions in bed format</B>&nbsp;";


print "<BR>\n";
print "<UL>\n";

print $query->textarea(-name=>'input',
		       -default=>$default{input},
		       -rows=>6,
		       -columns=>65);
### Option to upload a file with variant information (IDs, rsat-var file or 
### genomic regions in bed format)
print "<BR>Upload variants or regions<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);
print "</UL>\n";
print "<BR>\n";


### Input type
print "<B>Input format</B>&nbsp;";
print $query->popup_menu(-name=>'input_type',
			 -Values=>['rsat-var','id','bed'],
			 -default=>$default{input_type});
print "<\p>";
### Lenght of the sequences surranding the variant
print "<B>Length of sequence around the variant</B>&nbsp;\n";
print $query->textfield(-name=>'mml',
			-default=>$default{mml},
			-size=>5);
print "<BR>\n";

### send results by email or display on the browser
&SelectOutput("server");

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


################
## Data for demo

print $query->start_multipart_form(-action=>"retrieve-variation-seq_form.cgi");
## Data for demo
$demo_rsat_var_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set.rsat-var";
$demo_rsat_var=`cat $demo_rsat_var_file` ;


print "<TD><B>";
print $query->hidden(-name=>'organism',-default=>"Homo_sapiens_GRCh37");
print $query->hidden(-name=>'input',-default=>"$demo_rsat_var");
print $query->hidden(-name=>'input_type',-default=>"rsat-var");
print $query->hidden(-name=>'mml',-default=>"30");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.retrieve-variation-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_retrieve-variation-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";


print $query->end_html;

exit(0);

