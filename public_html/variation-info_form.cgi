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
$default{demo_descr1} = "";
$default{organism} = "Homo_sapiens_GRCh37";
$default{input_type}="gvf";
$default{out_type}="rsat-var";
$default{mml}=30 ; ## Length of the sequence sorounding the variant, 
                   ## has to be consistent with the longest matrix to be used
### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("Variation information", 'form');

### head
print "<CENTER>";
print "Get information (position, variants) about genomic variations, given a set of genomic regions (return overlapping variations) or a list of variation  IDS (return the list of variations maching query IDs).<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.epernicus.com/am27'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", Jeremy Delerce<sup>ci</sup>\n";
print ", Yvon Mbouamboua<sup>t</sup>\n";
print "</CENTER>";

print $default{demo_descr1};

print $query->start_multipart_form(-action=>"variation-info.cgi");


#print "<FONT FACE='Helvetica'>";

#### Select organims to retrieve variants sequences from

print "&nbsp;"x0, &OrganismPopUpString();
print "<p>\n";


### Query variants
### Variants can be input as a list of rs numbers, rsa variation file or bed regiones
### from where variants annotated in ensembl variations will be extracted.

print "<p>";
print "<B>Variations IDs or genomic regions of interest</B>&nbsp;";


print "<BR>\n";
print "<UL>\n";

print $query->textarea(-name=>'input',
		       -default=>$default{input},
		       -rows=>6,
		       -columns=>65);
### Option to upload a file with variant information (IDs  or 
### genomic regions in bed format)
print "<BR>Upload variants or regions<BR>\n";
print $query->filefield(-name=>'uploaded_file',
			-default=>'',
			-size=>45,
			-maxlength=>200);
print "</UL>\n";
print "<BR>\n";

 ## Option to fetch sequence file from an URL
print "&nbsp;"x3;
print "URL of variants file or genomic regions available on a Web server (e.g. Galaxy).<BR>\n";
print $query->textfield(-name=>'variants_url',
			-default=>"",
			-size=>62);
print "<br>\n";

### Input type
print "<B>Input format</B>&nbsp;";
print $query->popup_menu(-name=>'input_type',
			 -Values=>['bed', 'id'],
			 -default=>$default{input_type});
print "<\p>";


print "<\p>";
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
$descr1 .= "<blockquote class ='demo'>";
$descr1 .= "<p>In this demonstration, we retrieve information about a set of variations specified by providing a list of IDs as query. </p>\n\n";
$descr1 .= "<p>The genetic variations used in this example were taken from Weireauch, et al (Cell, 2014). ";
$descr1 .= "These authors collected from the litterature a set of variants shown to affect transcription factor binding. </p>\n";
$descr1 .= "</blockquote>";

print $query->start_multipart_form(-action=>"variation-info_form.cgi");
## Data for demo
$demo_gvf_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt";
$demo_gvf_var=`cat $demo_gvf_file` ;


print "<TD><B>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'organism',-default=>"Homo_sapiens_GRCh37");
print $query->hidden(-name=>'input',-default=>"$demo_gvf_var");
print $query->hidden(-name=>'input_type',-default=>"id");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.variation-info.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";


print $query->end_html;

exit(0);

