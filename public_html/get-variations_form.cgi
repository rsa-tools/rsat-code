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
$default{organism} = "";
$default{input_type}="gvf";
$default{out_type}="rsat-var";
$default{mml}=30 ; ## Length of the sequence surrounding the variant, 
                   ## has to be consistent with the longest matrix to be used
### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("Get variations", 'form');

### head
print "<CENTER>";
print "Get information (position, variants) about genomic variations, given a set of genomic regions (return overlapping variations) or a list of variation  IDS (return the list of variations maching query IDs).<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.epernicus.com/am27'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", Jeremy Delerce<sup>ci</sup>\n";
print ", Yvon Mbouamboua<sup>t</sup>\n";
print "</CENTER>";

#print $default{demo_descr1};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"get-variations.cgi");


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

print $query->textarea(-name=>'input', -id=>'input',
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
print $query->popup_menu(-name=>'input_type', -id=>'input_type',
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

print "<TD><B>";

### data for demo
$demo_gvf_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt";
$demo_gvf = "";
open(my $fh, $demo_gvf_file);
while (my $row = <$fh>){
    chomp $row;
    $demo_gvf .= $row . "\\n";
}

my $demo_org = "Homo sapiens GRCh37";
my $org = $demo_org;
$org =~ s/\ /_/g;
print '<script>
function setDemo(demo_gvf){
    $("#reset").trigger("click");
    var descr1 = "<blockquote class =\'demo\'> \
    <p>In this demonstration, we retrieve varian information using a list of IDs\n \
    <p> The genetic variants used in this example were collected by Weireauch, et al (Cell, 2014), these variants were reported in previous publications as affecting transcription factor binding. </p> \
    </blockquote>";
    
    demo_descr.innerHTML = descr1;
    
    $("#organism_name").val("'. $demo_org . '");
    $("#organism").val("' . $org . '");
    $("#input").val(demo_gvf);
    $("#input_type").val("id");
}
</script>';
print '<button type="button" onclick="setDemo('. "'$demo_gvf'" .')">DEMO</button>';


print "<TD><B><A HREF='help.get-variations.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_get-variations.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";


print $query->end_html;

exit(0);

