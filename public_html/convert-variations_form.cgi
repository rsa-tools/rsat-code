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
$default{organism} = "Homo sapiens GRCh37";
$default{input_type}="gvf";
$default{out_type}="varBed";
$default{mml}=30 ; ## Length of the sequence surrounding the variant, 
                   ## has to be consistent with the longest matrix to be used
### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("Convert variation formats", 'form');

### head
print "<CENTER>";
print "Convert between different file formats that store genetic variation information. The most commonly used formats are:<a href='http://en.wikipedia.org/wiki/Variant_Call_Format'> VCF </a> and <a href='http://www.sequenceontology.org/resources/gvf_1.00.html'>GVF</a>, varBed format format presents several advantages for scanning variations with  matrices using <a href='variation-scan_form.cgi'>variation-scan</a> .<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/index.html'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/people.html'>Walter Santana</a><sup>cit</sup>\n<P>";
print "A new version of this tool replace a previous prototype version developed by Jeremy Delerce and Jacques van Helden";
print "</CENTER>";

## demo description
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"convert-variations.cgi");


#print "<FONT FACE='Helvetica'>";

#### Select organims to retrieve variants sequences from

print "&nbsp;"x0, &OrganismPopUpString();
print "<p>\n";


### Query variants
### Variants can be input as a list of rs numbers, rsa variation file or bed regiones
### from where variants annotated in ensembl variations will be extracted.

print "<p>";
print "<B>Variants to be converted</B>&nbsp;";


print "<BR>\n";
print "<UL>\n";

print $query->textarea(-name=>'input',-id=>'input',
		       -default=>$default{input},
		       -rows=>6,
		       -columns=>65);
### Option to upload a file with variant information (IDs, varBed file or 
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
print "URL of variants file available on a Web server (e.g. Galaxy).<BR>\n";
print $query->textfield(-name=>'variants_url',
			-default=>"",
			-size=>62);
print "<br>\n";

### Input type
print "<B>Input format</B>&nbsp;";
print $query->popup_menu(-name=>'input_type',-id=>'input_type',
			 -Values=>['varBed','vcf','gvf'],
			 -default=>$default{input_type});
print "<\p>";

### Out type
print "<B>Output format</B>&nbsp;";
print $query->popup_menu(-name=>'out_type',-id=>'out_type',
			 -Values=>['varBed','vcf','gvf'],
			 -default=>$default{out_type});
print "<\p>";
print "<BR>\n";

### send results by email or display on the browser
&SelectOutput("server");

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;


################
## Data for demo
$demo_gvf_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.gvf";
$demo_gvf_var= "";

open(my $fh, $demo_gvf_file);
while(my $row = <$fh>){
    chomp $row;
    $demo_gvf_var .= $row;
    $demo_gvf_var .= "\\n";
}

my $org = $default{organism};
$org =~ s/\ /_/g;
print '<script>
function setDemo(demo_gvf_var){
    $("#reset").trigger("click");
    $("#organism_name").val("' . $default{organism} .'");
    $("#organism").val("'. $org .'");
    descr = "<blockquote class =\'demo\'>";
    
    descr = descr + "<p>In this demonstration, we convert variants in <a href=\'http://www.sequenceontology.org/resources/gvf_1.00.html\'>GVF</a> format to varBed format.</p>\n \
    <p> The genetic variants used in this example were collected by Weireauch, et al (Cell, 2014), these variants were reported in previous publications as affecting transcription factor binding. </p>\n";
    
    descr = descr + "</blockquote>";
    demo_descr.innerHTML = descr;
    demo.value = descr;
    
    $("#organism").val("Homo_sapiens_GRCh37").trigger("chosen:updated");
    $("#input").val(demo_gvf_var);
    $("#input_type").val("gvf");
    $("#out_type").val("varBed");
}
</script>';


print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_gvf_var'" .')">DEMO</button>';
print "</B></TD>\n";


print "<TD><B><A class='iframe' HREF='help.convert-variations.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "<br><br><font size=1 color=\"grey\" ><small>AMR and WS are supported by a PAPIIT-UNAM (IA206517) grant.</small></font>";


print $query->end_html;

exit(0);

