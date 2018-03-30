#!/usr/bin/perl
#### this cgi script fills the HTML form for the program retrieve-variation-seq
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;

require RSAT::organism;

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
print "Get information (position, variants) about genomic variations, given a set of genomic regions (return overlapping variations) or a list of variation  IDS (return the list of variations maching query IDs).<P>
Locally installed variants were directly downloaded from <a target='_blank' href='http://www.ensembl.org/index.html'>ensembl</a>, specific versions for each installed genome can be consulted here, updates can be done on demand. <\P>\n";

print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/index.html'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/people.html'>Walter Santana</a><sup>cit</sup>\n";
print ", Jeremy Delerce<sup>ci</sup>\n";
print ", Yvon Mbouamboua<sup>t</sup>\n";
print "</CENTER>";
 
################################################################
## Display the form only if this RSAT instance supports variation
## analysis.
&check_variation_tools();


################################################################
### formheader

#print $default{demo_descr1};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"variation-info.cgi");


#print "<FONT FACE='Helvetica'>";

#my @org_variation=();

my $data_rsat=join("/",$ENV{RSAT},"data") ;
my $supported_variation_organims_file=join ("/",$data_rsat,"supported_organisms_variation.tab");

if (-e $supported_variation_organims_file){
    print "&nbsp;"x0, &OrganismPopUpString('supported'=>'variations');
}
else {
  &RSAT::message::Warning("This RSAT site does not contain any genome with variations");
  exit();
}
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
print $query->textfield(-name=>'input_url', -id=>'input_url',
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
&SelectOutput('server');

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

print " <BR>\n<\p>";
################
$demo1_gvf_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt";
$demo1_gvf_var= "";

open(my $fh, $demo1_gvf_file);

while (my $row = <$fh>){
    chomp $row;
   $demo1_gvf_var .= $row . "\\n";
}


$demo2_url= $ENV{rsat_www}."/demo_files/Ballester_etal_elife_2014_module_beyondprimates_conserved_hg18_lift_to_hg19.bed";
my $demo_org = "Homo sapiens GRCh37";
my $org = $demo_org;
$org =~ s/\ /_/g;

print '<script>
function setDemo1(){
    $("#reset").trigger("click");
    
    var descr1 = "<blockquote class =\"demo\"> \
    <p>In this demonstration, we retrieve information about a set of variations specified by providing a list of IDs as query. </p>\n\n\
    <p>The genetic variations used in this example were taken from Weireauch, et al (Cell, 2014). \
    These authors collected from the litterature a set of variants shown to affect transcription factor binding. </p>\n \
    </blockquote>";
    
    demo_descr.innerHTML = descr1;
    
    $("#organism_name").val("'. $demo_org . '");
    $("#organism").val("' . $org . '");
    $("#input").val ("' . $demo1_gvf_var. '") ;
    $("#input_type").val("id");
}

function setDemo2(){
    $("#reset").trigger("click");
    
    var descr2 = "<blockquote class =\"demo\"> \
    <p>Genomic regions corresponding to cis-regulatory modules \
     characterized by chip-seq peaks for 4 transcription factors (HFN6, FOXA1, CBPalpha, HNF4alpha), conserved across 5 mammalian species (Human, macaque, rat, mouse, doc). Source: <a target=\'_blank\' href=\'http://www.ncbi.nlm.nih.gov/pubmed/25279814\'>Ballester et al. (2013). eLife.</a></p>\
     <p><font color=\"red\"><b>Warning</b>: this demo takes several minutes because it searches for variants in 1600 genomic regions. Email output is recommended.</font></b> \
    </blockquote>";
    
    demo_descr.innerHTML = descr2;
    
    $("#organism_name").val("'. $demo_org . '");
    $("#organism").val("' . $org . '");
    $("#input").val("");
    $("#input_url").val("'.$demo2_url.'");
    $("#input_type").val("bed");
}

</script>';
print '<button type="button" onclick="setDemo1('. "'$demo_gvf_var'" .')">DEMO 1: by SNP IDs</button>';
print '<button type="button" onclick="setDemo2()">DEMO 2: by genomic regions (peaks)</button>';

print "<td><b><a href='sample_outputs/variation-info_demo20180321.varBed'>[Sample Output]</a></B></TD>\n";
print "<TD><B><A HREF='help.variation-info.html'>MANUAL</A></B></TD>\n";

print "</TR></TABLE></UL></UL>\n";

print "<br><br><font size=1 color=\"grey\" ><small>AMR and WS are supported by a PAPIIT-UNAM (IA206517) grant.</small></font>";


print $query->end_html;

exit(0);

