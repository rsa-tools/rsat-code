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

print $default{demo_descr1};

print $query->start_multipart_form(-action=>"variation-info.cgi");


#print "<FONT FACE='Helvetica'>";

my @org_variation=();

my $data_rsat=join("/",$ENV{RSAT},"data") ;
my $supported_variation_organims_file=join ("/",$data_rsat,"supported_organisms_variation.tab");

if (-e $supported_variation_organims_file){
    my ($var_org) = &OpenInputFile($supported_variation_organims_file);
    while(<$var_org>){
	chomp;
	next unless (/\S/) ; # skip empty rows
	next if (/^;/); # skip comment lines
	next if (/^\#/); # Skip header line	
	my $org=$_ ;
	push (@org_variations, $org) ;
	
    }
    print "&nbsp;"x0, &OrganismPopUpString(@org_variations);
    
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
print $query->textfield(-name=>'input_url',
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
## Description and data for demo 1: by SNP IDs
print $query->start_multipart_form(-action=>"variation-info_form.cgi");
$descr1 .= "<blockquote class ='demo'>";
$descr1 .= "<p>In this demonstration, we retrieve information about a set of variations specified by providing a list of IDs as query. </p>\n\n";
$descr1 .= "<p>The genetic variations used in this example were taken from Weireauch, et al (Cell, 2014). ";
$descr1 .= "These authors collected from the litterature a set of variants shown to affect transcription factor binding. </p>\n";
$descr1 .= "</blockquote>";
$demo1_gvf_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt";
$demo1_gvf_var=`cat $demo1_gvf_file` ;
print "<TD><B>";
print $query->hidden(-name=>'demo1_descr1',-default=>$descr1);
print $query->hidden(-name=>'organism',-default=>"Homo_sapiens_GRCh37");
print $query->hidden(-name=>'input',-default=>"$demo1_gvf_var");
print $query->hidden(-name=>'input_type',-default=>"id");
print $query->submit(-label=>"DEMO 1: by SNP IDs");
print "</B></TD>\n";
print $query->end_form;


################################################################
## Data for demo 2
print $query->start_multipart_form(-action=>"variation-info_form.cgi");
$descr2 = "<p>Genomic regions corresponding to cis-regulatory modules";
$descr2 .= " characterized by chip-seq peaks for 4 transcription factors (HFN6, FOXA1, CBPalpha, HNF4alpha),";
$descr2 .= " conserved across 5 mammalian species (Human, macaque, rat, mouse, doc). ";
$descr2 .= " Source: <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/25279814'>Ballester et al. (2013). eLife.</a></p>";
$descr2 .= "<p><font color='red'><b>Warning</b>: this demo takes several minutes because it searches for variants ni 1600 genomic regions. Email output is recommended.</font></b>";
$demo2_url= $ENV{rsat_www}."/demo_files/Ballester_etal_elife_2014_module_beyondprimates_conserved_hg18_lift_to_hg19.bed";
print "<TD><B>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr2);
print $query->hidden(-name=>'organism',-default=>"Homo_sapiens_GRCh37");
print $query->hidden(-name=>'input_url',-default=>"$demo2_url");
print $query->hidden(-name=>'input_type',-default=>"bed");
print $query->hidden(-name=>'output',-default=>"email"); ## Email required for this demo because it takes a while
print $query->submit(-label=>"DEMO 2: by genomic regions (peaks)");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.variation-info.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "<br><br><font size=1 color=\"grey\" ><small>AMR and WS are supported by a PAPIIT-UNAM (IA206517) grant.</small></font>";


print $query->end_html;

exit(0);

