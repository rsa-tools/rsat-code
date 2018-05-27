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
print "Given a set of polymorphic variations, dbSNP IDs or varBed variants file, retrieve the corresponding variants and their flanking sequences. Optionally, the user can input a set of genomic coordinates in bed format, and the tool will retrieve dbSNP annotated variants within them, together with their flanking sequences. The output of this tool can then be used in  <a href='variation-scan_form.cgi'>variation-scan</a> to assess the impact of genetic variants on trasncription factor binding sites. <\P>
Locally installed variants were directly downloaded from <a target='_blank' href='http://www.ensembl.org/index.html'>ensembl</a>, specific versions for each installed genome can be consulted <a target='_blank' href='data/supported_variation_info.tab'>here</a>, updates can be done on demand. <\P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr//'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/index.html'>Alejandra Medina-Rivera</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://liigh.unam.mx/amedina/people.html'>Walter Santana</a><sup>cit</sup>\n";
print ", Jeremy Delerce<sup>ci</sup>\n";
print "</CENTER>";

################################################################
### display the form only if the organisms on the curent server
###  are coherent with this tool, otherwise, display an info message

if ($ENV{variations_tools} == 0){

print "<font color='#DD0000'>Sorry, this tool is not compatible with the organisms supported on this server.</font>\n";

print $query->end_html;

exit(0);

}

################################################################
### formheader

#print $default{demo_descr1};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"retrieve-variation-seq.cgi");

#### Select organims to retrieve variants sequences from
#print "<div id='organisms_popup' style='display:none'>";
#print "&nbsp;"x0, &OrganismPopUpString();
#print "<p>\n</div>";
#print "<div id='variations_popup'>";
#print &OrganismPopUpString("supported"=>"variations", "bg_org"=>"1");
#print "<p></div>";
### Query variants
### Variants can be input as a list of rs numbers, rsa variation file or bed regiones
### from where variants annotated in ensembl variations will be extracted.

print "<p>";
print "<B>Input a list of dbSNP variation IDs (rsID), a set of variants in varBed format, or genomic regions in bed format</B>&nbsp;";


print "<BR>\n";
print "<UL>\n";
 if ($variants_file = $query->param('variants_file')) {
    ## Variants file is already on the server machine
    ## (piped from a previous script)
    $variants_url = $variants_file;
    $variants_url =~ s|$ENV{RSAT}/public_html|$ENV{rsat_www}|;
    $variantsChoiceString .=  "<a href=$variants_url>";
    $variantsChoiceString .=  " transferred from previous query<BR>\n";
    $variantsChoiceString .=  "</a>";
    $variants_format = $query->param(variants_format);
    $variantsChoiceString .=  "<INPUT type='hidden' NAME='variants_format' VALUE='$variants_format'>\n";
    $variantsChoiceString .=  "<INPUT type='hidden' NAME='variants_file' VALUE='$variants_file'>\n";
    print $variantsChoiceString ;

}else{

    print $query->textarea(-name=>'input', -id=>'input',
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

}


### Input type
print "<B>Input format</B>&nbsp;";
    print $query->popup_menu(-name=>'input_type', -id=>'input_type',
			     -Values=>['varBed','id','bed'],
			     -default=>$default{input_type});
print "<\p>";

#### Select organims to retrieve variants sequences from
print "<div id='organisms_popup' style='display:none'>";
print "&nbsp;"x0, &OrganismPopUpString();
print "<p>\n</div>";
print "<div id='variations_popup'>";
print &OrganismPopUpString("supported"=>"variations", "bg_org"=>"1");
print "<p></div>";

### Lenght of the sequences surranding the variant
print "<B>Length of flanking sequence on each side of the variant</B>&nbsp;\n";
print $query->textfield(-name=>'mml', -id=>'mml',
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

print "<TD><B>";
###############

## Data for demo
$demo_rsat_var_file=$ENV{RSAT}."/public_html/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.varBed";
$demo_rsat_var = "";
open(my $fh, $demo_rsat_var_file);
while (my $row = <$fh>){
    chomp $row;
    $demo_rsat_var .= $row . "\\n";
}

print '<script>

$(document).ready(function(){
    $("#input_type").change(function(){
        if($("#input_type").val() == "varBed"){
            document.getElementById("organisms_popup").style.display = "block";
            document.getElementById("variations_popup").style.display = "none";
            $("#organism_bg").val("");
        }else{
            document.getElementById("organisms_popup").style.display = "none";
            document.getElementById("variations_popup").style.display = "block";
            $("#organism").val("");
        }
    });
});
function setDemo(demo_rsat_var){
    $("#reset").trigger("click");

    var descr1 = "<H4>Comment on the demonstration :</H4> \
            <blockquote class =\"demo\"> \
            <p>In this demonstration, we retrieve the sequence of genetic variants.</p> \
            <p> The genetic variants used in this example were collected by Weirauch et al (2014, Cell 158:1431-1443).\
            These variants had been reported in previous publications as affecting transcription factor binding. </p>\n \
            </blockquote>";

    demo_descr.innerHTML = descr1;
    document.getElementById("organisms_popup").style.display = "block";
    document.getElementById("variations_popup").style.display = "none";
    $("#organism_name").val("Homo sapiens GRCh37");
    $("#organism").val("Homo_sapiens_GRCh37");
    $("#input").val(demo_rsat_var);
    $("#input_type").val("varBed");
    $("#mml").val("30");
}
</script>';
print '<button type="button" onclick="setDemo('. "'$demo_rsat_var'" .')">DEMO</button>';

print "</B></TD>\n";

print $query->end_form;

print "<td><b><a href='sample_outputs/retrieve-variation-seq_demo20180226.varSeq'>[Sample Output]</a></B></TD>\n";
print "<TD><B><A HREF='help.retrieve-variation-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "<br><br><font size=1 color=\"grey\" ><small>AMR and WS are supported by a PAPIIT-UNAM (IA206517) grant.</small></font>";

print $query->end_html;

exit(0);
