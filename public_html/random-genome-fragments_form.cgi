#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program random-genome-fragments
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

### Read the CGI query
$default{demo_descr1} = "";
$default{demo_descr2} = "";

### default values for filling the form
$default{organism} = "";
$default{organism_ens} = "Homo_sapiens";
$default{frag_nb} = 20;
$default{frag_length} = 100;
$default{org_select}="rsat_org";
$default{fragment_sizes}="fixed";
$default{template_format}="bed";
$default{outputformat}="outputcoord";
$default{coord_format} = "bed";
$default{'sequence_url'.1} = "";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

## radio button checked values (to be placed after changing default values !!)
$checked{$default{org_select}}="CHECKED";
$checked{$default{outputformat}}="CHECKED";
$checked{$default{fragment_sizes}} = "CHECKED";

################################################################
### print the form ###
&RSA_header("random genome fragments", "form");
&ListParameters() if ($ENV{rsat_echo} >=2);

### head
print "<CENTER>";
print "Select a set of fragments with random positions in a given genome, and return their coordinates and/or sequences.<P>
Program developed by <a href='http://morgane.bardiaux.fr/'>Morgane Thomas-Chollier</a>\n";
print "</CENTER>";


## demo description
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"random-genome-fragments.cgi", -id=>"form");


#### fragments
print "<fieldset><legend><b><a class='iframe' href='help.random-genome-fragments.html#fragments'>Random fragments </a></b></legend>";

################################################################
## Option 1: specify a number of fragments of fixed size
print "<input type='radio' name='fragment_sizes' id='fragment_sizes' value='fixed' $checked{'fixed'}/>";

# Length of fragments
print "<b> <A class='iframe' HREF='help.random-genome-fragments.html#l_sequence_length'>Fixed fragment lengths</A>&nbsp;</B>\n ";
print $query->textfield(-name=>'frag_length', -id=>'frag_length',
			-default=>$default{frag_length},
			-size=>5);
print "bases\n";
print "&nbsp;"x5;

# number of fragments
print "<B><A class='iframe' HREF='help.random-genome-fragments.html#r_repetitions'>Number of fragments</A>&nbsp;</B>\n";
print $query->textfield(-name=>'frag_nb', -id=>'frag_nb',
			-default=>$default{frag_nb},
			-size=>5);

################################################################
## Option 2: use a file (bed or fasta) as template
print "<p><input type='radio' name='fragment_sizes' value='template' $checked{'file'}/>\n";

print "<b><A class='iframe' HREF='help.random-genome-fragments.html#lf_length_file'>Use template file (bed coordinates or fasta sequences): </a></b><br/> \n";
print "<div style='padding-left:30px'>";
&MultiSequenceChoice("Paste template data",1, "(bed coordinates or fasta sequences).");
print "<br><b>Template format</b> (bed coordinates or fasta sequences)&nbsp;\n";
print $query->popup_menu(-name=>'template_format', -id=>'template_format',
			 -Values=>['bed',
				   'fasta',
				   'lengths',
			 ],
			 -default=>$default{template_format});
print "</div>";
print "</fieldset><p/>";


#### Organisms
print "<fieldset>";
print "<legend><b><a class='iframe' href='help.random-genome-fragments.html#organism'>Organism </a></b></legend>";


print "<P/>\n";
print "<INPUT TYPE='radio' id='org_select_rsat' NAME='org_select' VALUE='rsat_org' $checked{'rsat_org'} style='display:none'/>";
#print "<b>Local RSAT </b>";
&OrganismPopUp();


#print "<INPUT TYPE='radio' id='org_select_ensembl' NAME='org_select' VALUE='ensembl_org' $checked{'ensembl_org'}/>";
#print "<b>Ensembl </b>"; print &OrganismPopUpEnsembl();
print "<P/>\n";

print "</fieldset><p/>";

#### Output
print "<fieldset>
<legend><b><a class='iframe' href='help.random-genome-fragments.html#output_format'>Output</a></b></legend>";
print "<P/>\n";

print "<INPUT TYPE='radio' id='outputformat_seq' NAME='outputformat' VALUE='outputseq' $checked{'outputseq'}/>";
print "<b>Sequences in fasta format</b>&nbsp;&nbsp;";

### Repeat masking
print $query->checkbox(-name=>'rm',
  		       -checked=>$default{rm},
  		       -label=>'');
print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#rm'><B>Mask repeats</B></A>";
print "&nbsp;<A class='iframe' HREF='help.retrieve-seq.html#rm_list'><B>(only valid for organisms with annotated repeats)</B></A>";
print "<BR>\n";
print "<P/>\n";

### Coordinates
print "<INPUT TYPE='radio' id='outputformat_coord' NAME='outputformat' VALUE='outputcoord' $checked{'outputcoord'}/>";

print "<b>Genomic coordinates. <a class='iframe' href='help.random-genome-fragments.html#output_format'>Output format</a> </B>&nbsp;\n";
print $query->popup_menu(-name=>'coord_format', -id=>'coord_format',
			 -Values=>['ft',
				   'bed',
				   'bed3col',
				   'great'
			 ],
			 -default=>$default{coord_format});


print "</fieldset><p/>";


### send results by email or display on the browser
print "<P>\n";
&SelectOutput("server");

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
## data for the demo with RSAT organism

my $demo_seq = "";
open(my $fh, "demo_files/MET_up800-noorf.fasta");
while(my $row = <$fh>){
    chomp $row;
    $demo_seq .= $row;
    $demo_seq .= "\\n";
}

print '<script>
function setDemo1(demo_seq){
    $("#form")[0].reset();
    descr = "<H4>Comment on the demonstration example for RSAT organism : </H4><blockquote class =\'demo\'>\
    In this demonstration, we calculate random fragments in the genome sequence of Saccharomyces cerevisiae.\
    We use a set of template sequences, to produce the same number of random fragments, of the same lengths as in the \template.<p/>\
    The program will return the sequences of these fragments, in fasta format.\
    </blockquote>";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;
    sequence1.value = demo_seq;
    $("input[name=fragment_sizes]").val(["template"]);    
    template_format.value = "fasta";
    $("#org_select_rsat").prop("checked", true);
    $("#organism").val("Saccharomyces_cerevisiae");
    $("#organism_name").val("Saccharomyces cerevisiae");
    $("#outputformat_seq").prop("checked", true);
}
</script>';


print "<TD><B>";

print '<button type="button" onclick="setDemo1('. "'$demo_seq'" .')">DEMO</button>';

print "</B></TD>\n";

################################################################
## data for the demo with a bed file as template

#print '<script>
#function setDemo3(demo_url){
#    $("#form")[0].reset();
#    descr = "<H4>Comment on the demonstration example for bed template : </H4><blockquote class =\'demo\'> In this demonstration, we use as template a bed file specifying genomic coordinates of ChIP-seq peaks.The coordinates correspond to ChIP-seq peas for the transcription factor CEBPA in the mm9 assembly of the Mus musculus genome (Schmidt et al, 2010). The program random-genome-fragments selects random genomic coordinates having the same sizes as the template peaks.</blockquote>";
    
#    demo_descr.innerHTML = descr;
#    demo.value = descr;
#    sequence_url1.value = demo_url;
#    $("input[name=fragment_sizes]").val(["template"]);
#    template_format.value = "bed";
#    frag_length.value = "100";
#    frag_nb.value = "20";
    
#    $("input[name=org_select]").val(["rsat_org"]);
    
#    $("#organism").val("Mus_musculus_GRCm37");
#    $("#organism_name").val("Mus musculus GRCm37");
#    $("#outputformat_seq").prop("checked", true);
#}
#</script>';

#my $demo_bed_url= $ENV{rsat_www}."/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed";

#print "<TD><B>";

#print '<button type="button" onclick="setDemo3('. "'$demo_bed_url'" .')">DEMO bed ChIP-seq peaks</button>';

#print "</B></TD>\n";

################################################################
## data for the demo with Ensembl organism
#print '<script>
#function setDemo2(demo_seq){
#    $("#form")[0].reset();
#    descr = "<H4>Comment on the demonstration example for Ensembl organism : </H4><blockquote class =\'demo\'>\
#    In this demonstration, we calculate the coordinates of randomly-chosen fragments in the genome sequence of Homo sapiens. We \would like 10 fragments of 100bp. <p/>\
#    The program will return the coordinates of these fragments, in BED format, that can be then used to extract the sequences \with tools of\
#    sequence providers (UCSC, Galaxy, Ensembl).</blockquote>";
#    $("input[name=fragment_sizes]").val(["template"]);
#    demo_descr.innerHTML = descr;
#    demo.value = descr;
#    sequence1.value = demo_seq;
    
#    frag_length.value = "100";
#    frag_nb.value = "20";
#    $("input[name=org_select]").val(["ensembl_org"]);
#	$("#organism").val("Saccharomyces_cerevisiae");
#    $("#organism_name").val("Saccharomyces cerevisiae");
#    $("#outputformat_coord").prop("checked", true);
#    $("#coord_format").val("bed");
#    $("#organism_ens").val("homo_sapiens");
#}
#</script>';


#print "<TD><B>";

#print '<button type="button" onclick="setDemo2('. "'$demo_seq'" .')">DEMO Ensembl organism</button>';

#print "</B></TD>\n";

print "<TD><B><A class='iframe' HREF='help.random-genome-fragments.html'>MANUAL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";


print $query->end_html;

exit(0);

