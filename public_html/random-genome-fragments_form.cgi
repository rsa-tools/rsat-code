#!/usr/bin/perl
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
$default{organism} = "Saccharomyces cerevisiae";
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
Program developed by <a href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</a>\n";
print "</CENTER>";


## demo description
print $default{demo_descr1};
print $default{demo_descr2};

print $query->start_multipart_form(-action=>"random-genome-fragments.cgi");


#### fragments
print "<fieldset><legend><b><a href='help.random-genome-fragments.html#fragments'>Fragment size and lengths</a></b></legend>";


################################################################
## Option 1: specify a number of fragments of fixed size
print "<input type='radio' name='fragment_sizes' value='fixed' $checked{'fixed'}/>";

# Length of fragments
print "<b> <A HREF='help.random-genome-fragments.html#l_sequence_length'>Fixed fragment lengths</A>&nbsp;</B>\n ";
print $query->textfield(-name=>'frag_length',
			-default=>$default{frag_length},
			-size=>5);
print "bases\n";
print "&nbsp;"x5;

# number of fragments
print "<b><A HREF='help.random-genome-fragments.html#r_repetitions'>Number of fragments</A>&nbsp;</B>\n";
print $query->textfield(-name=>'frag_nb',
			-default=>$default{frag_nb},
			-size=>5);

################################################################
## Option 2: use a file (bed or fasta) as template
print "<p><input type='radio' name='fragment_sizes' value='template' $checked{'file'}/>\n";
#print "<b><A HREF='help.random-genome-fragments.html#lf_length_file'>Template file (same nb of fragments, same lengths)</b></a>";

# print "<ul>";

# ## Text area to copy-paste the sequence
# $sequenceChoiceString .=  "Paste your sequence in the box below<BR>\n";
# $sequenceChoiceString .=  $query->textarea(-name=>'sequence',
# 					   -default=>$default{sequence},
# 					   -rows=>4,
# 					   -columns=>55);
# $sequenceChoiceString .=  "<BR>\n";

# print "<br>", $query->filefield(-name=>'uploaded_file',
# 				-default=>'',
# 				-size=>45,
# 				-maxlength=>200);

# print "</ul>";

# # print "<ul>\n";
# # print "<p><input type='radio' name='template_format' value='bed' $checked{'ODfragment_sizes_file'}/><b><A HREF='help.random-genome-fragments.html#lf_length_file'>Genomic coordinates (bed format)\n</b></a>";
# # print "<p><input type='radio' name='template_format' value='fasta' $checked{'ODfragment_sizes_file'}/><b><A HREF='help.random-genome-fragments.html#lf_length_file'>Sequences (fasta format)\n</b></a>";
# # print "</ul>\n";
# # print "<p>\n";

# Template file (bed coordinates or fasta sequences)
#print "<p/><b> OR </b> <p/>";
print "<b><A HREF='help.random-genome-fragments.html#lf_length_file'>Use template file (bed coordinates or fasta sequences): </a></b><br/> \n";
print "<div style='padding-left:30px'>";
&MultiSequenceChoice("Paste template data",1, "(bed coordinates or fasta sequences).");

print "<br><b>Template format</b> (bed coordinates or fasta sequences)&nbsp;\n";
print $query->popup_menu(-name=>'template_format',
			 -Values=>['bed',
				   'fasta',
				   'lengths',
			 ],
			 -default=>$default{template_format});

print "</div>";

print "</fieldset><p/>";


#### Organisms
print "<fieldset>";
print "<legend><b><a href='help.random-genome-fragments.html#organism'>Organism </a></b></legend>";


print "<P/>\n";
print "<input type='radio' name='org_select' value='rsat_org' $checked{'rsat_org'}/>";
print "<b>Local RSAT </b>"; &OrganismPopUp();


################################################################
## DEBUGGING TO BE DONE
## 
## JvH temporarily disactivated this option (2016-11-06) because the
## popup menu is empty on rsat-tagc.univ-amu.fr.
##
# print "<input type='radio' name='org_select' value='ensembl_org' $checked{'ensembl_org'}/>";
# print "<b>Ensembl </b>"; &OrganismPopUpEnsembl();
# print "<P/>\n";

print "</fieldset><p/>";

#### Output
print "<fieldset>
<legend><b><a href='help.random-genome-fragments.html#output_format'>Output</a></b></legend>";
print "<P/>\n";

print "<INPUT TYPE='radio' NAME='outputformat' VALUE='outputseq' $checked{'outputseq'}/>";
print "<b>Sequences in fasta format (only for RSAT organisms)</b>&nbsp;&nbsp;";

### Repeat masking
print $query->checkbox(-name=>'rm',
  		       -checked=>$default{rm},
  		       -label=>'');
print "&nbsp;<A HREF='help.retrieve-seq.html#rm'><B>Mask repeats</B></A>";
print "&nbsp;<A HREF='help.retrieve-seq.html#rm_list'><B>(only valid for organisms with annotated repeats)</B></A>";
print "<BR>\n";
print "<P/>\n";

### Coordinates
print "<INPUT TYPE='radio' NAME='outputformat' VALUE='outputcoord' $checked{'outputcoord'}/>";

print "<b>Genomic coordinates. <a href='help.random-genome-fragments.html#output_format'>Output format</a> </B>&nbsp;\n";
print $query->popup_menu(-name=>'coord_format',
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
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
## data for the demo with RSAT organism
my $descr1="<H4>Comment on the demonstration example for RSAT organism : </H4><blockquote class ='demo'>
In this demonstration, we calculate random fragments in the genome sequence of Saccharomyces cerevisiae.
We use a set of template sequences, to produce the same number of random fragments, of the same lengths as in the template.<p/>
The program will return the sequences of these fragments, in fasta format.
</blockquote>";
my $demo_seq=`cat demo_files/MET_up800-noorf.fasta`;

print $query->start_multipart_form(-action=>"random-genome-fragments_form.cgi");
print "<TD><B>";
$query->delete_all();
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'sequence1',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format1',-default=>'fasta');
print $query->hidden(-name=>'fragment_sizes',-default=>'file');
print $query->hidden(-name=>'template_format',-default=>'fasta');
print $query->hidden(-name=>'org_select',-default=>'rsat_org');
print $query->hidden(-name=>'organism',-default=>'Saccharomyces_cerevisiae');
print $query->hidden(-name=>'outputformat',-default=>'outputseq');
print $query->submit(-label=>"DEMO RSAT organism");
print "</B></TD>\n";
print $query->end_form;

################################################################
## Demo with a bed file as template
my $descr3="<H4>Comment on the demonstration example for bed template : </H4><blockquote class ='demo'>
In this demonstration, we use as template a bed file specifying genomic coordinates of ChIP-seq peaks. 
The coordinates correspond to ChIP-seq peas for the transcription factor CEBPA in the mm9 assembly of the Mus musculus genome (Schmidt et al, 2010). 
The program random-genome-fragments selects random genomic coordinates having the same sizes as the template peaks.
</blockquote>";
my $demo_bed_url= $ENV{rsat_www}."/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed";
print $query->start_multipart_form(-action=>"random-genome-fragments_form.cgi");
print "<TD><B>";
$query->delete_all();
print $query->hidden(-name=>'demo_descr3',-default=>$descr3);
print $query->hidden(-name=>'sequence_url1',-default=>$demo_bed_url);
print $query->hidden(-name=>'sequence_format1',-default=>'bed');
print $query->hidden(-name=>'fragment_sizes',-default=>'file');
print $query->hidden(-name=>'template_format',-default=>'bed');
#print $query->hidden(-name=>'org_select',-default=>'ensembl_org');
#print $query->hidden(-name=>'organism',-default=>'Mus_musculus'); ## PROBLEM HERE: we need mm9, whereas Ensembl does not provide support for the older versions
print $query->hidden(-name=>'org_select',-default=>'rsat_org');
print $query->hidden(-name=>'organism',-default=>'Mus_musculus_GRCm37'); ## PROBLEM HERE: on metazoa  the version GRCm37 is currently not supported
print $query->hidden(-name=>'outputformat',-default=>'outputseq');
print $query->submit(-label=>"DEMO bed ChIP-seq peaks");
print "</B></TD>\n";
print $query->end_form;


################################################################
## data for the demo with Ensembl organism
##
## BUG: JvH temporarily disactivates this demo because
## supported-ensembl-organisms does not work when called from
## RSA2.cgi.lib, whereas it does work from the command line.
##
# my $descr2="<H4>Comment on the demonstration example for Ensembl organism : </H4><blockquote class ='demo'>
# In this demonstration, we calculate the coordinates of randomly-chosen fragments in the genome sequence of Homo sapiens. We would like 10 fragments of 100bp. <p/>
# The program will return the coordinates of these fragments, in BED format, that can be then used to extract the sequences with tools of
# sequence providers (UCSC, Galaxy, Ensembl).
# </blockquote>";
# print $query->start_multipart_form(-action=>"random-genome-fragments_form.cgi");
# print "<TD><B>";
# $query->delete_all();
# print $query->hidden(-name=>'demo_descr2',-default=>$descr2);
# print $query->hidden(-name=>'frag_length',-default=>'100');
# print $query->hidden(-name=>'frag_nb',-default=>'10');
# print $query->hidden(-name=>'org_select',-default=>'ensembl_org');
# print $query->hidden(-name=>'outputformat',-default=>'outputcoord');
# print $query->hidden(-name=>'coord_format',-default=>'bed');
# print $query->hidden(-name=>'organism_ens',-default=>'Homo_sapiens');
# print $query->submit(-label=>"DEMO Ensembl organism");
# print "</B></TD>\n";
# print $query->end_form;

print "<TD><B><A HREF='help.random-genome-fragments.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:morgane\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";


print $query->end_html;

exit(0);

