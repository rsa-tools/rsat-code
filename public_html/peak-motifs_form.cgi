#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program peak-motifs
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

### Read the CGI query

################################################################
### default values for filling the form
$default{demo_descr} = "";
$default{lth_occ_sig}=0;
$default{uth_pval} = "1e-4";
$default{assembly} = "";
$default{oligo_length6}="checked";
$default{oligo_length7}="checked";
$default{oligo_length8}="";
$default{merge_lengths}="";
$default{"oligo-analysis"}="checked";
$default{"dyad-analysis"}="";
$default{"position-analysis"}="checked";
$default{'local-word-analysis'}="";
$default{'local-word-analysis_dyads'} ="";
$default{"position-analysis_dyads"} ="";
$default{'matrix-scan-quick'}="checked";
$default{compare_motif_db}="checked";
$default{title}="";
$default{max_seq_len}="500";
$default{top_sequences}="";
$default{nmotifs} = 5;
$default{origin} = "center";
$default{offset} = "0";
## Plot XY graph with R rather than GD
if ($ENV{R_supported}) {
    $default{r_plot} = "checked"; 
} else {
    $default{r_plot} = ""; 
}
$default{visualize}="none";


## motif database
$default{compare_motif_database}="jaspar_core_nonredundant_vertebrates";
$default{custom_motif_db_name}="custom_motif_collection";

### Replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
   if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
  if ($key eq "visualize"){
  	$checked{$query->param($key)} = "CHECKED";
  }
}

$checked{$default{visualize}} = "checked";



################################################################
### print the form ###

################################################################
### header
&RSA_header("peak-motifs", "form");

#print "<pre>default{visualize} ",$default{visualize}, "</pre>";
#print "<pre>checked{default{visualize}} ",$checked{$default{visualize}}, "</pre>";

print <<end_part_1;
<center>
<p>Pipeline for discovering motifs in massive ChIP-seq peak sequences.</p>
<!--
<p>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>:
<a target='_blank' href='http://jacques.van-helden.perso.luminy.univ-amu.fr/'>Jacques van Helden</a><sup>cit</sup>,
<a target='_blank' href='http://morgane.bardiaux.fr/'>Morgane Thomas-Chollier</a><sup>cit</sup>,
<a target='_blank' href='https://www.ulb.ac.be/rech/inventaire/chercheurs/4/CH10464.html'>Matthieu Defrance</a><sup>ci</sup>,
<a target='_blank' href='https://www.linkedin.com/in/olivier-sand-243570b5/'>Olivier Sand</a><sup>i</sup>,
<a target='_blank' href='http://www.ibens.ens.fr/spip.php?article26&lang=en'>Denis Thieffry</a><sup>ct</sup>,
and <a target='_blank' href='https://ibios.dkfz.de/tbi/index.php/about/eilslabs-people/128-carl-herrmann'>Carl Herrmann</a><sup>ct</sup>,
-->
</center>
<p><b>References</b>
<ol>
<li>Thomas-Chollier, M., Herrmann, C., Defrance, M., Sand,
  O., Thieffry, D. and van Helden, J. (2011). RSAT peak-motifs: motif
  analysis in full-size ChIP-seq datasets Nucleic Acids Research
  doi:10.1093/nar/gkr1104, 9.  [<a target='_blank'
  href='http://nar.oxfordjournals.org/content/early/2011/12/08/nar.gkr1104.full?keytype=ref&ijkey=zOvloLjtKzL73F8'>Open
  access</a>]
</li>
<li>
    Thomas-Chollier M, Darbo E, Herrmann C, Defrance M, Thieffry
    D, van Helden J. (2012). A complete workflow for the analysis
    of full-size ChIP-seq (and similar) data sets using
    peak-motifs. Nat Protoc 7(8): 1551-1568.  
    [<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/22836136'>PMID 22836136</a>]
</li>
</ol>
</p>

<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'105\')\" id=\"heading105\">

<font color='#0D73A7'>Information on the methods used in peak-motifs</font> </div>
 <div id=\"menu105\" class=\"menu_collapsible\">

The idea behind <i>peak-motifs</i> is that we detect <b>exceptional
words</b> based on <b>distinct and complementary criteria</b>:

<ul>

  <li> <b>global over-representation (oligo-analysis and
  dyad-analysis)</b>: a word/dyad is more frequent than expected from
  the background model. The over-repressentation is tested with a
  right-tailed binomial significance test.</li>

  <li> <b>positional bias (position-analysis)</b>: a word has a
  heterogeneous of occurrences in the test sequences, i.e. there are
  regions with higher frequency and other regions with lower
  frequencies than the average of the same word observed over the
  whole width of the sequences. Positional bias is tested with a
  chi-squared tests. </li>

  <li><b>local over-representation (local-word-analysis)</b>: the same
  test as for oligo-analysis (significance of the right tail of the
  binomial distribution) applied successfully to positional windows
  defined over the test set.</li>

</ul>

<br/>

For position-analysis and local-word-analysis, the <b>sequences</b>
are supposed to be <b>aligned</b> over some reference. For peaks, the
reference is the summit (or <b>center</b>) of each sequence.

</div></p>
end_part_1


&ListParameters() if ($ENV{rsat_echo} >=2);


## demo description
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

print $query->start_multipart_form(-action=>"peak-motifs.cgi", -onreset=>"resetHandler()");

################# Peak sequences
 &Panel1();

################# Restriction of data
 &Panel2();

################# Motif discovery parameters
 &Panel3();

################# Database comparison
 &Panel4();

################# sites and visualization
 &Panel5();

################# general parameters
 &Panel6();

################################################################
## Send results by email only
print "<p>\n";
#&SelectOutput('email', email_only=>1);
&SelectOutput('email');
print "<i>Note: email output is preferred for very large datasets or many comparisons with motifs collections</i>\n";

################################################################
## Action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;

################################################################
## Data for the demo single-set analysis
my $descr1 = "<H4>Comment on the demonstration example 1 :</H4>\n";
$descr1 .= "<blockquote class ='demo'>";
$descr1 .= "In this demonstration, we apply time- and memory-efficient
motif discovery algorithms to discover over-represented motifs in a
set of 1000 peak regions bound by the mouse transcription factor Oct4
(Chen et al., 2008)</p>\n";
$descr1 .= "</blockquote>";

$demo_url = $ENV{rsat_www}."/demo_files/peak-motifs_demo.fa";

print '<script>
function setDemo1(demo_url){
    $("#reset").trigger("click");
    $("#db_choice").val("' . $default{compare_motif_database} . '").change();
    descr = "<H4>Comment on the demonstration example 1 :</H4>\n";
    descr = descr + "<blockquote class =\'demo\'>";
    descr = descr + "In this demonstration, we apply time- and memory-efficient \
    motif discovery algorithms to discover over-represented motifs in a \
    set of 1000 peak regions bound by the mouse transcription factor Oct4 \
    (Chen et al., 2008)</p>\n</blockquote>";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;
    sequence_url1.value = demo_url;
    $("#title").val("Oct4 Chen2008 sites from Jaspar");
    max_seq_len.value = "";
    top_sequences.value = "";
    $("#visualize_galaxy").prop("checked", true);
    
}
function resetHandler(){
    $("#db_choice").val("").change();
}

$(function(){
    $("#db_choice").val("' . $default{compare_motif_database} . '").change();
});
</script>';

print "<TD><b>";
print '<button type="button" onclick="setDemo1('. "'$demo_url'".')">DEMO single</button>';
print "</B></TD>\n";

################################################################
## Data for the demo of differential analysis (test vs control)
$demo_url = $ENV{rsat_www}."/demo_files/peak-motifs_GSM559652_heart_p300_1000peaks.fa";
$ctrl_url = $ENV{rsat_www}."/demo_files/peak-motifs_GSM348066_limb_p300_1000peaks.fa";


print '<script>
function setDemo2(demo_url,ctrl_url){
    $("#reset").trigger("click");
    descr = "<H4>Comment on the demonstration example 2 :</H4>\n";
    descr = descr + "<blockquote class =\'demo\'>";
    descr = descr + "In this demonstration, we run a differential analysis \
    (test vs control) to discover the motifs that are over-represented in \
    one tissue (heart) compared to another tissue (limb), for a same \
    transcription factor (p300) (Blow et al, 2010)</p>\n</blockquote>";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;
    sequence_url1.value = demo_url;
    sequence_url2.value = ctrl_url;
    $("#title").val("p300 heart versus limb Blow2010");
    max_seq_len.value = "";
    $("#position-analysis").prop("checked",false);
    $("#oligo-analysis").prop("checked",true);
    $("#oligo_length6").prop("checked",true);
    $("#oligo_length7").prop("checked",true);
    $("#oligo_length8").prop("checked",false);
    $("#nmotifs").val(5);
    top_sequences.value = "";
    $("#visualize_galaxy").prop("checked", true);

};
</script>';

print "<TD><b>";
print '<button type="button" onclick="setDemo2('. "'$demo_url'".','."'$ctrl_url'".')">DEMO test vs ctrl</button>';
print "</B></TD>\n";

##print "<td><b><a href='tutorials/tut_peak_motif.html'>[TUTORIAL]</a></B></TD>\n";
print "<td><b><a class='iframe' href='help.peak-motifs.html'>[MANUAL]</a></B></TD>\n";
print "<td><b><a class='iframe' href='tutorials/tut_peak-motifs.html'>[TUTORIAL]</a></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



#################### Internal functions #######################

sub Panel1 {

  print "<fieldset>\n<legend><b><a class='iframe' href='help.formats.html'>Peak Sequences </a></b></legend>\n";
  print "<table>
  <tr><td colspan='2' style='text-align:center;'>";
  print "<b>Title</b> <font color='red'>(mandatory)</font>\n";
  print $query->textfield(-name=>'title', -id=>'title',
			  -default=>$default{title},
			  -size=>25);

 # print "</div>\n";

  print "</td>
  </tr>

  <tr><td style='padding-right:15px;border-right:1px solid #2D282E;'>";

  print "<span title=\"Provide here your peak sequences.This is the only mandatory input of the whole form\">";
  &MultiSequenceChoice("Peak sequences <font color='red'>(mandatory)</font>",1);
  print "</span>";
  print "<p/>\n";
  print "</td><td style='padding-left:15px;'>";
  print "<p><b>Optional:</b> <i>control dataset for differential analysis (test vs control)</i></p>\n";

  print "<p/>\n";
  print "<span title=\"Provide a second peak sequences set ONLY if you perform a differential analysis (mutant vs wild type, ...)\">";
  &MultiSequenceChoice("Control sequences",2);
  print "</span>";
  print "<p/>\n";

  print "</td></tr></table>";
  print '<b><font style="font-size:80%"><a class="iframe" href="tutorials/tut_peak-motifs.html#seq" target="_blank"> (I only have coordinates in a BED file, how to get sequences ?)</a></font></b></br>';
  print "</fieldset><p/>";
}

##########################################
sub Panel2 {
  print "<p class=\"clear\"></p>\n";

  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\">
     <span title=\"Panel for selecting only top peaks or cut sequences.\"><b>Reduce peak sequences</b></span></div>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible\">\n";
  print "";


  print "<p/><fieldset>\n";
  print "<legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Restrict the test dataset  </a></b></legend>\n";

  print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a class='iframe' href='help.peak-motifs.html#thresholds'>Number of top sequences to retain </a>&nbsp;</B>\n";
  print  $query->textfield(-name=>'top_sequences',-id=>'top_sequences',
			   -default=>$default{top_sequences},
			   -size=>3);

  print "<br/>";

  print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a class='iframe' href='help.peak-motifs.html#thresholds'>Cut peak sequences:</a> +/- &nbsp;</B>\n";
  print  $query->textfield(-name=>'max_seq_len',-id=>'max_seq_len',
			   -default=>$default{max_seq_len},
			   -size=>3);
  print "&nbsp;&nbsp;&nbsp;&nbsp;<b>bp on each side of peak centers</b>\n";


  print "</fieldset><p/>";
  print '</div><p class="clear"></p>';
}

##########################################
sub Panel3 {
print "
<div class=\"menu_heading_closed\" onclick=\"toggleMenu('102')\" id=\"heading102\">
 <span title=\"Select the algorithms and background here.\"><b>Motif discovery parameters</b></span> </div>\n";

print '
<div id="menu102" class="menu_collapsible">';

print "<p/><fieldset>
<legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Discover motifs </a></b></legend>";



################################################################
## Single words
print "<p> <b>Oligonucleotides (k-mers)</b>";

print "<ul>\n";

### oligo analysis
print "<br>";
print $query->checkbox(-name=>'oligo-analysis',-id=>'oligo-analysis',
		       -checked=>$default{"oligo-analysis"},
		       -label=>'');
print "&nbsp;<b>Discover over-represented words</b> <a class='iframe' href='help.oligo-analysis.html'> [oligo-analysis]</a>\n";


### position-analysis
print "<br>";
print $query->checkbox(-name=>'position-analysis',-id=>'position-analysis',
		       -checked=>$default{"position-analysis"},
		       -label=>'');
print "&nbsp;<b>Discover words with a positional bias</b> <a class='iframe' href='help.position-analysis.html'>[position-analysis]</a>\n";



### local-word-analysis
print "<br>", $query->checkbox(-name=>'local-word-analysis',
		       -checked=>$default{'local-word-analysis'},
		       -label=>'');
print "&nbsp;<b>Discover words with local over-representation</b> <a class='iframe' href='help.local-word-analysis.html'>[local-word-analysis]</a>\n";
print "&nbsp;"x2;
print "<br><i>Note:position-analysis and local-word-analysis will not run if a control set is provided</i>";


## Word size
print "<p><b><a class='iframe' href='help.oligo-analysis.html#oligo_length'>Oligomer lengths</a>&nbsp;</b> for the three programs above\n";

print $query->checkbox(-name=>'oligo_length6',-id=>'oligo_length6',
		       -checked=>$default{oligo_length6},
		       -label=>'6');
print "&nbsp;"x2;
print $query->checkbox(-name=>'oligo_length7',-id=>'oligo_length7',
		       -checked=>$default{oligo_length7},
		       -label=>'7');
print "&nbsp;"x2;
print $query->checkbox(-name=>'oligo_length8',-id=>'oligo_length8',
		       -checked=>$default{oligo_length8},
		       -label=>'8');
print "&nbsp;"x6;
print $query->checkbox(-name=>'merge_lengths',-id=>'merge_lengths',
		       -checked=>$default{merge_lengths},
		       -label=>'merge lengths for assembly');
print "<br><i>Note: motifs can be larger than word sizes (words are used as seed for building matrices)</i>";


## Markov  order (for oligo-analysis in single strand mode)
print "<br><p><b><a class='iframe' href='help.oligo-analysis.html'> Markov order (m) of the background model for oligo-analysis (k-mers)</a> </b><i>(only for single-dataset analysis, will be ignored if control set is provided)</i>\n";
$oligoPopup = "<br>";
$oligoPopup .=  "<select name='markov'>\n";
$oligoPopup .=  "<option value='auto'>automatic (adapted to sequence length)</option>\n";
$oligoPopup .=  "<option value='0'>m=0 (generally not ideal)</option>\n";
$oligoPopup .=  "<option value='1'>m=1 (more sensitive for small data sets, e.g. 100kb)</option>\n";
$oligoPopup .=  "<option value='-3'>m=k-3 (intermediate size data sets)</option>\n";
$oligoPopup .=  "<option value='-2'>m=k-2 (more stringent for large data sets e.g. > 1Mb)</option>\n";
$oligoPopup .=  "</select>";
print $oligoPopup;

print "</ul>\n";


################################################################
## dyad-analysis
print "<p> <b>Spaced word pairs (dyads)</b>";
print "<ul>\n";
print $query->checkbox(-name=>'dyad-analysis',-id=>'dyad-analysis',
			       -checked=>$default{"dyad-analysis"},
			       -label=>'');
print "&nbsp;<b>Discover over-represented spaced word pairs </b><a class='iframe' href='help.dyad-analysis.html'>[dyad-analysis] </a>\n";
print "</ul>\n";


### dyad sizes and spacer
#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.dyad-analysis.html#oligo_size'>Dyad length and spacer</a>&nbsp;</B>\n";
#	$oligoPopup = "";
#    $oligoPopup .=  "<select NAME='dyad-option'>\n";
#	$oligoPopup .=  "<option  value='4'>4 {0,20} 4</option>\n";
#	$oligoPopup .=  "<option  selected value='3'>3 {0,20} 3</option>\n";
#    $oligoPopup .=  "</select>";
#    print $oligoPopup;
#print "<br/>";

## Number of motifs per algorithm
print "<b>Number of motifs per algorithm</b>\n";
print $query->popup_menu(-name=>'nmotifs',-id=>'nmotifs',
			 -Values=>[1,2,3,4,5,6,7,8,9,10],
			 -default=>$default{nmotifs});

### 2str or 1str
print "<br/><b>Search on </b> ";
my $strandPopup =  "<select NAME='strand'>\n";
$strandPopup .=  "<option  selected value='-2str'>both strands</option>\n";
$strandPopup .=  "<option value='-1str'>single strand</option>\n";
$strandPopup .=  "</select>";
print $strandPopup;

### local-word-analysis => too slow to be put on the website until the program is optimized
#
#print $query->checkbox(-name=>'local-word-analysis_dyads',
#		       -checked=>$default{'local-word-analysis_dyads'},
#		       -label=>'');
#print "&nbsp;<b>Discover spaced word pairs with local over-representation</b> <a href='help.local-word-analysis.html'>[local-word-analysis]</a>\n";
#print "<br/>";

#
### position-analysis => still to be implemented
#
#print $query->checkbox(-name=>'position-analysis_dyads',
#		       -checked=>$default{"position-analysis_dyads"},
#		       -label=>'');
#print "&nbsp;<b>Discover words with a positional biais</b> <a href='help.local-word-analysis.html'>[position-analysis]</a>\n";
#print "<br/>";


### threshold (common to all programs)
#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.oligo-analysis.html#thresholds'>Lower threshold on significance</a>&nbsp;</B>\n";
#print  $query->textfield(-name=>'lth_occ_sig',
#							      -default=>$default{lth_occ_sig},
#							      -size=>3);
#print "<br/>";
print "</fieldset><p/>";

print '</div>';

}

################################################################
## Comparisons with motif databases
sub Panel4 {
  print '
<br/>
<div>';
  print "
<div class=\"menu_heading_closed\" onclick=\"toggleMenu('103')\" id=\"heading103\">

 <span title=\"Discovered motifs will be compared to known motifs. Select here the known motifs collections and/or provide your own.\"><b>Compare discovered motifs with databases (e.g. against Jaspar) or custom reference motifs</b></span> </div>\n
<div id=\"menu103\" class=\"menu_collapsible\">";

  ## Tasks
  print "<fieldset><legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Compare motifs </a></b></legend>";

  ### compare motifs
#  print $query->checkbox(-name=>'compare_motif_db',
#			 -checked=>$default{compare_motif_db},
#			 -label=>'');
#  print "&nbsp;<b>Compare discovered motifs with known motifs from databases</b> <a href=''>[compare-matrices]</a>\n";

  print "<p/> ";
#  print "<a href=''><b>Choose below the motif database(s):</b></a><br/>";
  print "<b>Compare discovered motifs with known motifs from databases</b><br/>";

&MotifSelection("mode" => "checkbox");
  ## Display supported motif databases
  # &DisplayMatrixDBchoice("mode"=>"checkbox");

  print "<p/> ";
  print "<b>Add your own motif database:</b><br/>";
  print  $query->textfield(-name=>'custom_motif_db_name',
			   -default=>$default{custom_motif_db_name},
			   -size=>20);
  print $query->filefield(-name=>'custom_motif_db',
			  -size=>10);
#  print "<br>Matrices should be in <b>Transfac format</b> (other formats can be converted with <a href='convert-matrix_form.cgi'><i>convert-matrix</i></a>).";

 ## Reference motifs
  print"</p>";
  print "<b>Add known reference motifs for this experiment:</b><br/>";
  print $query->filefield(-name=>'ref_motif',
			  -size=>10);
  print "<br>Database and reference motifs (matrices) should be in <b>Transfac format</b>";
  print "<br>(other formats can be converted with <a href='convert-matrix_form.cgi'><i>convert-matrix</i></a>).";


  print "</fieldset><p/>";

  print '</div><p class="clear"></p>';
}

##########################################

sub Panel5  {
  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu('104')\" id=\"heading104\">";
  print "<span title=\"Input peaks are scanned with the discovered motifs to obtain their exact position.\n";
  print "These putative binding sites can be visualized on genome browsers (Ensembl, UCSC genome browser,...)\">\n";
  print "<b>Locate motifs and export predicted sites as custom UCSC tracks</b></span>";

  print "</div>";
  print "<div id=\"menu104\" class=\"menu_collapsible\">";


  ### matrix-scan
  print "<fieldset><legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Locate motifs </a></b></legend>";
  print $query->checkbox(-name=>'matrix-scan-quick',
			 -checked=>$default{'matrix-scan-quick'},
			 -label=>'');
  print "&nbsp;<b>Search putative binding sites in the peak sequences</b> <a class='iframe' href='help.matrix-scan.html'>[matrix-scan]</a>\n";
  print "<br/>";

  ## Markov order for scanning (site prediction + motif enrichment)
  print "<b>Markov order</b> (m) of the background model for sequence scanning (site prediction + motif enrichment)\n";
  $scanMarkovPopup = "<br>";
  $scanMarkovPopup .=  "<select name='scan_markov'>\n";
  $scanMarkovPopup .=  "<option value='0'>m=0 (generally not ideal)</option>\n";
  $scanMarkovPopup .=  "<option selected value='1'>m=1</option>\n";
  $scanMarkovPopup .=  "<option value='2'>m=2 (slower)</option>\n";
  $scanMarkovPopup .=  "<option value='3'>m=3 (may be very slow)</option>\n";
  $scanMarkovPopup .=  "</select>";
  print $scanMarkovPopup;


  print "</fieldset><p/>";


  ################################################################
  ## Visualize UCSC custom track
  print "<fieldset><legend><b><a class='iframe' href='help.peak-motifs.html#tasks'>Visualize peaks and sites in genome browser </a></b></legend>";

  print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_none' value='none' $checked{'none'}>","<b>No</b>");

  print "<br/>";
  print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_getfasta' value='getfasta' $checked{'getfasta'}>",
	 "Peak coordinates specified in fasta headers in <a target='_blank' href='http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html'><b>bedtools getfasta</b></a> format (also for <b>retrieve-seq-bed</b> output).",
	 "<br>","&nbsp;"x7,"Fasta headers should be in the form: <tt>>3:81458-81806(.)</tt>");

  print "<br/>";
  print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_galaxy' value='galaxy' $checked{'galaxy'}>",
	 "Peak coordinates specified in fasta headers of the test sequence file (<a target='_blank' href='https://usegalaxy.org/'><b>Galaxy</b></a> format)",
	 "<br>","&nbsp;"x7,"Fasta headers should be in the form: <tt>>mm9_chr1_3473041_3473370_+</tt>");

  print "<br/>";
  print ("<INPUT TYPE='radio' NAME='visualize' id='visualize_bed_coord' value='bed_coord' $checked{'bed_coord'}>","Peak coordinates provided as a <b>custom BED file.</b>");
  print "&nbsp;"x7, "<br>The 4th column of the BED file (feature name) must correspond to the fasta headers of sequences</i><br/>";

  print "&nbsp;"x7, $query->filefield(-name=>'bed_file',
				      -size=>10);

  ### assembly
  print "&nbsp;&nbsp;&nbsp;&nbsp;<b>Assembly version (UCSC)&nbsp;</B>\n";
  print  $query->textfield(-name=>'assembly',
							      -default=>$default{assembly},
							      -size=>10);

  print "<br/></fieldset><p/>";

  ## Close divisions
  print "</div>\n";
  print "<p class='clear'></p>\n";

 }

##########################################

sub Panel6  {
  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu('106')\" id=\"heading106\">";
  print "<span title=\"Reporting options.\n";
  print "Display preferences and plotting options\">\n";
  print "<b>Reporting options</b></span>";

  print "</div>";
  print "<div id=\"menu106\" class=\"menu_collapsible\">";

  ### Plotting options
  print "<fieldset><legend><b>Plotting options</b></legend>";

  ################################################################
  ## Reference location for position-analysis and to scan sequences
  print "<br> <b>Reference position for position-analysis and sequence scanning</b>";

  ## Origin for calculating positions and scanning sequences
  print "&nbsp;"x4,  "<A class='iframe' HREF='help.position-analysis.html#origin'><B>Origin</B></A>\n";
  print $query->popup_menu(-name=>'origin',
			   -Values=>['start',
				     'center',
				     'end'],
			   -default=>$default{origin});
  
  ## Offset for calculating positions
  print "&nbsp;"x4,  "<A class='iframe' HREF='help.position-analysis.html#offset'><B>Offset</B></A>\n";
  print $query->textfield(-name=>'offset',
			  -default=>$default{offset},
			  -size=>8);


  
  ## Use R to generate XY plots.
  if ($ENV{R_supported}) {
      ## This option is displayed only if R is supported on the server. 
      print "<br>";
      print $query->checkbox(-name=>'r_plot',
			     -checked=>$default{"r_plot"},
			     -label=>'');
      print "&nbsp;<b>Use R to generate plots</b> (only works for servers with R installed).\n";
      
  }

  print "</fieldset><p/>";

  ## Close divisions
  print "</div>\n";
  print "<p class='clear'></p>\n";
 }
