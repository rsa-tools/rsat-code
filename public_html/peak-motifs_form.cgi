#!/usr/bin/perl
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
$default{demo_descr1} = "";
$default{demo_descr2} = "";

################################################################
### default values for filling the form
$default{lth_occ_sig}=0;
$default{uth_pval} = "1e-4";
$default{assembly} = "";
#$default{oligo_length6}="checked";
$default{oligo_length7}="checked";
$default{"oligo-analysis"}="checked";
$default{"dyad-analysis"}="";
$default{"position-analysis"}="checked";
$default{'local-word-analysis'}="";
$default{'local-word-analysis_dyads'} ="";
$default{"position-analysis_dyads"} ="";
$default{matrix-scan-quick}="checked";
$default{compare_motif_db}="checked";
$default{title}="title for this analysis";
$default{max_seq_len}="";
$default{top_sequences}="";


$default{visualize}="none";
$checked{$default{visualize}} = "CHECKED";	


## motif database
$default{compare_motif_database}="jaspar_core_vertebrates";
$default{perso_motif_name}="title for this collection";


### replace defaults by parameters from the cgi call, if defined
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



################################################################
### print the form ###

################################################################
### header
&RSA_header("peak-motifs", "form");
print "<CENTER>";
print "Pipeline for discovering motifs in massive ChIP-seq peak sequences.<P>\n";
print "<br>Conception<sup>c</sup>, implementation<sup>i</sup> and testing<sup>t</sup>: ";
print "<a target='_blank' href='http://www.bigre.ulb.ac.be/Users/jvanheld/'>Jacques van Helden</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/morgane/'>Morgane Thomas-Chollier</a><sup>cit</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/defrance/'>Matthieu Defrance</a><sup>ci</sup>\n";
print ", <a target='_blank' href='http://www.bigre.ulb.ac.be/Users/oly/'>Olivier Sand</a><sup>i</sup>\n";
print ", <a target='_blank' href='http://www.ibens.ens.fr/spip.php?article26&lang=en'>Denis Thieffry</a><sup>ct</sup>\n";
print "and <a target='_blank' href='http://biologie.univ-mrs.fr/view-data.php?id=202'>Carl Herrmann</a><sup>ct</sup>\n";
print "</CENTER>";
print "</BLOCKQUOTE>\n";
print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'105\')\" id=\"heading105\"><font color='#0D73A7'>Information on the methods used in peak-motifs</font> </div>\n";
 print "<div id=\"menu105\" class=\"menu_collapsible\">\n";
print "<BLOCKQUOTE>\n";
print "The idea behind <i>peak-motifs</i> is that we detect <b>exceptional words</b> based on <b>distinct and complementary criteria</b>:
<ul>
<li> <b>global over-representation (oligo-analysis and dyad-analysis)</b>: a word/dyad is more frequent than expected from the background model. The over-repressentation is tested with a right-tailed binomial significance test.</li>

<li> <b>positional bias (position-analysis)</b>: a word has a heterogeneous of occurrences in the test sequences, i.e. there are regions with higher frequency and other regions with lower frequencies than the average of the same word observed over the whole width of the sequences. Positional bias is tested with a chi-squared tests. 
</li>
<li><b>local over-representation (local-word-analysis)</b>: the same test as for oligo-analysis (significance of the right tail of the binomial distribution) applied successfully to positional windows defined over the test set.
</li>
<br/>
For position-analysis and local-word-analysis, the <b>sequences</b> are supposed to be <b>aligned</b> over some reference. For peaks, the reference is the summit (or <b>center</b>) of each sequence. ";
print "</BLOCKQUOTE>\n";
print "</div></p>\n";




&ListParameters() if ($ENV{rsat_echo} >=2);


## demo description
print $default{demo_descr1};

print $query->start_multipart_form(-action=>"peak-motifs.cgi");

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
 
################################################################
### send results by email only
print "<p>\n";
#&SelectOutput('email', email_only=>1);
&SelectOutput();
print "<i>Note: email output is preferred for very large datasets or many comparisons with motifs collections</i>\n";

################################################################
### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

################################################################
### data for the demo single

my $descr1 = "<H4>Comment on the demonstration example 1 :</H4>\n";
$descr1 .= "<blockquote class ='demo'>";

$descr1 .= "In this demonstration, we apply time- and memory-efficient
motif discovery algorithms to discover over-represented motifs in a
set of 1000 peak regions bound by the mouse transcription factor Oct4
(Chen et al., 2008)</p>\n";

#$descr1 .= "Discovered motifs are compared to JASPAR vertebrate
#motifs, and sequences are scanned to predict binding sites.</p>\n";

$descr1 .= "</blockquote>";

#print $query->start_multipart_form(-action=>"peak-motifs_form.cgi");
#$demo_seq=`cat demo_files/chip-seq-analysis_demo.fa`;
print $query->start_multipart_form(-action=>"peak-motifs_form.cgi");
$demo_seq=`cat demo_files/peak-motifs_demo.fa`;
print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'sequence1',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format1',-default=>'fasta');
print $query->hidden(-name=>'title',-default=>'Oct4 Chen2008 sites from Jaspar');
print $query->hidden(-name=>'max_seq_len',-default=>'');
print $query->hidden(-name=>'top_sequences',-default=>'');
print $query->hidden(-name=>'visualize',-default=>"galaxy");
#print $query->hidden(-name=>'user_email',-default=>'nobody@nowhere');
print $query->submit(-label=>"DEMO single");
print "</B></TD>\n";
print $query->end_form;

################################################################
### data for the demo differential (test vs control)

my $descr2 = "<H4>Comment on the demonstration example 2 :</H4>\n";
$descr2 .= "<blockquote class ='demo'>";

$descr2 .= "In this demonstration, we run a differential analysis (test vs control)
to discover the motifs that are over-represented in one tissue (heart) compared to another tissue (limb), for a same transcription factor (p300) (Blow et al, 2010)</p>\n";

$descr2 .= "</blockquote>";


print $query->start_multipart_form(-action=>"peak-motifs_form.cgi");
$demo_seq=`cat demo_files/peak-motifs_GSM559652_heart_p300_peaks.fa`;
$ctrl_seq=`cat demo_files/peak-motifs_GSM348066_limb_p300_peaks.fa`;
print "<TD><b>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr2);
print $query->hidden(-name=>'sequence1',-default=>$demo_seq);
print $query->hidden(-name=>'sequence2',-default=>$ctrl_seq);
print $query->hidden(-name=>'sequence_format1',-default=>'fasta');
print $query->hidden(-name=>'sequence_format2',-default=>'fasta');
print $query->hidden(-name=>'title',-default=>'p300 heart versus limb Blow2010');
print $query->hidden(-name=>'max_seq_len',-default=>'');
print $query->hidden(-name=>'top_sequences',-default=>'');
print $query->hidden(-name=>'visualize',-default=>"galaxy");
#print $query->hidden(-name=>'user_email',-default=>'nobody@nowhere');
print $query->submit(-label=>"DEMO test vs ctrl");
print "</B></TD>\n";
print $query->end_form;


##print "<td><b><a href='tutorials/tut_peak_motif.html'>[TUTORIAL]</a></B></TD>\n";
print "<td><b><a href='help.peak-motifs.html'>[MANUAL]</a></B></TD>\n";
print "<td><b><a href='tutorials/tut_peak-motifs.html'>[TUTORIAL]</a></B></TD>\n";
print "<TD><b><a href='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</a></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



#################### Internal functions #######################

sub Panel1 {

  print "<fieldset>\n<legend><b><a href='help.formats.html'>Peak Sequences </a></b></legend>\n";
 
  print "<table>
  <tr><td colspan='2' style='text-align:center;'>";
   print "<b>Title</B>\n";
  print $query->textfield(-name=>'title',
			  -default=>$default{title},
			  -size=>25);

 # print "</div>\n";
  
  print "</td>
  </tr>
  
  <tr><td style='padding-right:15px;border-right:1px solid #2D282E;'>";
  
  print "<span title=\"Provide here your peak sequences.This is the only mandatory input of the whole form\">";
   &MultiSequenceChoice("Peak sequences",1);
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
  print '<b><font style="font-size:80%"><a href="tutorials/tut_peak-motifs.html#seq" target="_blank"> (I only have coordinates in a BED file, how to get sequences ?)</a></font></b></br>';
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
  print "<legend><b><a href='help.peak-motifs.html#tasks'>Restrict the test dataset  </a></b></legend>\n";

  print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.peak-motifs.html#thresholds'>Number of top sequences to retain </a>&nbsp;</B>\n";
  print  $query->textfield(-name=>'top_sequences',
			   -default=>$default{top_sequences},
			   -size=>3);

  print "<br/>";

  print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.peak-motifs.html#thresholds'>Cut peak sequences:</a> +/- &nbsp;</B>\n";
  print  $query->textfield(-name=>'max_seq_len',
			   -default=>$default{max_seq_len},
			   -size=>3);
  print "&nbsp;&nbsp;&nbsp;&nbsp;<b>bp around the center of each peak </b>\n";


  print "</fieldset><p/>";
  print '</div>
</div>
<p class="clear"></p>';
}

##########################################
sub Panel3 {
print "
<div class=\"menu_heading_closed\" onclick=\"toggleMenu('102')\" id=\"heading102\">
 <span title=\"Select the algorithms and background here.\"><b>Change motif discovery parameters</b></span> </div>\n";

print '
<div id="menu102" class="menu_collapsible">';

print "<p/><fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Discover motifs </a></b></legend>";


print "<br/> <b>Continuous words</b> <br/>";

#
### oligo analysis
#
print $query->checkbox(-name=>'oligo-analysis',
		       -checked=>$default{"oligo-analysis"},
		       -label=>'');  
print "&nbsp;<b>Discover over-represented words</b> <a href='help.oligo-analysis.html'> [oligo-analysis]</a>\n";
print "<br/>";


### dyad sizes and spacer
#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.dyad-analysis.html#oligo_size'>Dyad length and spacer</a>&nbsp;</B>\n";
#	$oligoPopup = "";
#    $oligoPopup .=  "<SELECT NAME='dyad-option'>\n";
#	$oligoPopup .=  "<OPTION  VALUE='4'>4 {0,20} 4</option>\n";
#	$oligoPopup .=  "<OPTION  SELECTED VALUE='3'>3 {0,20} 3</option>\n";
#    $oligoPopup .=  "</SELECT>";
#    print $oligoPopup;
#print "<br/>";

#
### local-word-analysis
#
print $query->checkbox(-name=>'local-word-analysis',
		       -checked=>$default{'local-word-analysis'},
		       -label=>'');  
print "&nbsp;<b>Discover words with local over-representation</b> <a href='help.local-word-analysis.html'>[local-word-analysis]</a>\n";
print "<br/>";

#
### position-analysis
#
print $query->checkbox(-name=>'position-analysis',
		       -checked=>$default{"position-analysis"},
		       -label=>'');  
print "&nbsp;<b>Discover words with a positional biais</b> <a href='help.local-word-analysis.html'>[position-analysis]</a>\n";


## Word size
print "<p><b><a href='help.oligo-analysis.html#oligo_length'>Oligomer length</a>&nbsp;</b> for the three programs above\n";

print $query->checkbox(-name=>'oligo_length6',
		       -checked=>$default{oligo_length6},
		       -label=>'6'); 
print "&nbsp;"x2;
print $query->checkbox(-name=>'oligo_length7',
		       -checked=>$default{oligo_length7},
		       -label=>'7'); 
print "&nbsp;"x2;
print "<br><i>Note: motifs can be larger than word sizes (words are used as seed for building matrices)</i>";
print "</p>";

################################################################
## dyad-analysis
print "<p> <b>Spaced words pairs</b>";
print "<br>", $query->checkbox(-name=>'dyad-analysis',
			       -checked=>$default{"dyad-analysis"},
			       -label=>'');
print "&nbsp;<b>Discover over-represented spaced word pairs </b><a href='help.dyad-analysis.html'>[dyad-analysis] </a>\n";
print "</p>";

### 2str or 1str

print "<br/> <b>Search on </b> ";
my $strandPopup =  "<SELECT NAME='strand'>\n";
$strandPopup .=  "<OPTION  SELECTED VALUE='-2str'>both strands</option>\n";
$strandPopup .=  "<OPTION VALUE='-1str'>single strand</option>\n";
$strandPopup .=  "</SELECT>";
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

### markov  order (common to all programs)
print "<p><b><a href='help.oligo-analysis.html'> Markov order of the background model</a> </b><i>(only for single-dataset analysis, will be ignored if control set is provided)</i>\n";
$oligoPopup = "<br>";
$oligoPopup .=  "<SELECT NAME='markov'>\n";
$oligoPopup .=  "<OPTION  SELECTED VALUE='0'>0 (generally not ideal)</option>\n";
$oligoPopup .=  "<OPTION  SELECTED VALUE='1'>1 (more sensitive for small data sets, e.g. 100kb)</option>\n";
$oligoPopup .=  "<OPTION  SELECTED VALUE='-3'>oligo length -3 (intermediate size data sets)</option>\n";
$oligoPopup .=  "<OPTION VALUE='-2'>oligo length -2 (more stringent for large data sets e.g. > 1Mb)</option>\n";
$oligoPopup .=  "</SELECT>";
print $oligoPopup;
print "</p>";

### threshold (common to all programs)
#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.oligo-analysis.html#thresholds'>Lower threshold on significance</a>&nbsp;</B>\n";
#print  $query->textfield(-name=>'lth_occ_sig',
#							      -default=>$default{lth_occ_sig},
#							      -size=>3);
#print "<br/>";
print "</fieldset><p/>";

print '	
</div></div>';

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
  print "<fieldset><legend><b><a href='help.peak-motifs.html#tasks'>Compare motifs </a></b></legend>";

  ### compare motifs
#  print $query->checkbox(-name=>'compare_motif_db',
#			 -checked=>$default{compare_motif_db},
#			 -label=>'');
#  print "&nbsp;<b>Compare discovered motifs with known motifs from databases</b> <a href=''>[compare-matrices]</a>\n";

  print "<p/> ";
#  print "<a href=''><b>Choose below the motif database(s):</b></a><br/>";
  print "<a href=''><b>Compare discovered motifs with known motifs from databases</b></a><br/>";

  ## load the various databases that can be compared against
  &MatrixDBcheckBox("choice_mode"=>"checkbox");

  print "<p/> ";
  print "<a href=''><b>Add your own motif database:</b></a><br/>";
  print  $query->textfield(-name=>'perso_motif_name',
			   -default=>$default{perso_motif_name},
			   -size=>20);
  print $query->filefield(-name=>'perso_motif',
			-size=>10);
  print "Matrices should be in <b>Transfac format</b> (other formats can be converted with <a href='convert-matrix_form.cgi'><i>convert-matrix</i></a>).";
 ##
 print"</p>";
   print "<a href=''><b>Add known reference motifs for this experiment:</b></a><br/>";
  print $query->filefield(-name=>'ref_motif',
			-size=>10);
  print "Matrices should be in <b>Transfac format</b> (other formats can be converted with <a href='convert-matrix_form.cgi'><i>convert-matrix</i></a>).";
 
 
  print "</fieldset><p/>";

  print '
</div>
</div>
<p class="clear"></p>';
}

##########################################
 sub Panel5  {
print "
<div class=\"menu_heading_closed\" onclick=\"toggleMenu('104')\" id=\"heading104\"> 
<span title=\"Input peaks are scanned with the discovered motifs to obtain their exact position. These putative binding sites can be visualized on genome browsers (Ensembl, UCSC genome browser,...)\"><b>Locate motifs and export as UCSC custom track</b></span> </div>\n

<div id=\"menu104\" class=\"menu_collapsible\">";



#
### matrix-scan
#
print "<fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Locate motifs </a></b></legend>";
print $query->checkbox(-name=>'matrix-scan-quick',
		       -checked=>$default{matrix-scan-quick},
		       -label=>'');  
print "&nbsp;<b>Search putative binding sites in the peak sequences</b> <a href='help.matrix-scan.html'>[matrix-scan]</a>\n";

print "<br/>";

#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.matrix-scan.html'>Background model: Markov order</a>&nbsp;</B>\n";
#	$oligoPopup = "";
#    $oligoPopup .=  "<SELECT NAME='markov-scan'>\n";
#	$oligoPopup .=  "<OPTION  SELECTED VALUE='1'>0</option>\n";
#	$oligoPopup .=  "<OPTION  SELECTED VALUE='2'>1</option>\n";
#	$oligoPopup .=  "<OPTION  SELECTED VALUE='3'>2</option>\n";
#    $oligoPopup .=  "</SELECT>";
#   print $oligoPopup;
    
#print "<br/>";

### threshold (common to all programs)
#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.matrix-scan.html'>Upper threshold on P-value</a>&nbsp;</B>\n";
#print  $query->textfield(-name=>'uth_pval',
#							      -default=>$default{uth_pval},
#							      -size=>4);



print "</fieldset><p/>";


################################################################
## UCSC custom track
##
print "<fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Visualize motifs in genome browser </a></b></legend>";

print ("<INPUT TYPE='radio' NAME='visualize' VALUE='none' $checked{'none'}>","<b>No</b>");
print "<br/>";

print ("<INPUT TYPE='radio' NAME='visualize' VALUE='galaxy' $checked{'galaxy'}>",
       "<b>Yes; sequences fetched from <a href=''>Galaxy</a></b>",
       "(fasta headers should be in the form: <tt>>mm9_chr1_3473041_3473370_+ </tt>)");

print "<br/>";
print ("<INPUT TYPE='radio' NAME='visualize' VALUE='bed_coord' $checked{'bed_coord'}>","<b>Yes; use the following BED file.</b>");
print "<br/>";


print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.peak-motifs.html'>BED file with peak coordinates</a>&nbsp;</B>\n";
print $query->filefield(-name=>'bed_file',
			-size=>10);

### assembly
print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.peak-motifs.html'>Assembly version (UCSC)</a>&nbsp;</B>\n";
print  $query->textfield(-name=>'assembly',
							      -default=>$default{assembly},
							      -size=>10);
print "<br/>&nbsp;&nbsp;&nbsp;&nbsp;<i>The 4th column of the BED file (feature name) corresponds to the fasta headers of sequences</i>";

print "<br/>
</fieldset><p/>";

#########
print '	
</div>
</div>
<p class="clear"></p>';

 }
