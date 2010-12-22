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
$default{oligo_length6}="checked";
$default{oligo_length7}="checked";
$default{"oligo-analysis"}="checked";
$default{"dyad-analysis"}="";
$default{"position-analysis"}="checked";
$default{'local-word-analysis'}="";
$default{compare_motif_db}="checked";
$default{title}="title for this dataset";
$default{max_seq_len}="";
$default{top_sequences}="";

## motif database
$default{compare_motif_database}="jaspar_core_vertebrates";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}


################################################################
### print the form ###


################################################################
### header
&RSA_header("peak-motifs", "form");
print "<CENTER>";
print "Pipeline for discovering motifs in ChIP-seq peak sequences.<P>\n";
#print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

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
&SelectOutput('email', email_only=>1);

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
set of 1370 peak regions bound by the mouse transcription factor Otc4
(Chen et al., 2008)</p>\n";

$descr1 .= "Discovered motifs are compared to JASPAR vertebrate
motifs, and sequences are scanned to predict binding sites.</p>\n";

$descr1 .= "</blockquote>";

#print $query->start_multipart_form(-action=>"peak-motifs_form.cgi");
#$demo_seq=`cat demo_files/chip-seq-analysis_demo.fa`;
print $query->start_multipart_form(-action=>"peak-motifs_form.cgi");
$demo_seq=`cat demo_files/peak-motifs_demo.fa`;
print "<TD><B>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'sequence1',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format1',-default=>'fasta');
print $query->hidden(-name=>'title',-default=>'Oct4 Chen2008 sites from Jaspar');
print $query->hidden(-name=>'max_seq_len',-default=>'250');
print $query->hidden(-name=>'top_sequences',-default=>'5000');
print $query->submit(-label=>"DEMO 1");
print "</B></TD>\n";
print $query->end_form;



print "<TD><B><A HREF='help.peak-motifs.html'>[MANUAL]</A></B></TD>\n";
print "<TD><B><A HREF='http://www.bigre.ulb.ac.be/forums/' target='_top'>[ASK A QUESTION]</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



#################### Internal functions #######################

sub Panel1 {

  print "<fieldset>\n<legend><b><a href='help.formats.html'>Peak Sequences </a></b></legend>\n";
  print "<B>Title</B>\n";
  print $query->textfield(-name=>'title',
			  -default=>$default{title},
			  -size=>15);

  print "<p/>\n";

  &MultiSequenceChoice("Peak sequences",1);

  print '<b><font style="font-size:80%"><a href=""> (I only have coordinates in a BED file, how to get sequences ?)</a></font></b></br>';
  print "</fieldset><p/>";
}

##########################################
sub Panel2 {
  print "<p class=\"clear\"></p>\n";
  print "<div class=\"menu_heading_closed\" onclick=\"toggleMenu(\'101\')\" id=\"heading101\"><b>Reduce input peak sequences </b> </div>\n";
  print "<div id=\"menu101\" class=\"menu_collapsible\">\n";

  print "<p/><fieldset>\n";
  print "<legend><b><a href='help.peak-motifs.html#tasks'>Restrict the input dataset  </a></b></legend>\n";

  print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.peak-motifs.html#thresholds'>Number of top sequences to retain </A>&nbsp;</B>\n";
  print  $query->textfield(-name=>'top_sequences',
			   -default=>$default{top_sequences},
			   -size=>3);

  print "<br/>";

  print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.peak-motifs.html#thresholds'>Cut peak sequences:</A> +/- &nbsp;</B>\n";
  print  $query->textfield(-name=>'max_seq_len',
			   -default=>$default{max_seq_len},
			   -size=>3);
  print "&nbsp;&nbsp;&nbsp;&nbsp;<B>bp around the center of each peak </b>\n";


  print "</fieldset><p/>";
  print '</div>
</div>
<p class="clear"></p>';
}

##########################################
sub Panel3 {
print '
<div class="menu_heading_closed" onclick="toggleMenu(\'102\')" id="heading102"><b>Change motif discovery parameters </b> </div>
<div id="menu102" class="menu_collapsible">';

print "<p/><fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Discover motifs </a></b></legend>";

#
### oligo analysis
#
print $query->checkbox(-name=>'oligo-analysis',
		       -checked=>$default{"oligo-analysis"},
		       -label=>'');  
print "&nbsp;<b>Discover over-represented words</b> <a href='help.oligo-analysis.html'> [oligo-analysis]</a>\n";
print "<br/>";
	       
		       
#
### dyad-analysis
#   
print $query->checkbox(-name=>'dyad-analysis',
		       -checked=>$default{"dyad-analysis"},
		       -label=>''); 
print "&nbsp;<b>Discover over-represented spaced word pairs </b><a href='help.dyad-analysis.html'>[dyad-analysis] </a>\n";

print "<br/>";	
### dyad sizes and spacer
#print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.dyad-analysis.html#oligo_size'>Dyad length and spacer</A>&nbsp;</B>\n";
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
print "<br/>";


#print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.local-word-analysis.html#window_width'>Fixed window of width</A>&nbsp;</B>\n";
#	$oligoPopup = "";
#    $oligoPopup .=  "<SELECT NAME='orm_window'>\n";
#	$oligoPopup .=  "<OPTION VALUE='50'>50</option>\n";
#	$oligoPopup .=  "<OPTION SELECTED VALUE='20'>20</option>\n";
#	$oligoPopup .=  "<OPTION VALUE='10'>10</option>\n";
#    $oligoPopup .=  "</SELECT>";
#    print $oligoPopup;

print "<br/>";	
		       

### markov  order (common to all programs)
print "<br/> <b>Common options for above programs</b> <br/>";

print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-analysis.html#oligo_length'>Oligomer length</A>&nbsp;</B>\n";

print $query->checkbox(-name=>'oligo_length6',
		       -checked=>$default{oligo_length6},
		       -label=>'6'); 
print "&nbsp;"x2;
print $query->checkbox(-name=>'oligo_length7',
		       -checked=>$default{oligo_length7},
		       -label=>'7'); 
print "&nbsp;"x2;
print $query->checkbox(-name=>'oligo_length8',
		       -checked=>$default{oligo_length8},
		       -label=>'8'); 
print "&nbsp;"x2;

print "<br/>";

print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-analysis.html'>Background model: Markov order</A>&nbsp;</B>\n";
 	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='markov'>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='-2'>oligo length -2 </option>\n";
	$oligoPopup .=  "<OPTION VALUE='-3'>oligo length -3</option>\n";
    $oligoPopup .=  "</SELECT>";
    print $oligoPopup;
    
print "<br/>";

### threshold (common to all programs)
#print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-analysis.html#thresholds'>Lower threshold on significance</A>&nbsp;</B>\n";
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
  print '
<div class="menu_heading_closed" onclick="toggleMenu(\'103\')" id="heading103">
<b>Compare discovered motifs with databases (e.g. against Jaspar) or custom reference motifs</b></div>
<div id="menu103" class="menu_collapsible">';

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
  &MatrixDBcheckBox();

  print "<p/> ";
  print "<a href=''><b>Add known reference motifs or your own motif database:</b></a><br/>";
  print $query->filefield(-name=>'ref_motif',
			-size=>10);
  print "(matrices should be in <b>Transfac format</b>)";
  print "</fieldset><p/>";

  print '
</div>
</div>
<p class="clear"></p>';
}

##########################################
 sub Panel5  {
print '
<div class="menu_heading_closed" onclick="toggleMenu(\'104\')" id="heading104"> <b>Locate motifs and export as UCSC custom track</b> </div>
<div id="menu104" class="menu_collapsible">';



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

print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.matrix-scan.html'>Background model: Markov order</A>&nbsp;</B>\n";
	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='markov-scan'>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='1'>0</option>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='2'>1</option>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='3'>2</option>\n";
    $oligoPopup .=  "</SELECT>";
    print $oligoPopup;
    
print "<br/>";

### threshold (common to all programs)
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.matrix-scan.html'>Upper threshold on P-value</A>&nbsp;</B>\n";
print  $query->textfield(-name=>'uth_pval',
							      -default=>$default{uth_pval},
							      -size=>4);
print "<br/>";	



print "</fieldset><p/>";
#
### UCSC custom track
#
print "<fieldset>
<legend><b><a href='help.peak-motifs.html#tasks'>Visualize motifs </a></b></legend>";


print $query->checkbox(-name=>'bed_custom_track',
		       -checked=>$default{matrix-scan-quick},
		       -label=>'');  
print "&nbsp;<b>Visualize on UCSC genome browser</b> <a href='help.matrix-scan.html'></a>\n";
print "<br/>";


print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.peak-motifs.html'>BED file with peak coordinates</A>&nbsp;</B>\n";
print $query->filefield(-name=>'bed_file',
			-size=>10);

### threshold (common to all programs)
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.peak-motifs.html'>Assembly version (UCSC)</A>&nbsp;</B>\n";
print  $query->textfield(-name=>'assembly',
							      -default=>$default{assembly},
							      -size=>10);
print "<br/>";	

print "<br/>
</fieldset><p/>";

#########
print '	
</div>
</div>
<p class="clear"></p>';

 }
