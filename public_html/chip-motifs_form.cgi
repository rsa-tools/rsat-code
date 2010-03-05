#!/usr/bin/perl
#### this cgi script fills the HTML form for the program chip-seq-analysis
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
$default{oligo_length7}="checked";
$default{oligo-analysis}="checked";
$default{dyad-analysis}="checked";
$default{orm}="checked";


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
&RSA_header("ChIP-seq analysis", "form");
print "<CENTER>";
print "Pipeline for discovering motifs in ChIP-seq peak sequences.<P>\n";
print "<p><font color=red><b>Warning, this is still a prototype version</b></font>\n";
print "</CENTER>";
print "<BLOCKQUOTE>\n";

## demo description
print $default{demo_descr1};
print $default{demo_descr2};

print $query->start_multipart_form(-action=>"chip-seq-analysis.cgi");

################################################################
#### Start tab panel
print '

<div class="tabber">

     <div class="tabbertab">
	  <h2>Single dataset</h2>
'; 

## here is the content of the first panel
	  &Panel1();

# end of first panel
#start of second panel
print '</div>

<div class="tabbertab">
	  <h2>Two datasets (test + control)</h2>';
## here is the content of the second panel
	 &Panel2();

# end second panel 
# end the tab panel
print '</div>
</div>
';

################################################################
#### Change  parameters

print '
<br/>
<div>
<div class="menu_heading_closed" onclick="toggleMenu(\'96\')" id="heading96"> <b>Change default parameters</b> </div>
<div id="menu96" class="menu_collapsible">';
	

#### Tasks

##
print "<fieldset>
<legend><b><a href='help.chip-seq-analysis.html#tasks'>Compare motifs </a></b></legend>";

#
### tom-tom
#
print $query->checkbox(-name=>'tom-tom',
		       -checked=>$default{tom-tom},
		       -label=>'');  
print "&nbsp;<b>Compare discovered motifs with known motifs from databases</b> <a href='http://meme.nbcr.net/meme4_3_0/tomtom-intro.html'>[tom-tom from the MEME Suite]</a> <a href='http://genomebiology.com/2007/8/2/R24'> Gupta et al (2008) </a>\n";
print "<br/> ";

print "</fieldset><p/>";
#
### matrix-scan
#
print "<fieldset>
<legend><b><a href='help.chip-seq-analysis.html#tasks'>Locate motifs </a></b></legend>";
print $query->checkbox(-name=>'matrix-scan-quick',
		       -checked=>$default{matrix-scan-quick},
		       -label=>'');  
print "&nbsp;<b>Search putative binding sites in the peak sequences</b> <a href='help.matrix-scan.html'>[matrix-scan]</a>\n";

print "<br/>";

print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.matrix-scan.html'>Background model: Markov order</A>&nbsp;</B>\n";
	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='markov-scan'>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='1'>0</option>\n";
	$oligoPopup .=  "<OPTION VALUE='3'>3</option>\n";
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
### logo
#
print "<fieldset>
<legend><b><a href='help.chip-seq-analysis.html#tasks'>Visualize motifs </a></b></legend>";
print $query->checkbox(-name=>'logo',
		       -checked=>$default{logo},
		       -label=>'');  
print "&nbsp;<b>Visualize discovered motifs as logos</b> [seqlogo]\n";
print "<br/> ";

print $query->checkbox(-name=>'featuremap',
		       -checked=>$default{featuremap},
		       -label=>'');  
print "&nbsp;<b>Draw a feature-map to visualize putative binding sites</b> <a href='help.feature-map.html'>[feature-map]</a>\n";
print "<br/>";

print $query->checkbox(-name=>'avoir',
		       -checked=>$default{avoir},
		       -label=>'');  
print "&nbsp;<b>Draw positional profiles </b> [???]\n";
print "<br/>";

print $query->checkbox(-name=>'matrix-scan-quick',
		       -checked=>$default{matrix-scan-quick},
		       -label=>'');  
print "&nbsp;<b>Visualize on UCSC genome browser</b> <a href='help.matrix-scan.html'>[???]</a>\n";
print "<br/>";


print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.chip-seq-analysis.html'>BED file with peak coordinates</A>&nbsp;</B>\n";
print $query->filefield(-name=>'bed_file',
			-size=>10);

### threshold (common to all programs)
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.chip-seq-analysis.html'>Assembly version (UCSC)</A>&nbsp;</B>\n";
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

#print "<hr>";


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

my $descr1="<H4>Comment on the demonstration example 1 : </H4><blockquote class ='demo'>
In this demonstration, we .....<p/>
The program will return ...
</blockquote>";

print $query->start_multipart_form(-action=>"chip-seq-analysis_form.cgi");
$demo_seq=`cat demo_files/chip-seq-analysis_demo.fa`;
print "<TD><B>";
print $query->hidden(-name=>'demo_descr1',-default=>$descr1);
print $query->hidden(-name=>'sequence1',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format1',-default=>'fasta');
print $query->submit(-label=>"DEMO 1");
print "</B></TD>\n";
print $query->end_form;

my $descr2="<H4>Comment on the demonstration example 2 : </H4><blockquote class ='demo'>
In this demonstration, we .....<p/>
The program will return ...
</blockquote>";

### data for the demo with control
print $query->start_multipart_form(-action=>"chip-seq-analysis_form.cgi");
$demo_seq=`cat demo_files/chip-seq-analysis_demo.fa`;
print "<TD><B>";
print $query->hidden(-name=>'demo_descr2',-default=>$descr2);
print $query->hidden(-name=>'sequence1',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format1',-default=>'fasta');
print $query->hidden(-name=>'sequence2',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format2',-default=>'fasta');
print $query->submit(-label=>"DEMO 2");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.chip-seq-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);



#################### Internal functions #######################

sub Panel1 {

print "<b><font style='font-size:80%'><a href=''> (I only have coordinates in a BED file, how to get sequences ?)</a></font></b><p/> ";

&MultiSequenceChoice("Peak sequences",1);

print '
</br>
<div>

<div class="menu_heading_closed" onclick="toggleMenu(\'95\')" id="heading95"><b>Change motif discovery parameters </b> </div>
<div id="menu95" class="menu_collapsible">';
	

################# Motif discovery single input
print "<p/><fieldset>
<legend><b><a href='help.chip-seq-analysis.html#tasks'>Discover motifs </a></b></legend>";

#
### oligo analysis
#
print $query->checkbox(-name=>'oligo-analysis',
		       -checked=>$default{oligo-analysis},
		       -label=>'');  
print "&nbsp;<b>Discover over-represented words</b> <a href='help.oligo-analysis.html'> [oligo-analysis]</a>\n";
print "<br/>";

### oligo size
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
		       
#
### dyad-analysis
#   
print $query->checkbox(-name=>'dyad-analysis',
		       -checked=>$default{dyad-analysis},
		       -label=>''); 
print "&nbsp;<b>Discover over-represented spaced word pairs </b><a href='help.dyad-analysis.html'>[dyad-analysis] </a>\n";

print "<br/>";	
### dyad sizes and spacer
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.dyad-analysis.html#oligo_size'>Dyad length and spacer</A>&nbsp;</B>\n";
	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='dyad-option'>\n";
	$oligoPopup .=  "<OPTION  VALUE='4'>4 {0,20} 4</option>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='3'>3 {0,20} 3</option>\n";
    $oligoPopup .=  "</SELECT>";
    print $oligoPopup;
print "<br/>";

#
### ORM
#
print $query->checkbox(-name=>'orm',
		       -checked=>$default{orm},
		       -label=>'');  
print "&nbsp;<b>Discover words with positional biais</b> [ORM] [position-analysis]\n";
print "<br/>";

### orm oligo size
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.ORM.html#oligo_length'>Oligomer length</A>&nbsp;</B>\n";

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

print "&nbsp;&nbsp;&nbsp;&nbsp;<b><a href='help.ORM.html#window_width'>Fixed window of width</A>&nbsp;</B>\n";
	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='orm_window'>\n";
	$oligoPopup .=  "<OPTION VALUE='50'>50</option>\n";
	$oligoPopup .=  "<OPTION SELECTED VALUE='20'>20</option>\n";
	$oligoPopup .=  "<OPTION VALUE='10'>10</option>\n";
    $oligoPopup .=  "</SELECT>";
    print $oligoPopup;

print "<br/>";	
		       

### markov  order (common to all programs)
print "<br/> <b>Common options for above programs</b> <br/>";
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-analysis.html'>Background model: Markov order</A>&nbsp;</B>\n";
 	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='markov'>\n";
	$oligoPopup .=  "<OPTION  SELECTED VALUE='-2'>oligo length -2 </option>\n";
	$oligoPopup .=  "<OPTION VALUE='-3'>oligo length -3</option>\n";
    $oligoPopup .=  "</SELECT>";
    print $oligoPopup;
    
print "<br/>";

### threshold (common to all programs)
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-analysis.html#thresholds'>Lower threshold on significance</A>&nbsp;</B>\n";
print  $query->textfield(-name=>'lth_occ_sig',
							      -default=>$default{lth_occ_sig},
							      -size=>3);
print "<br/>";			 
print "</fieldset><p/>";

print '	
</div>
</div>
<p class="clear"></p>';
}

##########################################
sub Panel2 {
	
print '	<p/>';
	&MultiSequenceChoice("Peak sequences",1);

print '	<p/>';

	&MultiSequenceChoice("Control sequences",2);
print "<br/> ";
print '
<div>

<div class="menu_heading_closed" onclick="toggleMenu(\'97\')" id="heading97"><b>Change motif discovery parameters</b> </div>
<div id="menu97" class="menu_collapsible">';
	

	
print "<p/><fieldset>
<legend><b><a href='help.chip-seq-analysis.html#tasks'>Discover motifs </a></b></legend>";



print $query->checkbox(-name=>'oligo-diff',
		       -checked=>$default{oligo-diff},
		       -label=>'');  
print "&nbsp;<b>Discover words enriched in the peak sequences compared to the control sequences </b> <a href='help.oligo-diff.html'>[oligo-diff] </a>\n";
print "<br/> ";

## oligo size
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-diff.html#oligo_length'>Oligomer length</A>&nbsp;</B>\n";
my 	$oligoPopup = "";
    $oligoPopup .=  "<SELECT NAME='diff_length'>\n";
	$oligoPopup .=  "<OPTION  VALUE='6'>6</option>\n";
	$oligoPopup .=  "<OPTION VALUE='7'>7</option>\n";
	$oligoPopup .=  "<OPTION VALUE='8'>8</option>\n";
	$oligoPopup .=  "<OPTION SELECTED VALUE='6-8'>6 to 8</option>\n";
    $oligoPopup .=  "</SELECT>";
    print $oligoPopup;

print "<br/>";
### threshold 
print "&nbsp;&nbsp;&nbsp;&nbsp;<B><A HREF='help.oligo-analysis.html#thresholds'>Lower threshold on significance</A>&nbsp;</B>\n";
print  $query->textfield(-name=>'lth_occ_sig',
							      -default=>$default{lth_occ_sig},
							      -size=>3);
print "<br/>";


print "</fieldset><p/>";

print '	
</div>
</div>
<p class="clear"></p>';
}