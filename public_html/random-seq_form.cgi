#!/usr/bin/perl

## CVS
## added the possibility to specify the expected frequency for each nucleotide separately

#### this cgi script fills the HTML form for the program random-seq
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
$default{sequence_format} = "fasta";
$default{length} = 1000;
$default{repet} = 10;
$default{lw} = 50;

#$default{ATfreq} = 0.325;
#$default{CGfreq} = 0.175;
$default{Afreq} = 0.325;
$default{Tfreq} = 0.325;
$default{Cfreq} = 0.175;
$default{Gfreq} = 0.175;

$default{organism} = "Saccharomyces cerevisiae";
$default{oligo_size} = 6;
$default{oligopept_size} = 2;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

### print the form ###
&RSA_header("random sequence", "form");

### header
print "<center>";
print "Generate random DNA or protein sequences according to various probabilistic models<br>(independently distributed residues or Markov models)<P>\n";
print "</center>";


print $query->start_multipart_form(-action=>"random-seq.cgi");

print "<font face='Helvetica'>";

#### fragments
print "<fieldset><legend><b>Sequence number and sizes</b></legend>";

#print "<h2>Sequence number and sizes</h2>\n";

print "<UL>\n";

print " <A HREF='help.random-seq.html#length'>Sequence length</A> ";
print $query->textfield(-name=>'length',
			-default=>$default{length},
			-size=>7);

print " <A HREF='help.random-seq.html#repet'>Number of sequences</A> ";
print $query->textfield(-name=>'repet',
			-default=>$default{repet},
			-size=>5);

print "<BR>\n";

### sequence format
print "<A HREF='help.random-seq.html#formats'>Sequence format</A>&nbsp;";
print $query->popup_menu(-name=>'format',
			 -Values=>['fasta',
				   'IG',
				   'wconsensus',
				   'multi'],
			 -default=>$default{sequence_format});

print " <A HREF='help.random-seq.html#lw'>Line width</A> ";
print $query->textfield(-name=>'lw',
			-default=>$default{lw},
			-size=>5);

print "</UL>\n";


# File lengths
print "<p/><b> OR </b> <p/>";
print "<b><A HREF='help.random-genome-fragments.html#lf_length_file'>Use a set of sequences as template (same nb of fragments, same lengths): </a></b><br/> \n";
print "<div style='padding-left:30px'>";
&MultiSequenceChoice("Template sequences",1);
print "</div>";

print "</fieldset><p/>";

print "<fieldset><legend><a href='help.random-seq.html#alphabet'>Background model</a></legend>";
#print "<h2><a href='help.random-seq.html#alphabet'>Background model</a></h2>";

print "<ul>";

print "<b>Organism-specific Markov model</b> (Note: oligomer length = Markov order + 1)<br>";

print "<UL>";
&OrganismPopUp();

### oligo size
print "<br><INPUT TYPE='radio' NAME='bg_method' VALUE='upstream' checked>";
print "DNA sequences calibrated on non-coding upstream sequences)";
print "&nbsp"x3, "<b><a href='help.random-seq.html#oligo_size'>Oligonucleotide size</A>&nbsp;</b>\n";
print $query->popup_menu(-name=>'oligo_size',
			 -Values=>[1..8],
			 -default=>$default{oligo_size});

print "<br><INPUT TYPE='radio' NAME='bg_method' VALUE='protein'>";
print "Protein sequences calibrated on all proteins of this organism";
print "&nbsp"x3, "<b><a href='help.random-seq.html#oligopept_size'>Oligopeptide size</A>&nbsp;</b>\n";
print $query->popup_menu(-name=>'oligopept_size',
			 -Values=>[1..3],
			 -default=>$default{oligopept_size});

print "</UL>";

################################################################
## Independent nucleotides with distinct probabilities
print "<INPUT TYPE='radio' NAME='bg_method' VALUE='alphabet'>Independent nucleotides with distinct probabilities<BR>";

print "<UL>";
print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->td(["<B>A</B>",
					    $query->textfield(-name=>'Afreq',
							      -default=>$default{Afreq},
							      -size=>6),
					    "<B>T</B>",
					    $query->textfield(-name=>'Tfreq',
							      -default=>$default{Tfreq},
							      -size=>6)]),
				$query->td(["<B>C</B>",
					    $query->textfield(-name=>'Cfreq',
							      -default=>$default{Cfreq},
							      -size=>6),
					    "<B>G</B>",
					    $query->textfield(-name=>'Gfreq',
							      -default=>$default{Gfreq},
							      -size=>6)]),
			       ]));


print "</UL>";

print "<INPUT TYPE='radio' NAME='bg_method' VALUE='equi'>Independent and equiprobable nucleotides<BR>";
print "<p/>";
 print ("<b>Custom background model</b><br>");
 print "<input type='radio' NAME='bg_method' VALUE='file_upload'>";
    print "Upload your own background file (oligo-analysis format)", "&nbsp;"x5;
    print "<br/><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
  print $query->filefield(-name=>'upload_bgfile',
			    -default=>'starting value',
			    -size=>30,
			    -maxlength=>200);
  print "<p>\n";

print "</UL>";



print "</fieldset><p/>";


### send results by email or display on the browser
&SelectOutput("server");

### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.random-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_random-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

