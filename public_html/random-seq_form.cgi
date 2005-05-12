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
require "RSA.cgi.lib";
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

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

### print the form ###
&RSA_header("random sequence");

### head
print "<CENTER>";
print "Generate random DNA sequences according to various probabilistic models (Markov chains or independently distributed nucleotides)<P>\n";
print "</CENTER>";


print $query->start_multipart_form(-action=>"random-seq.cgi");

print "<FONT FACE='Helvetica'>";


print "<B>General options </B>\n";

print "<UL>\n";

print " <A HREF='help.random-seq.html#length'>Sequence length</A> ";
print $query->textfield(-name=>'length',
			-default=>$default{length},
			-size=>10);

print " <A HREF='help.random-seq.html#repet'>Number of sequences</A> ";
print $query->textfield(-name=>'repet',
			-default=>$default{repet},
			-size=>10);

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
			-size=>10);

print "</UL>\n";


print "<H4><A HREF='help.random-seq.html#alphabet'>Nucleotide probabilities</a></H4>";

print "<UL>";

print "<INPUT TYPE='radio' NAME='proba' VALUE='upstream' checked>Markov chain (calibrated on oligonucleotide frequencies in non-coding upstream sequences)<BR>";

print "<UL>";
&OrganismPopUp();

### oligo size
print "<B><A HREF='help.random-seq.html#oligo_size'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_size',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_size});

print "<font size=-1><b>Note:</b> Markov order = oligonucleotide length minus 1</font>";

print "</UL>";


################################################################
## Independent nucleotides with distrinct probabilities
print "<INPUT TYPE='radio' NAME='proba' VALUE='alphabet'>Independent nucleotides with distinct probabilities<BR>";

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


# print "<B>A:T</B>&nbsp;\n";
# print $query->textfield(-name=>'ATfreq',
# 			-default=>$default{ATfreq},
# 			-size=>6);

# print "<B>&nbsp;C:G</B>&nbsp;\n";
# print $query->textfield(-name=>'CGfreq',
# 			-default=>$default{CGfreq},
# 			-size=>6);

print "</UL>";

print "<INPUT TYPE='radio' NAME='proba' VALUE='equi'>Independent and equiprobable nucleotides<BR>";

print "</UL>";


### send results by email or display on the browser
&SelectOutput();

### action buttons
print "<UL><UL><TABLE>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.random-seq.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_random-seq.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";

print $query->end_html;

exit(0);

