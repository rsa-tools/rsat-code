#!/usr/bin/perl
#### this cgi script fills the HTML form for the program position-analysis
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
$default{output} = "display";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{oligo_length} = 6;
$default{class_interval} = 20;
$default{strand} = "single strand";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{purge} = 'checked';
$default{origin} = "-0";

#### return values
$default{return_chi} = 'checked';
$default{return_rank} = 'checked';
$default{return_distrib} = 'checked';
$default{return_exp} = '';
$default{return_graph} = '';

### thresholds and filtering
$default{sort} = 'checked';
$default{check} = 'checked';
$default{filter} = '';
$default{lth} = "0";
$default{oth} = "1";

### print the form ###
&RSA_header("position-analysis", "form");

print "<blockquote>";

print q {

   Calculates the positional distribution of oligonucleotides in a set
   of sequences, and detects those which significantly discard from a
   homogeneous distribution.

   <b>Warning</b>: this program is useful for large data sets (some
   hundreds or thousands of sequences), pre-aaligned on some signal
   (e.g. start codon).

};

print "<HR>";

#&ListParameters;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

print $query->start_multipart_form(-action=>"position-analysis.cgi");


print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td([&SequenceChoice()])
			]),
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td(["<B><A HREF='help.position-analysis.html#sequence_type'>Sequence type</A></B>".
			       $query->popup_menu(-name=>'sequence_type',
						  -Values=>["dna","protein","other"],
						  -default=>$default{sequence_type})
			       ])
			])

		 );

#### purge sequences
print $query->checkbox(-name=>'purge',
 		       -checked=>$default{purge},
 		       -label=>'');
print "&nbsp;<A HREF='help.position-analysis.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";

print "<HR width=550 align=left>\n";


### OLIGOy size
print "<B><A HREF='help.position-analysis.html#oligo_length'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_length',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_length});

### prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "&nbsp;<A HREF='help.position-analysis.html#noov'><B>prevent overlapping matches</B></A>";
print "<BR>\n";


### strand ###
print "<B><A HREF='help.position-analysis.html#count_strands'>Count on</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

#### group patterns by pairs of reverse complement
print $query->checkbox(-name=>'grouprc',
		       -checked=>$default{grouprc},
		       -label=>'');
print "&nbsp;<A HREF='help.position-analysis.html#grouprc'><B>return reverse complements together in the output</B></A>";
print "<BR>";


print "<HR width=550 align=left>\n";

print "<B><A HREF='help.position-analysis.html#class_grouping'>Class grouping</A>&nbsp;</B>\n";

### class interval
print "<B><A HREF='help.position-analysis.html#class_interval'>interval</A>&nbsp;</B>\n";
print $query->textfield(-name=>'class_interval',
			-default=>$default{class_interval},
			-size=>3);

#### origin
print "<B><A HREF='help.position-analysis.html#origin'>origin</A>&nbsp;</B>\n";
print $query->textfield(-name=>'origin',
			-default=>$default{origin},
			-size=>3);

#### table with all the statistics and thresholds
print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			 [
			  $query->th([" <A HREF='help.position-analysis.html#return'>Return</A> ",
				   " <A HREF='help.position-analysis.html#thresholds'>Lower<BR>Threshold</A> ",
#				   " <A HREF='help.position-analysis.html#thresholds'>Upper<BR>Threshold</A> "
				      ]),

			  ### occurrences
			  $query->td(["occurrences",
				   $query->textfield(-name=>'oth',
						     -default=>$default{oth},
						     -size=>5),
				   '']),

			  ### chi-square
			  $query->td([$query->checkbox(-name=>'return_chi',
						    -checked=>$default{return_chi},
						    -label=>' Chi2 '),
				   $query->textfield(-name=>'lth',
						     -default=>$default{lth},
						     -size=>5)]),

			  ### rank
			  $query->td([$query->checkbox(-name=>'return_rank',
						    -checked=>$default{return_rank},
						    -label=>' Rank '),
				      ''
				      ]),

			  ### position distribution
			  $query->td([$query->checkbox(-name=>'return_distrib',
						       -checked=>$default{return_distrib},
						       -label=>' Position distribution '),
				      '']),


			  ### expected distribution
			  $query->td([$query->checkbox(-name=>'return_exp',
						       -checked=>$default{return_exp},
						       -label=>' Expected distribution (homogeneous model) '),
				      '']),

			  ### graphs
#			  $query->td([$query->checkbox(-name=>'return_graph',
#						       -checked=>$default{return_graph},
#						       -label=>' Graphs '),
#				      '']),



			 ]
			)
		);
print "</BLOCKQUOTE>\n";

### check applicability condition for the chi2 test
print "<P>";
print "&nbsp;<A HREF='help.position-analysis.html#applicability'><B>Applicability condition for the chi2 test</B></A>&nbsp;";
print $query->checkbox(-name=>'check',
		       -checked=>$default{check},
		       -label=>' check');
print $query->checkbox(-name=>'filter',
		       -checked=>$default{filter},
		       -label=>' filter ');
print "<BR>\n";


#### sort the patterns accoring to the best score
print "<P>";
print "&nbsp;<A HREF='help.position-analysis.html#sort'><B>Sort patterns according to the score</B></A>&nbsp;";
print $query->checkbox(-name=>'sort',
		       -checked=>$default{sort},
		       -label=>'');

print "<HR width=550 align=left>\n";

print "<font color=red><B>Warning !</B>position-analysis is time-consuming. We recommend email output.</font><BR>\n";


### send results by email or display on the browser
&SelectOutput($default{output});


### action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
print $query->start_multipart_form(-action=>"position-analysis_form.cgi");
### $demo_seq_file = "$ENV{RSAT}/public_html/data/demo_files/all_yeast_downstream_200bp.fasta.gz";
$demo_seq_file = "$ENV{RSAT}/public_html/demo_files/Mycoplasma_genitalium_upstream_-30_+29.fasta.gz";

print "<TD><B>";
print $query->hidden(-name=>'sequence_file',-default=>$demo_seq_file);
print $query->hidden(-name=>'sequence_format',-default=>'fasta');
print $query->hidden(-name=>'output',-default=>'display');
print $query->hidden(-name=>'oligo_length',-default=>3);
print $query->hidden(-name=>'organism',-default=>'Saccharomyces cerevisiae');
print $query->hidden(-name=>'class_interval',-default=>'3');
print $query->hidden(-name=>'origin',-default=>'-30');
print $query->hidden(-name=>'chi',-default=>'50');
#print $query->hidden(-name=>'filter',-default=>'checked');
print $query->hidden(-name=>'strand',-default=>'single strand');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.position-analysis.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.position-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_position-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


