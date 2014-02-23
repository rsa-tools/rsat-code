#!/usr/bin/perl
## this cgi script fills the HTML form for the program position-analysis
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

## Read the CGI query
$query = new CGI;

## Default values for filling the form
$default{output} = "email";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{oligo_length} = 6;
$default{class_interval} = 20;
$default{strand} = "single strand";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{purge} = 'checked';
$default{origin} = "start";
$default{offset} = "0";
$default{center} = "";
$default{demo_descr} = "";

## Output fields
$default{return_chi} = 'checked';
$default{return_rank} = 'checked';
$default{return_distrib} = 'checked';
$default{return_exp} = '';
$default{return_clusters} = '';
$default{return_matrices} = 'checked';
$default{return_graphs} = '';

## Thresholds and filtering
$default{sort} = 'checked';
$default{check} = 'checked';
$default{filter} = '';
$default{lth_occ} = "1";
$default{lth_chi} = "none";
$default{lth_sig} = "0";
$default{clust_nb} = "5";
$default{max_asmb_nb} = "5";

## Print the form
&RSA_header("position-analysis", "form");

print "<blockquote>";

print q {

   Calculates the positional distribution of oligonucleotides in a set
   of sequences, and detects those which significantly discard from a
   homogeneous distribution.</p>

   <b>Warning</b>: this program is useful for large data sets (some
   hundreds or thousands of sequences), pre-aligned on some signal
   (e.g. start codon).</p>

   <b>Reference</b>: van Helden, J., del Olmo, M. and Perez-Ortin,
   J. E. (2000). Statistical analysis of yeast genomic downstream
   sequences reveals putative polyadenylation signals. Nucleic Acids
   Res 28, 1000-10.<p>

};
print "<hr>";

#&ListParameters;

## Replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

## Print the demo description
print $default{demo_descr};

## Form
print $query->start_multipart_form(-action=>"position-analysis.cgi");

#&MultiSequenceChoice("Sequences",1);
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

## Purge sequences
print $query->checkbox(-name=>'purge',
 		       -checked=>$default{purge},
 		       -label=>'');
print "&nbsp;<A HREF='help.position-analysis.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";

print "<HR width=550 align=left>\n";


## Oligo size
print "<B><A HREF='help.position-analysis.html#oligo_length'>Oligonucleotide size</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_length',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_length});

## Prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "&nbsp;<A HREF='help.position-analysis.html#noov'><B>prevent overlapping matches</B></A>";
print "<BR>\n";


## Strand
print "<B><A HREF='help.position-analysis.html#count_strands'>Count on</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

## Group patterns by pairs of reverse complement
print $query->checkbox(-name=>'grouprc',
		       -checked=>$default{grouprc},
		       -label=>'');
print "&nbsp;<A HREF='help.position-analysis.html#grouprc'><B>return reverse complements together in the output</B></A>";
print "<BR>";


print "<HR width=550 align=left>\n";

print "<B><A HREF='help.position-analysis.html#class_grouping'>Positions</A>&nbsp;</B>\n";

## Class interval
print "<B>","&nbsp"x5,"<A HREF='help.position-analysis.html#class_interval'>window size</A>&nbsp;</B>\n";
print $query->textfield(-name=>'class_interval',
			-default=>$default{class_interval},
			-size=>3);

## origin for calculating positions
print "&nbsp;"x4,  "<A HREF='help.position-analysis.html#origin'><B>Origin</B></A>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['start',
				   'center',
				   'end'],
			 -default=>$default{origin});

## Offset for calculating positions
print "&nbsp;"x4,  "<A HREF='help.position-analysis.html#offset'><B>Offset</B></A>\n";
print $query->textfield(-name=>'offset',
			-default=>$default{offset},
			-size=>8);

## Output fields and thresholds
print "<p><b>Output fields</b></p>\n";
print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			 [
			  $query->th([" <A HREF='help.position-analysis.html#return'>Return</A> ",
				   " <A HREF='help.position-analysis.html#thresholds'>Lower<BR>Threshold</A> ",
				   " <A HREF='help.position-analysis.html#thresholds'>Upper<BR>Threshold</A> "
				      ]),

			  ## occurrences
			  $query->td(["Occurrences",
				   $query->textfield(-name=>'lth_occ',
						     -default=>$default{lth_occ},
						     -size=>5),
				   '']),

			  ## chi-square
			  $query->td([$query->checkbox(-name=>'return_chi',
						    -checked=>$default{return_chi},
						    -label=>' Chi2 statistics '),
				   $query->textfield(-name=>'lth_chi',
						     -default=>$default{lth_chi},
						     -size=>5)]),

			  ## significance
			  $query->td([join("","&nbsp"x10,"chi2 significance"),
#				      $query->checkbox(-name=>'return_sig',
#						       -checked=>$default{return_sig},
#						       -label=>' Chi2 significance '),
				      $query->textfield(-name=>'lth_sig',
							-default=>$default{lth_sig},
							-size=>5)]),

			  ## rank
			  $query->td([$query->checkbox(-name=>'return_rank',
						    -checked=>$default{return_rank},
						    -label=>' Rank '),
				      ''
				      ]),

			  ## position distribution
			  $query->td([$query->checkbox(-name=>'return_distrib',
						       -checked=>$default{return_distrib},
						       -label=>' Position distribution '),
				      '']),


			  ## Expected occurrences
			  $query->td([$query->checkbox(-name=>'return_exp',
						       -checked=>$default{return_exp},
						       -label=>' Expected occurrences'),
				      '']),

			  ## clusters
			  $query->td([$query->checkbox(-name=>'return_clusters',
						       -checked=>$default{return_clusters},
						       -label=>' Clusters (words clustered by position profiles)'),
				      "&nbsp",
				      $query->textfield(-name=>'clust_nb',
							-default=>$default{clust_nb},
							-size=>5)]),

			  ## matrices
			  $query->td([$query->checkbox(-name=>'return_matrices',
						       -checked=>$default{return_matrices},
						       -label=>' Position-specific scoring matrices'),
				      "&nbsp",
				      $query->textfield(-name=>'max_asmb_nb',
							-default=>$default{max_asmb_nb},
							-size=>5)]),
#				      '']),

			  ## graphs
			  $query->td([$query->checkbox(-name=>'return_graphs',
						       -checked=>$default{return_graphs},
						       -label=>' Graphs '),
				      '']),



			 ]
			)
		);
print "</BLOCKQUOTE>\n";

## check applicability condition for the chi2 test
print "<P>";
print "&nbsp;<A HREF='help.position-analysis.html#applicability'><B>Applicability condition for the chi2 test</B></A>&nbsp;";
print $query->checkbox(-name=>'check',
		       -checked=>$default{check},
		       -label=>' check');
print $query->checkbox(-name=>'filter',
		       -checked=>$default{filter},
		       -label=>' filter ');
print "<BR>\n";


## Sort the patterns accoring to the best score
print "<P>";
print "&nbsp;<A HREF='help.position-analysis.html#sort'><B>Sort patterns by significance.</B></A>&nbsp;";
print $query->checkbox(-name=>'sort',
		       -checked=>$default{sort},
		       -label=>'');

print "<HR width=550 align=left>\n";

#print "<font color=red><B>Warning!</B> position-analysis is time-consuming. If the result is not displayed after 5 minutes, try email output.</font><BR>\n";

## Send results by email or display on the browser
&SelectOutput($default{output});

## Action buttons
print "<UL><UL><TABLE class = 'formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

## Data for the demo 
print $query->start_multipart_form(-action=>"position-analysis_form.cgi");
$demo_seq_file = "$ENV{RSAT}/public_html/demo_files/SWEMBL_mmus_HNF4A_vs_mmus_Input_peaks_R0.05_nof_200bp.fasta.gz";
##$demo_seq_file = "$ENV{RSAT}/public_html/data/demo_files/all_yeast_downstream_200bp.fasta.gz";
##$demo_seq_file = "$ENV{RSAT}/public_html/demo_files/Mycoplasma_genitalium_upstream_-30_+29.fasta.gz";
$demo_seq = `gunzip -c $demo_seq_file`;
##$demo_url= $ENV{rsat_www}."/demo_files/peak-motifs_GSM559652_heart_p300_1000peaks.fa";



################################################################
## Description of the demo test case
my $demo_descr = "<H4>Comment on the demonstration example:</H4>\n";
$demo_descr .= "<blockquote class ='demo'>";
$demo_descr .= "<p>In this demonstration, we run position-analysis to detect position-biased motifs in a peak set from a ChIP-seq experiment for the transcription factor HNF4a, in mouse liver cells.";
$demo_descr .= "<br>The peak set was obtained by running SWEMBL with stringent parameters (R=0.05) on the aligned reads. We filtered out the peaks smaller than 220bp, and clipped the remaining sequences to 110bp on each side of the peak center.\n";
$demo_descr .= "<br>Data source: Schmidt et al. Five-vertebrate ChIP-seq reveals the evolutionary dynamics of transcription factor binding. Science (2010) vol. 328 (5981) pp. 1036-40.";
$demo_descr .= "</p>";
$demo_descr .= "</blockquote>";

print "<TD><B>";
print $query->hidden(-name=>'demo_descr',-default=>$demo_descr);
print $query->hidden(-name=>'sequence',-default=>$demo_seq);
print $query->hidden(-name=>'sequence_format',-default=>$seq_format);
print $query->hidden(-name=>'oligo_length',-default=>'6');
print $query->hidden(-name=>'class_interval',-default=>'10');
print $query->hidden(-name=>'origin',-default=>'center');
print $query->hidden(-name=>'return_chi',-default=>'on');
print $query->hidden(-name=>'lth_sig',-default=>'5');
print $query->hidden(-name=>'return_clusters',-default=>'on');
print $query->hidden(-name=>'clust_nb',-default=>'2');
print $query->hidden(-name=>'return_matrices',-default=>'on');
print $query->hidden(-name=>'max_asmb_nb',-default=>'2');
print $query->hidden(-name=>'return_graphs',-default=>'on');
print $query->hidden(-name=>'strand',-default=>'both strands');
print $query->hidden(-name=>'output',-default=>'email');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


#print "<TD><B><A HREF='demo.position-analysis.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='help.position-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='tutorials/tut_position-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


