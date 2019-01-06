#!/usr/bin/env perl
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
$default{min_pos} = "none";
$default{max_pos} = "none";
$default{strand} = "single strand";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{purge} = 'checked';
$default{origin} = "start";
$default{offset} = "0";
$default{center} = "";
$default{demo_descr} = "";

## Background model
$default{bg_method} = "Homogeneity";
$default{markov_order} = "1";

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
$default{max_asmb_per_cluster} = "5";

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
#print $default{demo_descr};
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr'></div>";

## Form
print $query->start_multipart_form(-action=>"position-analysis.cgi", -id=>"form");

#&MultiSequenceChoice("Sequences",1);
print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td([&SequenceChoice()])
			]),
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td(["<b><a class='iframe' href='help.position-analysis.html#sequence_type'>Sequence type</a></b>".
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
print "&nbsp;<a class='iframe' href='help.position-analysis.html#purge'><b>purge sequences (highly recommended)</b></a>";
print "<BR>";

print "<HR width=550 align=left>\n";


## Oligo size
print "<b><a class='iframe' href='help.position-analysis.html#oligo_length'>Oligonucleotide size</a>&nbsp;</b>\n";
print $query->popup_menu(-name=>'oligo_length',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_length});

## Prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "&nbsp;<a class='iframe' href='help.position-analysis.html#noov'><b>prevent overlapping matches</b></a>";
print "<BR>\n";


## Strand
print "<b><a class='iframe' href='help.position-analysis.html#count_strands'>Count on</a>&nbsp;</b>\n";
print $query->popup_menu(-name=>'strand',-id=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

## Group patterns by pairs of reverse complement
print $query->checkbox(-name=>'grouprc',
		       -checked=>$default{grouprc},
		       -label=>'');
print "&nbsp;<a class='iframe' href='help.position-analysis.html#grouprc'><b>return reverse complements together in the output</b></a>";
print "<BR>";


print "<b><a class='iframe' href='help.position-analysis.html#class_grouping'>Positions</a>&nbsp;</b>\n";

## Class interval
print "<b>","&nbsp"x5,"<a class='iframe' href='help.position-analysis.html#class_interval'>window size</a>&nbsp;</b>\n";
print $query->textfield(-name=>'class_interval',
			-default=>$default{class_interval},
			-size=>3);

## Origin for calculating positions
print "&nbsp;"x4,  "<a class='iframe' href='help.position-analysis.html#origin'><b>Origin</b></a>\n";
print $query->popup_menu(-name=>'origin',
			 -Values=>['start',
				   'center',
				   'end'],
			 -default=>$default{origin});

## Offset for calculating positions
print "&nbsp;"x4,  "<a class='iframe' href='help.position-analysis.html#offset'><b>Offset</b></a>\n";
print $query->textfield(-name=>'offset',
			-default=>$default{offset},
			-size=>8);

## max and min positions for computing the chi2 statistics
print "<br><b>", "<a class='iframe' href='help.position-analysis.html#min_mos'>Position limits for chi2: </a>&nbsp;</b>\n";
print "&nbsp;"x5, "min:&nbsp";
print $query->textfield(-name=>'min_pos',
			-default=>$default{min_pos},
			-size=>4);

print "&nbsp;"x5, "max:&nbsp";
print $query->textfield(-name=>'max_pos',
			-default=>$default{min_pos},
			-size=>4);



################################################################
## Background model

print "<hr>";
print "<p><b>Background model</b></p>\n";
print "<blockquote>\n";

## Homogeneous distribution of the k-mers across the windows
print ("<br><input type='radio' NAME='bg_method' VALUE='Homogeneous' checked>", "Homogeneous distribution across windows.");

## Markov model
print ("<br><input type='radio' NAME='bg_method' VALUE='Markov model'>", "Window-specific background model");

print "&nbsp;"x10,"Markov order:&nbsp;";
print $query->popup_menu(-name=>'markov_order',
			 -Values=>[0..2],
			 -default=>$default{markov_order});
print "</blockquote>";

################################################################
## Output fields and thresholds
print "<hr>";
print "<p><b>Output fields</b></p>\n";
print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			 [
			  $query->th([" <a class='iframe' href='help.position-analysis.html#return'>Return</a> ",
				   " <a class='iframe' href='help.position-analysis.html#thresholds'>Lower<BR>Threshold</a> ",
				   " <a class='iframe' href='help.position-analysis.html#thresholds'>Upper<BR>Threshold</a> "
				      ]),

			  ## occurrences
			  $query->td(["Occurrences",
				   $query->textfield(-name=>'lth_occ',
						     -default=>$default{lth_occ},
						     -size=>5),
				   '']),

			  ## chi-square
			  $query->td([$query->checkbox(-name=>'return_chi', -id=>'return_chi',
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
			  $query->td([$query->checkbox(-name=>'return_clusters',-id=>'return_clusters',
						       -checked=>$default{return_clusters},
						       -label=>' Clusters (words clustered by position profiles)'),
				      "&nbsp",
				      $query->textfield(-name=>'clust_nb',
							-default=>$default{clust_nb},
							-size=>5)]),

			  ## matrices
			  $query->td([$query->checkbox(-name=>'return_matrices',-id=>'return_matrices',
						       -checked=>$default{return_matrices},
						       -label=>' Position-specific scoring matrices'),
				      "&nbsp",
				      $query->textfield(-name=>'max_asmb_per_cluster',
							-default=>$default{max_asmb_per_cluster},
							-size=>5)]),
#				      '']),

			  ## graphs
			  $query->td([$query->checkbox(-name=>'return_graphs',-id=>'return_graphs',
						       -checked=>$default{return_graphs},
						       -label=>' Graphs '),
				      '']),



			 ]
			)
		);
print "</BLOCKQUOTE>\n";

## check applicability condition for the chi2 test
print "<P>";
print "&nbsp;<a class='iframe' href='help.position-analysis.html#applicability'><b>Applicability condition for the chi2 test</b></a>&nbsp;";
print $query->checkbox(-name=>'check',
		       -checked=>$default{check},
		       -label=>' check');
print $query->checkbox(-name=>'filter',
		       -checked=>$default{filter},
		       -label=>' filter ');
print "<BR>\n";


## Sort the patterns accoring to the best score
print "<P>";
print "&nbsp;<a class='iframe' href='help.position-analysis.html#sort'><b>Sort patterns by significance.</b></a>&nbsp;";
print $query->checkbox(-name=>'sort',
		       -checked=>$default{sort},
		       -label=>'');

print "<HR width=550 align=left>\n";

#print "<font color=red><b>Warning!</b> position-analysis is time-consuming. If the result is not displayed after 5 minutes, try email output.</font><BR>\n";

## Send results by email or display on the browser
&SelectOutput($default{output});

## Action buttons
print "<ul><ul><table class = 'formbutton'>\n";
print "<tr VALIGN=MIDDLE>\n";
print "<td>", $query->submit(-label=>"GO"), "</td>\n";
print "<td>", $query->reset(-id=>"reset"), "</td>\n";
print $query->end_form;

## Data for the demo 

$demo_seq_file = "$ENV{RSAT}/public_html/demo_files/SWEMBL_mmus_HNF4A_vs_mmus_Input_peaks_R0.05_nof_200bp.fasta.gz";
$demo_seq_raw = `gunzip -c $demo_seq_file`;
$demo_seq = join("\\n", split(/\n/, $demo_seq_raw));

################################################################
## Description of the demo test case
print '<script>
function setDemo(demo_seq){
    $("#reset").trigger("click");
    descr = "<H4>Comment on the demonstration example:</H4>\n<blockquote class =\'demo\'><p>In this demonstration, we run position-analysis to detect position-biased motifs in a peak set from a ChIP-seq experiment for the transcription factor HNF4a, in mouse liver cells.<br>The peak set was obtained by running SWEMBL with stringent parameters (R=0.05) on the aligned reads. We filtered out the peaks smaller than 220bp, and clipped the remaining sequences to 110bp on each side of the peak center.\n<br>Data source: Schmidt et al. Five-vertebrate ChIP-seq reveals the evolutionary dynamics of transcription factor binding. Science (2010) vol. 328 (5981) pp. 1036-40.</p></blockquote>";
    
    demo_descr.innerHTML = descr;
    demo.value = descr;
    
    sequence.value = demo_seq;
    $("[name=\'oligo_length\']").val(6);
    $("[name=\'class_interval\']").val(10);
    $("[name=\'origin\']").val("center");
    $("#return_chi").prop("checked",true);
    $("[name=\'lth_sig\']").val(5);
    $("#return_clusters").prop("checked",true);
    $("[name=\'clust_nb\']").val(2);
    $("#return_matrices").prop("checked",true);
    $("[name=\'max_asmb_per_cluster\']").val(2);
    $("#return_graphs").prop("checked",true);
    $("#strand").val("both strands");
}
</script>';


print "<td><b>";
print '<button type="button" onclick="setDemo('."'$demo_seq'".')">DEMO</button>';
print "</b></td>\n";


#print "<td><b><a href='demo.position-analysis.html'>DEMO</a></b></td>\n";
print "<td><b><a class='iframe' href='help.position-analysis.html'>MANUAL</a></b></td>\n";
print "<td><b><a class='iframe' href='tutorials/tut_position-analysis.html'>TUTORIAL</a></b></td>\n";
print "<td><b><a href='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</a></b></td>\n";
print "</tr></table></ul></ul>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


