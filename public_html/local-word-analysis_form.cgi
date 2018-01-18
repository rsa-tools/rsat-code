#!/usr/bin/perl
#### this cgi script fills the HTML form for the program local-word-analysis
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
$default{organism} = "";
$default{title} = "";
$default{sequence} = "";
$default{sequence_format} = "fasta";
$default{sequence_file} = "";
$default{upload_freq_file} = "";
#$default{sequence_type} = "dna";
$default{oligo_length} = 6;
$default{background} = "upstream-noorf";
$default{markov_order} = 2;
$default{strand} = "both strands";
$default{noov} = 'checked';
$default{grouprc} = 'checked';
$default{window_group} = 'checked';
$default{purge} = 'checked';
#$default{side} = 'over-represented';
$default{align} = 'right';
$default{freq_estimate} = "background";
$default{bg_level} = "organism";


$default{rank} = 'checked';
$default{lth_rank} = "none";
$default{uth_rank} = 50;

$default{lth_w_rank} = "none";
$default{uth_w_rank} = 1;


$default{occ} = 'checked';
$default{lth_occ} = 2;
$default{uth_occ} = "none";

$default{proba} = 'checked';
$default{lth_occ_P} = "none";
$default{uth_occ_P} = "none";

$default{eval} = 'checked';
$default{lth_occ_E} = "none";
$default{uth_occ_E} = "none";

$default{lth_occ_sig} = "0";
$default{uth_occ_sig} = "none";

$default{freq} = '';
$default{lth_observed_freq} = "none";
$default{uth_observed_freq} = "none";


$default{window_width} = '50';
$default{bg_window_width} = '500';

#$default{return}="fields";


### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 
$checked{$default{freq_estimate}} = "CHECKED";

### print the form ###
&RSA_header("local-word-analysis", "form");
print '<style><!-- textarea {height: 100px; width: 550px;}--></style>';
print "<CENTER>";
print "Analysis of oligonucleotide occurrences in a set of DNA sequences.\n";
print "<br>Program developed by <a href='mailto:defrance\@bigre.ulb.ac.be (Matthieu Defrance)'>Matthieu Defrance</A><P>";
print "</center>";

&ListDefaultParameters() if ($ENV{rsat_echo} >= 2);

print $query->start_multipart_form(-action=>"local-word-analysis.cgi");

print "<hr>";
print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>0},
	       $query->Tr({-align=>left,-valign=>TOP},
		       [
		      $query->td([&SequenceChoice()])
			]),
	       $query->Tr({-align=>left,-valign=>TOP},
		       [

			])

		 );

#### purge sequences
print $query->checkbox(-name=>'purge',
		       -checked=>$default{purge},
		       -label=>'');
print "<a class='iframe' href='help.local-word-analysis.html#purge'><b>Purge sequences (highly recommended)</b></a>";
print "<br />";

print "<hr width=550 align=left />\n";


print "<b>Search parameters (motif)</b><br />\n";
print '<input type="radio" name="oligotype" value="oligo" checked="checked"/>';
### oligo size
print "<b><a class='iframe' href='help.local-word-analysis.html#oligo_length'>Oligonucleotides of length </a></b>\n";
print $query->popup_menu(-name=>'oligo_length',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_length});

print '<br />';
print '<input type="radio" name="oligotype" value="dyad" />';
print "<b><a class='iframe' href='help.local-word-analysis.html#oligo_length'>Dyads with monad of length </a> </b>\n";
print $query->popup_menu(-name=>'monad_length',
			 -Values=>[1,2,3],
			 -default=>3);
 print "<b><a class='iframe' href='help.local-word-analysis.html#oligo_length'> and spacing from </a></b>\n";
 print $query->popup_menu(-name=>'spacing_a',
			 -Values=>[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
 			 -default=>0);
 print "<b><a class='iframe' href='help.local-word-analysis.html#oligo_length'> to </a></b>\n";
 print $query->popup_menu(-name=>'spacing_b',
			 -Values=>[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
 			 -default=>20);



print '<br /><br />';

### prevent overlapping matches of the same pattern
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print '<a class="iframe" href="help.local-word-analysis.html#noov""><b>Prevent overlapping matches</b></a>';

### strand ###
print '&nbsp;&nbsp&nbsp;&nbsp<b><a class="iframe" href="help.local-word-analysis.html#count_strands">Count on </a></b>';
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});

### align ###
print '&nbsp;&nbsp&nbsp;&nbsp<b><a class="iframe" href="help.local-word-analysis.html#align">Align </a></b>';
print $query->popup_menu(-name=>'align',
			 -Values=>['right',
 				  'left'],
 			 -default=>$default{align});


#### group patterns by pairs of reverse complement
#print $query->checkbox(-name=>'grouprc',
#		       -checked=>$default{grouprc},
#		       -label=>'');
#print "&nbsp;<A HREF='help.local-word-analysis.html#grouprc'><B>return reverse complements together in the output</B></A>";
#print "<BR>";
print '<br /><br />';
print '<input type="radio" name="windowtype" value="no" />';
print '<b><a class="iframe" href="help.local-word-analysis.html#window_width">No window (like oligo-analysis or dyad-analysis)</a></b>';

print '<br />';

print '<input type="radio" name="windowtype" value="fixed" checked="checked"/>';
print '<b><a class="iframe" href="help.local-word-analysis.html#window_width">Fixed window of width</a></b>';
print $query->textfield(-name=>'window_width',-id=>'window_width',
			-default=>$default{window_width},
			-size=>5);
print ' ';
print $query->checkbox(-name=>'window_group',
		       -checked=>$default{window_group},
		       -label=>'');
print '<b><a class="iframe" href="help.local-word-analysis.html#window_width">Group windows</a></b>';


print '<br />';
print '<input type="radio" name="windowtype" value="variable" />';
print '<b><a class="iframe" href="help.local-word-analysis.html#window_width">Variable window width</a></b> (Warning ! this can be time consuming)';

print '<hr width="550" align="left">';

################################################################
## Background model
#print '<br /><br />';
print '<b><a class="iframe" href="help.local-word-analysis.html#bgwindow_width">Background window width</a></b>';
print $query->textfield(-name=>'bg_window_width',-id=>'bg_window_width',
			-default=>$default{bg_window_width},
			-size=>5);

print '<br />';
&PrintOligoBackgroundOptions();

################################################################
print "<HR width=550 align=left>\n";

#print "<A HREF='help.local-word-analysis.html#exp_freq'><B>Expected frequency</B></A>&nbsp;";
#print $query->radio_group(-name=>'freq_estimate',
#			  -Values=>['Equiprobable residues',
#				    'Residue frequencies from input sequence',
#				    'Markov Model (higher order dependencies)',
#				    'Lexicon partitioning',
#				    'Oligo frequencies from all intergenic regions'],
#			  -default=>$default{freq_estimate});
#print "<BR>";

&ReturnTable();

print "<HR width=550 align=left>\n";

### send results by email or display on the browser
&SelectOutput('email');

### action buttons
print "<ul><ul><table class = 'formbutton'>\n";
print "<tr valign=middle>\n";
print "<td>", $query->submit(-label=>"GO"), "</td>\n";
print "<td>", $query->reset(-id=>"reset"), "</td>\n";
print $query->end_form;

### data for the demo 
$demo_sequence = "";

open($fh, "demo_files/local-word-analysis_demo_seq.fa");
while(my $row = <$fh>){
    chomp $row;
    $demo_sequence .= $row;
    $demo_sequence .= "\\n";
}
print '<script>
function setDemo(demo_sequence){
    $("#reset").trigger("click");
    lth_occ.value = 2;
    window_width.value = 800;
    $("#bg_level_organism").prop("checked",true);
    bg_window_width.value = 800;
    sequence.value = demo_sequence;
    $("#organism").val("Saccharomyces_cerevisiae");
    $("#organism_name").val("Saccharomyces cerevisiae");
}
</script>';

print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_sequence'" .')">DEMO</button>';
print "</B></TD>\n";



print "<TD><B><A class='iframe' HREF='help.local-word-analysis.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A class='iframe' HREF='tutorials/tut_local-word-analysis.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "<hr />";

print $query->end_html;

exit(0);


################################################################
## Table with all the supported statistics and thresholds
sub ReturnTable {

print '<b><a class="iframe" href="help.local-word-analysis.html#thresholds">Thresholds</a><br />', "\n";

#print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>CENTER},
			 [
			  $query->th(["Fields ",
				   " <a class='iframe' href='help.local-word-analysis.html#thresholds'>Lower<br>Thresholds</a>\n",
				   " <a class='iframe' href='help.local-word-analysis.html#thresholds'>Upper<br>Thresholds</a>\n"]),

			  ### occurrences
			  $query->td(['Occurrences',
				      $query->textfield(-name=>'lth_occ',-id=>'lth_occ',
							-default=>$default{lth_occ},
							-size=>5),
				      $query->textfield(-name=>'uth_occ',
							-default=>$default{uth_occ},
							-size=>5)
				   ]),

			  ### binomial proba
#			  $query->td(["Probability",
#				      $query->textfield(-name=>'lth_occ_P',
#							-default=>$default{lth_occ_P},
#							-size=>5),
#				      $query->textfield(-name=>'uth_occ_P',
#							-default=>$default{uth_occ_P},
#							-size=>5),
#				      $query->popup_menu(-name=>'side',
#							 -Values=>['over-represented','under-represented','both'],
#							 -default=>$default{side})
#				     ]),

			  ### binomial E-value
			  $query->td(["E-value",
#				      $query->checkbox(-name=>'proba',
#						       -checked=>$default{proba},
#						       -label=>' Binomial E-value '),
				      $query->textfield(-name=>'lth_occ_E',
							-default=>$default{lth_occ_E},
							-size=>5),
				      $query->textfield(-name=>'uth_occ_E',
							-default=>$default{uth_occ_E},
							-size=>5),
				     ]),

			  ### significance index
			  $query->td(["Significance",
#				      $query->checkbox(-name=>'proba',
#						    -checked=>$default{proba},
#						    -label=>' Significance '),
				   $query->textfield(-name=>'lth_occ_sig',
						     -default=>$default{lth_occ_sig},
						     -size=>5),
				   $query->textfield(-name=>'uth_occ_sig',
						     -default=>$default{uth_occ_sig},
						     -size=>5)
				   ]),


			  ### frequencies
#			  $query->td(["Frequencies",
#				   $query->textfield(-name=>'lth_observed_freq',
#						     -default=>$default{lth_observed_freq},
#						     -size=>5),
#				   $query->textfield(-name=>'uth_observed_freq',
#						     -default=>$default{uth_observed_freq},
#						     -size=>5)
#				      ]),


			  ### rank
			  $query->td(["Rank",
				      $query->textfield(-name=>'lth_rank',
							-default=>$default{lth_rank},
							-size=>5),
				      $query->textfield(-name=>'uth_rank',
							-default=>$default{uth_rank},
							-size=>5)
				      ]),
			  ### w_rank
			  $query->td(["Window rank",
				      $query->textfield(-name=>'lth_w_rank',
							-default=>$default{lth_w_rank},
							-size=>5),
				      $query->textfield(-name=>'uth_w_rank',
							-default=>$default{uth_w_rank},
							-size=>5)
				      ]),



			  
			  ]
			       )
		    );
#print "</BLOCKQUOTE>\n";

}
