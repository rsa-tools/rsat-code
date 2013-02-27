#!/usr/bin/perl
#### this cgi script fills the HTML form for the program oligo-diffi
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
$default{test_seq} = "";
$default{upload_test_seq} = "";
$default{ctrl_seq} = "";
$default{upload_ctrl_seq} = "";
$default{lth_occ} = "3";
$default{uth_occ} = "none";
$default{lth_occ_sig} = "0";
$default{uth_occ_sig} = "none";
$default{lth_occ_Pval} = "none";
$default{uth_occ_Pval} = "none";
$default{lth_occ_Eval} = "none";
$default{uth_occ_Eval} = "none";
$default{lth_ratio} = "none";
$default{uth_ratio} = "none";
$default{oligo_len} = "6";
$default{noov} = "checked";
$default{strand} = "both strands";
$default{purge} = 'checked';

&MatrixFromPatterns_defaults();

### print the form ###
&RSA_header("oligo-diff", 'form');
print "<CENTER>";
print "Compare oligonucleotide occurrences between two input sequence files, and return oligos that are significantly enriched in one of the files respective to the other one.<P>\n";
print "Program developed by <A HREF='mailto:jvhelden\@ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A>\n";
print "</CENTER>";
print "<HR>";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

print $query->start_multipart_form(-action=>"oligo-diff.cgi");

################################################################
## Sequence sets
print "<h2>Input sequences</h2>\n";

################################################################
## Test sequence set

## Textarea for the test sequence set
print "<p><b><a href='help.oligo-diff.html#upload_test_seq'>Test sequence set</a></b><br>";
print "Paste your sequence in fasta format in the box below<BR>\n";
print $query->textarea(-name=>'test_seq',
		       -default=>$default{test_seq},
		       -rows=>8,
		       -columns=>80);

#### Upload test sequence file
print  "<br>Upload file","&nbsp;"x3;
print $query->filefield(-name=>'upload_test_seq',
			-default=>$default{upload_test_seq},
			-size=>60,
			-maxlength=>200);
print "</p>\n";

################################################################
## Control sequence set
print "<p><b><a href='help.oligo-diff.html#upload_ctrl_seq'>Control sequence set</a></b><br>";
print "Paste your sequence in fasta format in the box below<BR>\n";


## Textarea for the control sequence set
print $query->textarea(-name=>'ctrl_seq',
		       -default=>$default{ctrl_seq},
		       -rows=>8,
		       -columns=>80);

#### Upload secnd sequence file
print  "<br>Upload file","&nbsp;"x3;
print $query->filefield(-name=>'upload_ctrl_seq',
			-default=>$default{upload_ctrl_seq},
			-size=>60,
			-maxlength=>200);
print "</p>";


#### purge sequences
print $query->checkbox(-name=>'purge',
		       -checked=>$default{purge},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-diff.html#purge'><B>purge sequences (highly recommended)</B></A>";
print "<BR>";


################################################################
## Oligonucleotide counting options
print "<hr>\n";
print "<h2>Oligonucleotide countint options</h2>\n";

## oligo size
print "<B><A HREF='help.oligo-diff.html#oligo_len'>Oligomer length</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'oligo_len',
			 -Values=>[1,2,3,4,5,6,7,8],
			 -default=>$default{oligo_len});

## prevent overlapping matches of the same pattern
#print "<br>\n";
print "&nbsp;"x5;
print $query->checkbox(-name=>'noov',
		       -checked=>$default{noov},
		       -label=>'');
print "&nbsp;<A HREF='help.oligo-diff.html#noov'><B>prevent overlapping matches</B></A>";

## strand
print "<br>\n";
#print "&nbsp;"x5;
print "<B><A HREF='help.oligo-diff.html#count_strands'>Count on</A>&nbsp;</B>\n";
print $query->popup_menu(-name=>'strand',
			 -Values=>['single strand',
				  'both strands'],
			 -default=>$default{strand});
print "<hr>\n";



################################################################
## table with all the thresholds
print "<h2>Thresholds</h2>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.oligo-diff.html#return_fields'>Fields</A> ",
					    " <A HREF='help.oligo-diff.html#thresholds'>Lower<BR>Threshold</A> ",
					    " <A HREF='help.oligo-diff.html#thresholds'>Upper<BR>Threshold</A> ",
					    ]),

				### Query class size
				$query->td([' Occurrences ',
					    $query->textfield(-name=>'lth_occ',
							      -default=>$default{lth_occ},
							      -size=>5),
					    $query->textfield(-name=>'uth_occ',
							      -default=>$default{uth_occ},
							      -size=>5),
					    ]),

				### Significance 
				$query->td([' Significance ',
					    $query->textfield(-name=>'lth_occ_sig',
							      -default=>$default{lth_occ_sig},
							      -size=>5),
					    $query->textfield(-name=>'uth_occ_sig',
							      -default=>$default{uth_occ_sig},
							      -size=>5),
					    ]),

				### P-value 
				$query->td([' P-value ',
					    $query->textfield(-name=>'lth_occ_Pval',
							      -default=>$default{lth_occ_Pval},
							      -size=>5),
					    $query->textfield(-name=>'uth_occ_Pval',
							      -default=>$default{uth_occ_Pval},
							      -size=>5),
					    ]),

				### E-value 
				$query->td([' E-value ',
					    $query->textfield(-name=>'lth_occ_Eval',
							      -default=>$default{lth_occ_Eval},
							      -size=>5),
					    $query->textfield(-name=>'uth_occ_Eval',
							      -default=>$default{uth_occ_Eval},
							      -size=>5),
					    ]),
# 				### Jaccard index
# 				$query->td([' Jaccard index ',
# 					    $query->textfield(-name=>'lth_jac',
# 							      -default=>$default{lth_jac},
# 							      -size=>5),
# 					    $query->textfield(-name=>'uth_jac',
# 							      -default=>$default{uth_jac},
# 							      -size=>5),
# 					    ]),

			 ]
			)
		);

#### Convert patterns to matrix
&MatrixFromPatterns_print_form();

### send results by email or display on the browser
print "<HR width=550 align=left>\n";
&SelectOutput();

### action buttons
print "<ul><ul><table class='formbutton'>\n";
print "<tr valign=middle>\n";
#print "<td>", $query->submit(-label=>"DEMO"), "</td>\n";
print "<td>", $query->submit(-label=>"GO"), "</td>\n";
print "<td>", $query->reset, "</td>\n";
print $query->end_form;


################################################################
### data for the demo
print $query->start_multipart_form(-action=>"oligo-diff_form.cgi");

$demo_test = $ENV{RSAT}."/public_html/demo_files/MET_up800-noorf.fasta";
#$demo_test = $ENV{RSAT}."/public_html/demo_files/peak-motifs_GSM348066_limb_p300_1000peaks.fa";
$demo_ctrl = $ENV{RSAT}."/public_html/demo_files/PHO_up800-noorf.fasta";
#$demo_ctrl = $ENV{RSAT}."/public_html/demo_files/peak-motifs_GSM559652_heart_p300_1000peaks.fa";
$demo_test_seq=`cat $demo_test`;
$demo_ctrl_seq=`cat $demo_ctrl`;

print "<TD><B>";
print $query->hidden(-name=>'test_seq',-default=>$demo_ctrl_seq);
print $query->hidden(-name=>'ctrl_seq',-default=>$demo_test_seq);
print $query->hidden(-name=>'side',-default=>'both');
#print $query->hidden(-name=>'to_matrix',-default=>'0');
print $query->hidden(-name=>'ratio',-default=>'none');
print $query->hidden(-name=>'oligo_len',-default=>'5');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


print "<TD><B><A HREF='help.oligo-diff.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_oligo-diff.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);



