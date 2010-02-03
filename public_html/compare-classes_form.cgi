#!/usr/bin/perl
#### this cgi script fills the HTML form for the program compare-classes
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
$default{query_classes} = "";
$default{upload_query_classes} = "";
$default{ref_classes} = "";
$default{upload_ref_classes} = "";


$default{occ} = "checked";
$default{lth_occ} = 1;
$default{uth_occ} = "none";
$default{freq} = "checked";
$default{sig} = "checked";
$default{lth_sig} = 0;
$default{uth_sig} = "none";
$default{proba} = "checked";
$default{freq} = "checked";
$default{jac} = "checked";
$default{entropy} = "checked";
$default{members} = "";
$default{sort_key} = "sig";
$default{pop_size} = "auto";
$default{entropy} = "checked";
$default{jac} = "checked";

### print the form ###
&RSA_header("compare-classes", 'form');
print "<CENTER>";
print "Compare two classifications (clustering results, functional classes, ...), and assess the statistical significance of common members between each pair of classes.<P>\n";
print "Program developed by <A HREF='mailto:jtran\@bigre.ulb.ac.be (Joseph Tran)'>Joseph Tran</A>\n";
print "and <A HREF='mailto:jvanheld\@bigre.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A>\n";
print "</CENTER>";
print "<HR>";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
} 

print $query->start_multipart_form(-action=>"compare-classes.cgi");


#### upload query classifcation file
print "<a href='help.compare-classes.html#upload_query_classes'>Query classification file</a><BR>";
print $query->filefield(-name=>'upload_query_classes',
			-default=>$default{upload_query_classes},
			-size=>30,
			-maxlength=>200);
print "<p>";

#### upload reference classifcation file
print "<a href='help.compare-classes.html#upload_ref_classes'>Reference classification file</a><BR>";
print $query->filefield(-name=>'upload_ref_classes',
			-default=>$default{upload_ref_classes},
			-size=>30,
			-maxlength=>200);
print "<p>";

#### table with all the statistics and thresholds
print "<h4>Return</h4>\n";

print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-classes.html#return_fields'>Fields</A> "]),

				### occurrences
				$query->td([$query->checkbox(-name=>'occ',
							     -checked=>$default{occ},
							     -label=>' Occurrences ')
					    ]),

				### Frequencies
				$query->td([$query->checkbox(-name=>'freq',
							     -checked=>$default{freq},
							     -label=>' Frequencies ')
					    ]),

				### Probabilities
				$query->td([$query->checkbox(-name=>'proba',
							     -checked=>$default{proba},
							     -label=>' Probabilities ')
					    ]),


				### Jaccard index
				$query->td([$query->checkbox(-name=>'jac',
							     -checked=>$default{jac},
							     -label=>' Jaccard index ')
					    ]),
				### Entropy
				$query->td([$query->checkbox(-name=>'entropy',
							     -checked=>$default{entropy},
							     -label=>' Entropy ')
					    ]),

				### Members
				$query->td([$query->checkbox(-name=>'members',
							     -checked=>$default{members},
							     -label=>' Members '),
					    ]),

			 ]
			)
		);

print "<h4>Thresholds</h4>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-classes.html#return_fields'>Fields</A> ",
					    " <A HREF='help.compare-classes.html#thresholds'>Lower<BR>Threshold</A> ",
					    " <A HREF='help.compare-classes.html#thresholds'>Upper<BR>Threshold</A> ",
					    ]),
				
				### Query class size
				$query->td([' Query size ',
					    $query->textfield(-name=>'lth_q',
							      -default=>$default{lth_q},
							      -size=>5),
					    $query->textfield(-name=>'uth_q',
							      -default=>$default{uth_q},
							      -size=>5),
					    ]),
				### Reference class size
				$query->td([' Reference size ',
					    $query->textfield(-name=>'lth_r',
							      -default=>$default{lth_r},
							      -size=>5),
					    $query->textfield(-name=>'uth_r',
							      -default=>$default{uth_r},
							      -size=>5),
					    ]),
				### Intersection size
				$query->td([' Intersection size ',
					    $query->textfield(-name=>'lth_qr',
							      -default=>$default{lth_qr},
							      -size=>5),
					    $query->textfield(-name=>'uth_qr',
							      -default=>$default{uth_qr},
							      -size=>5),
					    ]),
				### Significance 
				$query->td([' Significance ',
					    $query->textfield(-name=>'lth_sig',
							      -default=>$default{lth_sig},
							      -size=>5),
					    $query->textfield(-name=>'uth_sig',
							      -default=>$default{uth_sig},
							      -size=>5),
					    ]),

				### P-value 
				$query->td([' P-value ',
					    $query->textfield(-name=>'lth_pval',
							      -default=>$default{lth_pval},
							      -size=>5),
					    $query->textfield(-name=>'uth_pval',
							      -default=>$default{uth_pval},
							      -size=>5),
					    ]),

				### E-value 
				$query->td([' E-value ',
					    $query->textfield(-name=>'lth_eval',
							      -default=>$default{lth_eval},
							      -size=>5),
					    $query->textfield(-name=>'uth_eval',
							      -default=>$default{uth_eval},
							      -size=>5),
					    ]),
				### Jaccard index
				$query->td([' Jaccard index ',
					    $query->textfield(-name=>'lth_jac',
							      -default=>$default{lth_jac},
							      -size=>5),
					    $query->textfield(-name=>'uth_jac',
							      -default=>$default{uth_jac},
							      -size=>5),
					    ]),
				$query->td([' Mutual information ',
					    $query->textfield(-name=>'lth_mi',
							      -default=>$default{lth_mi},
							      -size=>5),
					    $query->textfield(-name=>'uth_mi',
							      -default=>$default{uth_mi},
							      -size=>5),
					    ]),

			 ]
			)
		);



################################################################
## sort key
print "<b><a href='help.compare-classes.html#sort_key'>Sort key </a></b>";
print  $query->popup_menu(-name=>'sort_key',
			  -Values=>['sig',
				    'E_val', 
				    'P_val',
				    'Jaccard index',
				    'Mutual information',
				    'names'
				    ],
			  -default=>$sequence_format);

################################################################
## population size
print "&nbsp"x8, "<b><a href='help.compare-classes.html#pop_size'>Population size </a></b>";
print $query->textfield(-name=>'pop_size',
			-default=>$default{pop_size},
			-size=>5);

### send results by email or display on the browser
print "<HR width=550 align=left>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;

### data for the demo 
# print $query->start_multipart_form(-action=>"compare-classes_form.cgi");
# print "<TD><B>";
# print $query->hidden(-name=>'sort_key',-default=>"sig");
# print $query->submit(-label=>"DEMO");
# print "</B></TD>\n";
# print $query->end_form;


print "<TD><B><A HREF='help.compare-classes.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_compare-classes.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "<HR>";

print $query->end_html;

exit(0);


