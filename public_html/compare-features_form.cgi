#!/usr/bin/perl
#### this cgi script fills the HTML form for the program compare-features
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
$default{featQ} = "";
$default{upload_query_features} = "";
$default{featRef} = "";
$default{upload_ref_features} = "";


$default{stats} = "checked";
$default{diff} = "checked";
$default{inter} = "checked";
$default{inter_len} = "none";
$default{inter_cov} = 0.8;

### print the form ###
&RSA_header("compare-features", 'form');
print "<CENTER>";
print "Compare two or more sets of features. This program takes as input several feature files (two or more), and calculates the intersection, union and difference between features. It also computes contingency tables and comparison statistics.<P>\n";
#print "Program developed by <A HREF='mailto:jtran\@scmbb.ulb.ac.be (Joseph Tran)'>Joseph Tran</A>\n";
#print "and <A HREF='mailto:jvanheld\@scmbb.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A>\n";
print "</CENTER>";
print "<HR>";
print "<blockquote>";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

print $query->start_multipart_form(-action=>"compare-features.cgi");

#### features text areas

print ("<table border=0 cellspacing=0 cellpadding=4 align=center>\n");

print ("<tr align = 'center'><th>Query features</th><th>Reference features</th></tr>\n");

################################################################
## Query input classes textarea
print "<tr><td>\n";
print $query->textarea(-name=>'featQ',
		       -default=>$default{featQ},
		       -rows=>10,
		       -columns=>60);


#print ("<textarea name='featQ' rows='10' cols='60'>$demo_featQ</textarea>");
print "</td><td>";

################################################################
## Reference input classes textarea

print $query->textarea(-name=>'featRef',
		       -default=>$default{featRef},
		       -rows=>10,
		       -columns=>60);

#print ("<td><textarea name='featRef' rows='10' cols='60'>$demo_featRef</textarea></td></tr>");
print ("<tr align = 'center'></td>");


#### upload query classifcation file
print ("<tr><td>");
print "<a href='help.compare-features.html#upload_query_features'>Query feature file</a><BR>";
print $query->filefield(-name=>'upload_query_features',
			-default=>$default{upload_query_features},
			-size=>30,
			-maxlength=>200);
print ("</td>");
#print "<p>";
print ("<td>");
#### upload reference feature file
print "<a href='help.compare-features.html#upload_ref_features'>Reference feature file</a><BR>";
print $query->filefield(-name=>'upload_ref_features',
			-default=>$default{upload_ref_features},
			-size=>30,
			-maxlength=>200);
print ("</td></tr>");
print ("</table>");
print "<p>";

#### table with all the statistics and thresholds

print ("<table border=0 cellspacing=0 cellpadding=4 align=center>\n");

print ("<tr align = 'center'><td>\n");

print "<h4>Return</h4>\n";

print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-features.html#return_fields'>Fields</A> "]),

				### occurrences
				$query->td([$query->checkbox(-name=>'stats',
							     -checked=>$default{stats},
							     -label=>' Statistics ')
					    ]),

				### Frequencies
				$query->td([$query->checkbox(-name=>'inter',
							     -checked=>$default{inter},
							     -label=>' Intersections ')
					    ]),

				### Probabilities
				$query->td([$query->checkbox(-name=>'diff',
							     -checked=>$default{diff},
							     -label=>' Differences ')
					    ]),

			 ]
			)
		);
print "</BLOCKQUOTE>\n";
print ("</td><td>");
print "<h4>Thresholds</h4>\n";
print "<BLOCKQUOTE>\n";
print $query->table({-border=>0,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-features.html#return_fields'>Fields</A> ",
					    " <A HREF='help.compare-features.html#thresholds'>Lower<BR>Threshold</A> ",
#					    " <A HREF='help.compare-features.html#thresholds'>Upper<BR>Threshold</A> ",
					    ]),
				
				### Query class size
				$query->td([' Intersection size (nt) ',
					    $query->textfield(-name=>'inter_len',
							      -default=>$default{inter_len},
							      -size=>5),
#					    $query->textfield(-name=>'uth_q',
#							      -default=>$default{uth_q},
#							      -size=>5),
					    ]),
				### Reference class size
				$query->td([' Intersection coverage (0-1)',
					    $query->textfield(-name=>'inter_cov',
							      -default=>$default{inter_cov},
							      -size=>5),
#					    $query->textfield(-name=>'uth_r',
#							      -default=>$default{uth_r},
#							      -size=>5),
					    ]),

			 ]
			)
		);
print "</BLOCKQUOTE>\n";
print ("</td></tr>");
print ("</table>");
print "<p>";



################################################################
## output format
# print "<b><a href='help.compare-features.html#sort_key'>Sort key </a></b>";
# print  $query->popup_menu(-name=>'sort_key',
# 			  -Values=>['sig',
# 				    'E_val', 
# 				    'P_val',
# 				    'Jaccard index',
# 				    'Mutual information',
# 				    'names'
# 				    ],
# 			  -default=>$sequence_format);

################################################################
## population size
# print "&nbsp"x8, "<b><a href='help.compare-features.html#pop_size'>Population size </a></b>";
# print $query->textfield(-name=>'pop_size',
# 			-default=>$default{pop_size},
# 			-size=>5);

### send results by email or display on the browser
print "<HR width=550 align=left>\n";
&SelectOutput();

### action buttons
print "<UL><UL><TABLE class='formbutton'>\n";
print "<TR VALIGN=MIDDLE>\n";
#print "<TD>", $query->submit(-label=>"DEMO"), "</TD>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"compare-features_form.cgi");
$demoQuery="eve	TFBS	eve_mas	chr2R	-5303	-5198	5	105
eve	TFBS	eve_proximal_promoter_inc._TATA	chr2R	-127	80	9	207
eve	TFBS	eve_stripe2	chr2R	-1480	-996	19	484
eve	TFBS	eve_stripe_3+7	chr2R	-3741	-3230	18	511
";

$demoRef="eve	TFBS	eve_mas	chr2R	-5303	-5198	5	105
eve	TFBS	eve_proximal_promoter_inc._TATA	chr2R	-127	80	9	207
eve	TFBS	eve_stripe2	chr2R	-1480	-996	19	484
eve	TFBS	eve_stripe_3+7	chr2R	-3741	-3230	18	511
";

print "<TD><B>";
print $query->hidden(-name=>'featQ',-default=>$demoQuery);
print $query->hidden(-name=>'featRef',-default=>$demoRef);
# print $query->hidden(-name=>'input_format',-default=>'tab');
# print $query->hidden(-name=>'info',-default=>"on");
# print $query->hidden(-name=>'weights',-default=>"on");
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;


### data for the demo 
# print $query->start_multipart_form(-action=>"compare-features_form.cgi");
# print "<TD><B>";
# print $query->hidden(-name=>'sort_key',-default=>"sig");
# print $query->submit(-label=>"DEMO");
# print "</B></TD>\n";
# print $query->end_form;


print "<TD><B><A HREF='help.compare-features.html'>MANUAL</A></B></TD>\n";
#print "<TD><B><A HREF='tutorials/tut_compare-features.html'>TUTORIAL</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);



