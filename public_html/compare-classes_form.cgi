#!/usr/bin/perl
#### this cgi script fills the HTML form for the program compare-classes
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
$default{query_classes} = "";
$default{upload_query_classes} = "";
$default{ref_classes} = "";
$default{upload_ref_classes} = "";


$default{occ} = "checked";
$default{lth_occ} = 1;
$default{uth_occ} = "none";
$default{percent} = "checked";
$default{lth_percent} = "none";
$default{uth_percent} = "none";
$default{sig} = "checked";
$default{lth_sig} = 0;
$default{uth_sig} = "none";
$default{proba} = "checked";
$default{lth_proba} = "none";
$default{uth_proba} = "none";
$default{members} = "";
$default{sort_key} = "occ";
$default{pop_size} = "auto";

### print the form ###
&RSA_header("compare-classes");
print "<CENTER>";
print "Compare two classifications (clustering results, functional classes, ...), and assess the statistical significance of common members between each pair of classes.<P>\n";
print "Program developed by <A HREF='mailto:jtran\@scmbb.ulb.ac.be (Joseph Tran)'>Joseph Tran</A>\n";
print "and <A HREF='mailto:jvanheld\@scmbb.ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A>\n";
print "</CENTER>";
print "<HR>";
print "<blockquote>";

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
print "<h3>Return</h3>\n";

print "<BLOCKQUOTE>\n";
print $query->table({-border=>1,-cellpadding=>0,-cellspacing=>0},
		    $query->Tr({-align=>left,-valign=>TOP},
			       [
				$query->th([" <A HREF='help.compare-classes.html#return_fields'>Fields</A> ",
					    " <A HREF='help.compare-classes.html#thresholds'>Lower<BR>Threshold</A> ",
					    " <A HREF='help.compare-classes.html#thresholds'>Upper<BR>Threshold</A> ",
					    ]),
				
				### occurrences
				$query->td([$query->checkbox(-name=>'occ',
							     -checked=>$default{occ},
							     -label=>' Occurrences '),
					    $query->textfield(-name=>'lth_occ',
							      -default=>$default{lth_occ},
							      -size=>5),
					    $query->textfield(-name=>'uth_occ',
							      -default=>$default{uth_occ},
							      -size=>5),
					    ]),

				### percentages
				$query->td([$query->checkbox(-name=>'percent',
							     -checked=>$default{percent},
							     -label=>' Percentages '),
					    $query->textfield(-name=>'lth_percent',
							      -default=>$default{lth_percent},
							      -size=>5),
					    $query->textfield(-name=>'uth_percent',
							      -default=>$default{uth_percent},
							      -size=>5),
					    ]),

				### Probabilities
				$query->td([$query->checkbox(-name=>'proba',
							     -checked=>$default{proba},
							     -label=>' Probabilities '),
					    $query->textfield(-name=>'lth_proba',
							      -default=>$default{lth_proba},
							      -size=>5),
					    $query->textfield(-name=>'uth_proba',
							      -default=>$default{uth_proba},
							      -size=>5),
					    ]),

				### Significance
				$query->td([$query->checkbox(-name=>'Significance',
							     -checked=>$default{sig},
							     -label=>' Significance '),
					    $query->textfield(-name=>'lth_sig',
							      -default=>$default{lth_sig},
							      -size=>5),
					    $query->textfield(-name=>'uth_sig',
							      -default=>$default{uth_sig},
							      -size=>5),
					    ]),

				### Members
				$query->td([$query->checkbox(-name=>'members',
							     -checked=>$default{members},
							     -label=>' Members '),
					    '',
					    '',
					    ]),

			 ]
			)
		);
print "</BLOCKQUOTE>\n";


################################################################
## sort key
print "<b><a href='help.compare-classes.html#sort_key'>Sort key </a></b>";
print  $query->popup_menu(-name=>'sort_key',
			  -Values=>[
				    'occ',
				    'E_val', 
				    'P_val',
				    'sig',
				    'name',
				    'query_percent',
				    'ref_percent',
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
print "<UL><UL><TABLE>\n";
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
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</FONT>\n";
print "</blockquote>";
print "<HR>";

print $query->end_html;

exit(0);


