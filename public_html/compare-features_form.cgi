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


#local @supported_input_formats = sort(keys(%RSAT::feature::supported_input_format));
local @supported_input_formats = qw(bed dnapat ft galaxy_seq gft gff gff3bed swembl ucsc_seq);

### default values for filling the form
$default{featQ} = "";
$default{upload_query_features} = "";
$default{featRef} = "";
$default{upload_ref_features} = "";

$default{stats} = "checked";
$default{diff} = "";
$default{inter} = "on";
$default{inter_len} = "1";
$default{inter_cov} = "none";

$default{input_format} = "bed";


### print the form ###
&RSA_header("compare-features", 'form');
print "<CENTER>";
print "Compare two or more sets of features. The web-based program takes as input two feature files and computes the intersection, union and difference between features, as well as  contingency tables and comparison statistics. Note: the command-line version of this tool can take more than 2 files.<P>\n";
print "Program developed by <A HREF='mailto:jvhelden\@ulb.ac.be (Jacques van Helden)'>Jacques van Helden</A>\n";
print "and <A HREF='mailto:jturatsi\@ulb.ac.be (Jean-Valéry Turatsinze)'>Jean-Valéry Turatsinze</A>\n";
print "</CENTER>";
print "<hr>";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
}

print $query->start_multipart_form(-action=>"compare-features.cgi");

#### features text areas

print ("<table border='0' cellspacing='0' cellpadding='4'>\n");

print ("<tr align='center'><th>Query features</th><th>Reference features</th></tr>\n");

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
print ("<tr align='center'></td>");


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

#### feature format (pop-up menu)
print "<A HREF='help.convert-features.html'><B>Input feature format</B></a>&nbsp;";
print  $query->popup_menu(-name=>'feature_format',
			 -values=>[@supported_input_formats],
			 -default=>$default{input_format});
print "<br/>";


#print "<hr color=\"#BBBBBB\"><p>";

#### table with all the statistics and thresholds


print "<p><a href='help.compare-features.html#return_fields'><b>Output fields</b></a>\n", "&nbsp"x5;

## Return matching statistics
print $query->checkbox(-name=>'stats', 
		       -checked=>$default{stats},
		       -label=>' Statistics ');

### Return intersections
print $query->checkbox(-name=>'inter',
		       -checked=>$default{inter},
		       -label=>' Intersections ');

### Return differences
print $query->checkbox(-name=>'diff',
		       -checked=>$default{diff},
		       -label=>' Differences ');

print "</p>";

print "<p><b>Thresholds</b></p>\n";
print $query->table({-border=>0,-cellpadding=>3,-cellspacing=>3},
		    $query->Tr({-align=>'left',-valign=>'middle'},
			       [
				$query->th([" <A HREF='help.compare-features.html#return_fields'>Field</A> ",
					    " <A HREF='help.compare-features.html#thresholds'>Lower threshold</A> ",
					    ]),
				
				### Query class size
				$query->td(['Minimal overlap (bp)',
					    $query->textfield(-name=>'inter_len',
							      -default=>$default{inter_len},
							      -size=>5),
					    ]),
				### Reference class size
				$query->td([' Intersection coverage (0-1)',
					    $query->textfield(-name=>'inter_cov',
							      -default=>$default{inter_cov},
							      -size=>5),
					    ]),

			 ]
			)
		);

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
print "<tr valign=middle>\n";
#print "<TD>", $query->submit(-label=>"DEMO"), "</TD>\n";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print "<TD>", $query->reset, "</TD>\n";
print $query->end_form;


################################################################
### data for the demo 
print $query->start_multipart_form(-action=>"compare-features_form.cgi");
$demo_query_file = "demo_files/Blow_2010_GSM348064_forebrain_p300_peaks.bed";
$demoQuery=`cat $demo_query_file`;
# $demoQuery="
# eve	CRER	crer_1	DR	-5395	-5283	2	1.618	2.4e-02	2.0e-09	16.3	113
# eve	CRER	crer_2	DR	-4379	-4202	2	1.223	6.0e-02	6.3e-09	14.7	178
# eve	CRER	crer_3	DR	-3875	-3727	2	1.373	4.2e-02	5.9e-10	14.3	149
# eve	CRER	crer_4	DR	-3741	-3703	2	2.763	1.7e-03	3.1e-10	17.4	39
# eve	CRER	crer_5	DR	-3875	-3703	3	2.157	7.0e-03	2.4e-15	24.4	173
# eve	CRER	crer_6	DR	-3717	-3688	2	3.145	7.2e-04	3.7e-10	17.8	30
# eve	CRER	crer_7	DR	-3741	-3688	3	3.864	1.4e-04	2.8e-14	25.1	54
# eve	CRER	crer_8	DR	-3875	-3688	4	3.042	9.1e-04	2.2e-19	32.1	188
# eve	CRER	crer_9	DR	-3702	-3636	2	2.129	7.4e-03	1.4e-09	14.6	67
# eve	CRER	crer_10	DR	-3717	-3636	3	3.194	6.4e-04	5.5e-15	24.7	82
# eve	CRER	crer_11	DR	-3741	-3636	4	4.082	8.3e-05	4.2e-19	32	106
# eve	CRER	crer_12	DR	-3650	-3582	2	2.098	8.0e-03	1.1e-10	16.7	69
# eve	CRER	crer_13	DR	-3702	-3582	3	2.633	2.3e-03	9.6e-15	24.4	121
# eve	CRER	crer_14	DR	-3717	-3582	4	3.616	2.4e-04	3.9e-20	34.5	136
# eve	CRER	crer_15	DR	-3741	-3582	5	4.484	3.3e-05	3.0e-24	41.8	160
# eve	CRER	crer_16	DR	-3596	-3495	2	1.712	1.9e-02	5.5e-10	15.9	102
# eve	CRER	crer_17	DR	-3650	-3495	3	2.291	5.1e-03	8.2e-15	22.8	156
# eve	CRER	crer_18	DR	-1377	-1199	2	1.218	6.0e-02	2.9e-09	15	179
# eve	CRER	crer_19	DR	-1213	-1194	2	3.991	1.0e-04	2.7e-09	15.1	20
# eve	CRER	crer_20	DR	-1377	-1194	3	2.078	8.4e-03	2.2e-13	22.3	184
# eve	CRER	crer_21	DR	-1208	-1073	2	1.452	3.5e-02	3.4e-09	15.4	136
# eve	CRER	crer_22	DR	-1213	-1073	3	2.425	3.8e-03	1.2e-13	23.2	141
# eve	CRER	crer_23	DR	-1087	-945	2	1.408	3.9e-02	9.9e-10	16.7	143
# eve	CRER	crer_24	DR	-959	-780	2	1.214	6.1e-02	1.2e-09	16.4	180
# eve	CRER	crer_25	DR	-794	-687	2	1.659	2.2e-02	4.6e-09	15.2	108
# eve	CRER	crer_26	DR	-701	-676	2	3.392	4.1e-04	8.1e-09	14.6	26
# eve	CRER	crer_27	DR	-794	-676	3	2.656	2.2e-03	4.6e-13	22.4	119
# eve	CRER	crer_28	DR	-369	-341	2	3.200	6.3e-04	4.8e-10	17.2	29
# eve	CRER	crer_29	DR	-355	-299	2	2.304	5.0e-03	4.8e-10	17.2	57
# eve	CRER	crer_30	DR	-369	-299	3	3.415	3.8e-04	1.1e-14	25.8	71
# eve	CRER	crer_31	DR	-313	-276	2	2.798	1.6e-03	3.1e-11	19.5	38
# eve	CRER	crer_32	DR	-355	-276	3	3.231	5.9e-04	6.8e-16	28.1	80
# eve	CRER	crer_33	DR	-369	-276	4	4.315	4.8e-05	1.5e-20	36.7	94
# eve	CRER	crer_34	DR	-290	-220	2	2.069	8.5e-03	5.0e-14	23.6	71
# eve	CRER	crer_35	DR	-313	-220	3	2.991	1.0e-03	1.1e-18	32.2	94
# eve	CRER	crer_36	DR	-355	-220	4	3.616	2.4e-04	2.4e-23	40.8	136
# eve	CRER	crer_37	DR	-369	-220	5	4.629	2.3e-05	5.4e-28	49.4	150
# ";

$demo_ref_file = "demo_files/Blow_2010_GSM559653_midbrain_p300_peaks.bed";
$demoRef = `cat $demo_ref_file`;
# $demoRef="
# eve	TFBS	eve_mas	chr2R	-5303	-5198	5	105
# eve	TFBS	eve_proximal_promoter_inc._TATA	chr2R	-127	80	9	207
# eve	TFBS	eve_stripe2	chr2R	-1480	-996	19	484
# eve	TFBS	eve_stripe_3+7	chr2R	-3741	-3230	18	511
# ";

print "<TD><B>";
print $query->hidden(-name=>'featQ',-default=>$demoRef);
print $query->hidden(-name=>'featRef',-default=>$demoQuery);
print $query->hidden(-name=>'inter',-default=>"on");
print $query->hidden(-name=>'stats',-default=>"on");
print $query->hidden(-name=>'inter_cov',-default=>"0.25");
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
print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>\n";

print "</font>\n";
print "<hr>";

print $query->end_html;

exit(0);



