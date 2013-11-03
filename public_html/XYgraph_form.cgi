#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### intialization

### Read the CGI query
$query = new CGI;

### read values to fill the form ###
$default{title} = "";
$default{title2} = "";
$default{pointsize} = "2";
$default{bg} = "white";
$default{legend} = 'checked';
$default{htmap} = 'checked';
$default{label_col} = 1;
$default{symbols} = 'checked';
$default{format} = $ENV{rsat_img_format} || "png";

## X axis options
$default{xcol} = 1;
$default{xleg1} = "";
$default{xleg2} = "";
$default{xgstep1} = "auto";
$default{xgstep2} = "auto";
$default{xsize} = 500;
$default{xlog} = "none";
$default{xmax} = "auto";
$default{xmin} = "auto";

## Y axis options
$default{ycol} = 2;
$default{yleg1} = "";
$default{yleg2} = "";
$default{ygstep1} = "auto";
$default{ygstep2} = "auto";
$default{ysize} = 500;
$default{ylog} = "none";
$default{ymax} = "auto";
$default{ymin} = "auto";

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
    if ($query->param($key)) {
	$default{$key} = $query->param($key);
    }
} 


#### specific treatment for internal XYgraph file (piped from dna-patttern)
if (-e $query->param('XYgraph_file')) {
    $file = $query->param('XYgraph_file');
    $default{data} = `cat $file`;
} else {
    $default{data} = $query->param('data');
}
$default{data} =~ s/\"//g; #### remove quotes for security reasons (avoid imbedded command)


### remove quotes from the title and data
$default{data} =~ s/\"//g;
$default{title} =~ s/\"//g;


### print the form 
&RSA_header("XYgraph","form");


print "<CENTER>";
print "<P>Draws a XY graph from a table of numeric data\n";
print "</CENTER>";

print $query->start_multipart_form(-action=>"XYgraph.cgi");

## data
print "<B>Data</B><br>\n";
print $query->textarea(-name=>'data',
		       -default=>$default{data},
		       -rows=>6,
		       -columns=>60);
print "<br>File ";
print "</A></B>&nbsp;";
print $query->filefield(-name=>'uploaded_file',
			-size=>45);

print "<hr>\n";

print "<P><B>Graph options</B><br>\n";

## first title
print "<br>First Title", "&nbsp;"x2;
print $query->textfield(-name=>'title',
			-default=>$default{title},
			-size=>50,
			-maxlength=>80);

## second title
print "<br>Second Title", "&nbsp;"x2;
print $query->textfield(-name=>'title2',
			-default=>$default{title2},
			-size=>50,
			-maxlength=>80);

## point size
print "<br>Point size", "&nbsp;"x2;
print $query->textfield(-name=>'pointsize',
			-default=>$default{pointsize},
			-size=>50,
			-maxlength=>80);

## background color
print "<br>Background color";
print $query->popup_menu(-name=>'bg',
			 -Values=>["white","black","blue","gray"],
			 -default=>$default{bg});

## Image format
print "Format";
print $query->popup_menu(-name=>'format',
			 -Values=>["jpg","png","eps","pdf"],
			 -default=>$default{format});

print "<br>";
## draw lines to join points
print $query->checkbox(-name=>'lines',
		       -checked=>$default{lines},
		       -label=>' Lines between points ');

## legend
print $query->checkbox(-name=>'legend',
		       -checked=>$default{legend},
		       -label=>' First row contains legend');

print "<br>";
## symbols
print $query->checkbox(-name=>'symbols',
		       -checked=>$default{symbols},
		       -label=>' Use symbols');


## dynamic map
print $query->checkbox(-name=>'htmap',
		       -checked=>$default{htmap},
		       -label=>' Dynamic map (information in the status bar))');



################################################################
## X axis options
print "<hr>\n";
print "<br><B>X axis</B><br>\n";

## xcol
print "<br>Data column for X axis", "&nbsp;"x2;
print $query->textfield(-name=>'xcol',
			-default=>$default{xcol},
			-size=>5,
			);

## xsize
print "&nbsp;"x2, "Size (pixels)", "&nbsp;"x2;
print $query->textfield(-name=>'xsize',
			-default=>$default{xsize},
			-size=>5,
			);

## Logarithmic scale for the X axis
print "&nbsp;"x2, "Log base", "&nbsp;"x2;
print $query->textfield(-name=>'xlog',
			-default=>$default{xlog},
			-size=>4,
			);

## X labels
print "<br>First label", "&nbsp;"x2;
print $query->textfield(-name=>'xleg1',
			-default=>$default{xleg1},
			-size=>40,
			);
print "<br>Second label", "&nbsp;"x2;
print $query->textfield(-name=>'xleg2',
			-default=>$default{xleg2},
			-size=>40,
			);

## xmin
print "<br>Min value", "&nbsp;"x2;
print $query->textfield(-name=>'xmin',
			-default=>$default{xmin},
			-size=>5,
			);
## xmax
print "&nbsp;"x2, "Max value", "&nbsp;"x2;
print $query->textfield(-name=>'xmax',
			-default=>$default{xmax},
			-size=>5,
			);
## X grid
print "<br>First grid step", "&nbsp;"x2;
print $query->textfield(-name=>'xgstep1',
			-default=>$default{xgstep1},
			-size=>5,
			);
print "&nbsp;"x2,"Second grid step", "&nbsp;"x2;
print $query->textfield(-name=>'xgstep2',
			-default=>$default{xgstep2},
			-size=>5,
			);


################################################################
## Y axis options
print "<hr>\n";
print "<br><B>Y axis</B><br>\n";

## ycol
print "<br>Data columns for Y axis (e.g. 4,5,7)", "&nbsp;"x2;
print $query->textfield(-name=>'ycol',
			-default=>$default{ycol},
			-size=>5,
			);

## ysize
print "&nbsp;"x2, "Size (pixels)", "&nbsp;"x2;
print $query->textfield(-name=>'ysize',
			-default=>$default{ysize},
			-size=>5,
			);

## Logarithmic scale for the Y axis
print "&nbsp;"x2, "Log base", "&nbsp;"x2;
print $query->textfield(-name=>'ylog',
			-default=>$default{ylog},
			-size=>4,
			);

## Y labels
print "<br>First label", "&nbsp;"x2;
print $query->textfield(-name=>'yleg1',
			-default=>$default{yleg1},
			-size=>40,
			);
print "<br>Second label", "&nbsp;"x2;
print $query->textfield(-name=>'yleg2',
			-default=>$default{yleg2},
			-size=>40,
			);

## ymin
print "<br>Min value", "&nbsp;"x2;
print $query->textfield(-name=>'ymin',
			-default=>$default{ymin},
			-size=>5,
			);
## ymax
print "&nbsp;"x2, "Max value", "&nbsp;"x2;
print $query->textfield(-name=>'ymax',
			-default=>$default{ymax},
			-size=>5,
			);
## Y grid
print "<br>First grid step", "&nbsp;"x2;
print $query->textfield(-name=>'ygstep1',
			-default=>$default{ygstep1},
			-size=>5,
			);
print  "&nbsp;"x2,"Second grid step", "&nbsp;"x2;
print $query->textfield(-name=>'ygstep2',
			-default=>$default{ygstep2},
			-size=>5,
			);

print "<hr>\n";

print "<UL><UL><TABLE class = 'formbutton'><TR>";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print"<TD>", $query->reset, "</TD>\n";
print $query->end_form;


### data for the demo 
print $query->start_multipart_form(-action=>"XYgraph_form.cgi");
$demo_data = ";seq	markov_order_2	markov_order_3	markov_order_4	markov_order_5	chi2
tatata	54.09	21.00	15.44	4.06	889.136
atatat	48.37	16.99	12.67	1.68	604.254
tatgta	23.81	14.51	14.19	4.57	414.300
tacata	23.68	13.58	13.66	4.31	324.961
taaata	7.76	7.62	10.38	2.50	417.557
tgtata	25.32	10.95	7.65	3.28	387.945
tattta	5.25	6.39	6.92	2.93	242.951
aaataa	11.31	11.16	10.32	2.21	190.985
aataaa	7.96	7.82	8.81	2.87	206.891
atgtat	15.67	7.70	5.31	-0.30	225.927
taagta	5.51	5.29	8.99	3.73	194.527
acatat	17.40	6.12	5.75	0.67	222.117
catata	20.20	4.76	4.40	-0.46	261.653
atacat	14.60	6.14	6.57	-0.28	140.709
ttattt	6.33	7.62	8.01	0.26	145.632
tacgta	12.36	7.80	10.86	4.02	112.762
atatac	19.83	4.22	4.26	-0.02	228.029
atatgt	14.81	4.27	5.43	1.72	165.831
aaaaaa	24.86	13.92	7.47	2.19	137.726
tagata	1.40	2.58	7.27	2.02	108.966
ataaat	1.76	1.64	3.94	-1.60	163.664
tttatt	7.36	8.67	8.39	0.32	88.791
ttaaat	-6.29	2.60	5.93	2.78	116.107
tataca	14.22	3.22	3.66	0.50	142.831
tgtgta	9.48	4.64	4.78	2.77	173.857
attatt	-0.45	7.22	4.99	0.54	110.796
taatta	-5.36	7.80	3.89	0.82	156.250
gtatgt	11.95	6.37	7.98	1.84	63.260
atctat	3.65	3.78	6.65	1.21	86.762
aattaa	-7.36	5.40	5.20	4.15	113.195
ttaatt	-6.81	7.77	5.68	1.22	98.071
tttctt	20.93	10.42	6.28	3.64	126.385
tatcta	3.08	3.21	5.89	1.29	77.140
ttttct	22.37	8.68	5.63	3.04	111.101
acatac	8.82	3.39	4.76	0.54	63.318
gtatac	10.97	3.00	3.74	2.72	106.534
aaatag	2.62	3.67	5.18	1.53	83.667
tcattt	4.87	4.86	3.76	1.39	93.809
tagtta	-0.92	8.00	5.68	1.63	61.023
ttcatt	6.13	6.86	6.50	1.63	51.757
ttatta	-1.53	6.82	4.67	0.46	51.300
tgtaca	9.47	5.24	3.61	2.13	74.210
ttttgt	9.09	3.52	5.35	3.12	63.699
tagttt	0.44	5.54	3.82	1.78	81.174
ttattc	2.15	3.57	5.09	0.30	64.085
atgtgt	6.76	3.13	4.08	1.01	61.967
tattct	8.28	3.74	4.09	3.06	74.961
atagat	-0.37	0.74	3.95	0.10	56.987
gtaaat	-0.49	1.33	3.68	-0.50	64.505
aataga	1.29	1.96	3.55	0.95	76.647
tctatt	2.93	2.89	4.55	-1.04	53.530
ttacgt	3.96	3.70	4.00	0.28	51.598
taattg	-2.01	6.60	4.04	0.31	53.376
agttaa	-3.97	3.19	3.61	1.03	56.641
";
print "<TD><B>";

print $query->hidden(-name=>'data',-default=>$demo_data);
print $query->hidden(-name=>'title',-default=>'XYgraph demo');
print $query->hidden(-name=>'title2',-default=>'score comparisons');
print $query->hidden(-name=>'pointsize',-default=>'5');
print $query->hidden(-name=>'bg',-default=>'white');
print $query->hidden(-name=>'htmap',-default=>'checked');
print $query->hidden(-name=>'legend',-default=>'checked');

print $query->hidden(-name=>'xcol',-default=>'6');
print $query->hidden(-name=>'xleg1',-default=>'positional bias');
print $query->hidden(-name=>'xleg2',-default=>'(chi-square value)');
print $query->hidden(-name=>'xmin',-default=>'0');
print $query->hidden(-name=>'xmax',-default=>'auto');
print $query->hidden(-name=>'xsize',-default=>'500');
print $query->hidden(-name=>'xgstep1',-default=>'auto');
print $query->hidden(-name=>'xgstep2',-default=>'none');

print $query->hidden(-name=>'ycol',-default=>'2,3,4,5');
print $query->hidden(-name=>'yleg1',-default=>'z-score');
print $query->hidden(-name=>'yleg2',-default=>'(Markov models, different orders)');
print $query->hidden(-name=>'ymin',-default=>'auto');
print $query->hidden(-name=>'ymax',-default=>'auto');
print $query->hidden(-name=>'ysize',-default=>'400');
print $query->hidden(-name=>'ygstep1',-default=>'auto');
print $query->hidden(-name=>'ygstep2',-default=>'none');

print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;
#print "<TD><B><A HREF='help.XYgraph.html'>MANUAL</A></B></TD>\n";



print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>";


## credits
print "<hr>";
print <<EndCredits;
<H4>Credits</H4>
<UL>

<P> XYgraph has been written by Jacques van Helden (<A
HREF="mailto:jvanheld\@bigre.ulb.ac.be (Jacques van
Helden)">jvanheld\@bigre.ulb.ac.be</A>). This program can be used
freely by academic users via its web interface. For commercial users,
please read our <A HREF="disclaimer.html">disclaimer</A>.

<P> XYgraph uses a graphical library developed by Thomas Boutell (<A
HREF="mailto:boutell\@boutell.com (Thomas
Boutell)">boutell\@boutell.com</A>).

<p> Copyright statement for the graphical library GD: Portions
copyright 1994, 1995, 1996, 1997, 1998, by Cold Spring Harbor
Laboratory. Funded under Grant P41-RR02188 by the National Institutes
of Health.

<p> Portions copyright 1996, 1997, 1998, by Boutell.Com, Inc.

<p>GIF decompression code copyright 1990, 1991, 1993, by David Koblas
(koblas\@netcom.com).

<p> Non-LZW-based GIF compression code copyright 1998, by Hutchison
Avenue Software Corporation (http://www.hasc.com/, info\@hasc.com).

<p> Permission has been granted to copy and distribute gd in any
context, including a commercial application, provided that this notice
is present in user-accessible supporting documentation.

<p> This does not affect your ownership of the derived work itself,
and the intent is to assure proper credit for the authors of gd, not
to interfere with your productive \u\s\e of gd. If you have questions,
ask.  "Derived works" includes all programs that utilize the library.
Credit must be given in user-accessible documentation.

<p> Permission to \u\s\e, copy, modify, and distribute this software
and its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all copies
and that both that copyright notice and this permission notice appear
in supporting documentation. This software is provided "as is" without
express or implied warranty.  </UL>

EndCredits

## close the form
print $query->end_html;
exit();

