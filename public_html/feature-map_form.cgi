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
$feature_map_command = "$SCRIPTS/feature-map";
$tmp_file_name = sprintf "feature-map.%s", &AlphaDate();

# $features_from_swissprot_cmd = "$SCRIPTS/features-from-swissprot";
# $features_from_msf_cmd = "$SCRIPTS/features-from-msf";
# $features_from_gibbs_cmd = "$SCRIPTS/features-from-gibbs";
# $features_from_fugue_cmd = "$SCRIPTS/features-from-fugue";
# $features_from_dssp_cmd = "$SCRIPTS/features-from-dssp";
# $features_from_matins_cmd = "$SCRIPTS/features-from-matins";
# $features_from_sigscan_cmd = "$SCRIPTS/features-from-sigscan";
# $features_from_dnapat_cmd = "$SCRIPTS/features-from-dnapat";
# $features_from_tffact_cmd = "$SCRIPTS/features-from-tffact";
# $features_from_tfgene_cmd = "$SCRIPTS/features-from-tfgene";
# $features_from_patser_cmd = "$SCRIPTS/features-from-patser";


### Read the CGI query
$query = new CGI;

### read values to fill the form ###
$default{title} = "";
$default{data} = "";
$default{format} = 'feature map';
$default{img_format} = $ENV{rsat_img_format} || "jpg";
$default{from} = 'auto';
$default{to} = 'auto';
$default{handle} = 'none';
$default{origin} = '0';
$default{map_len} = 500;
$default{spacing} = 2;
$default{thickness} = 25;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
    if ($query->param($key)) {
	$default{$key} = $query->param($key);
    }
} 

#### specific treatment for internal feature file (piped from dna-patttern)
if (-e $query->param('feature_file')) {
    $file = $query->param('feature_file');
    $default{data} = `cat $file`;
} else {
    $default{data} = $query->param('data');
}
$default{data} =~ s/\"//g; #### remove quotes for security reasons (avoid imbedded command)


### remove quotes from the title and data
$default{data} =~ s/\"//g;
$default{title} =~ s/\"//g;


### print the form ###
&RSA_header("feature map", "form");

print "<CENTER>";
print "Generates a graphical map of features localized on one or several sequences.<P>\n";
print "</CENTER>";

print "<FONT FACE='Helvetica'>";

print $query->start_multipart_form(-action=>"feature-map.cgi");

print "<B>Feature list</B>&nbsp;&nbsp;&nbsp;&nbsp;(";
print "<A HREF='help.feature-map.html#formats'>";
print "Format</A>&nbsp;";

print $query->popup_menu(-name=>'format',
			 -Values=>['feature map', 
				   'dna-pattern',
				   'Patser',
				   'Matinspector',
				   'Signal scan',
				   'Swissprot',
				   'Transfac factor',
				   'Transfac gene',
				   'GCG msf',
				   'DSSP',
				   'Gibbs sampler',
				   'Fugue'],
			 -default=>$default{format});
print ")<BR>\n";


print $query->textarea(-name=>'data',
		       -default=>$default{data},
		       -rows=>6,
		       -columns=>60);
print "<BR>\n";

print "<B><A HREF='help.feature-map.html#file'>";
print "File ";
print "</A></B>&nbsp;";
print $query->filefield(-name=>'uploaded_file',
			-size=>45);


print "<B><A HREF='help.feature-map.html#title'>";

print "<P>";
print "Title";
print "</A></B>&nbsp;";
print $query->textfield(-name=>'title',
			-default=>$default{title},
			-size=>50,
			-maxlength=>80);

print "<BR>\n";
print "<div id='grey_background'>\n";
print $query->checkbox(-name=>'legend',
		       -checked=>'checked',
		       -label=>'');
print "<A HREF='help.feature-map.html#legend'><b>Legend</b></a>";

print "&nbsp"x5;
print $query->checkbox(-name=>'scalebar',
		       -checked=>'checked',
		       -label=>'');
print "<A HREF='help.feature-map.html#scalebar'><b>Scalebar</b></a>";

print "&nbsp;"x2, "step ";
print $query->textfield(-name=>'scalestep',
			-default=>'auto',
			-size=>5);
print "\n";
print "&nbsp;&nbsp;&nbsp;";

print "<br>\n";
print $query->checkbox(-name=>'seq_names',
		       -checked=>'checked',
		       -label=>'');
print "<A HREF='help.feature-map.html#seq_names'><b>Sequence names</b></a>";

print "&nbsp;"x5, "<B><A HREF='help.feature-map.html#orientation'>Orientation</B></A>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'orientation',
			 -Values=>['horizontal','vertical'],
			 -default=>'horizontal');

print "<BR>\n";
print "<B><A HREF='help.feature-map.html#limits'>Display limits</A></B>&nbsp;";
print "&nbsp;From&nbsp;";
print $query->textfield(-name=>'from',
			-default=>$default{from},
			-size=>5);

print "&nbsp;To&nbsp;";
print $query->textfield(-name=>'to',
			-default=>$default{to},
			-size=>5);
print "&nbsp;origin&nbsp;";
print $query->textfield(-name=>'origin',
			-default=>$default{origin},
			-size=>5);




print "<BR>\n";
print "<A HREF='help.feature-map.html#dimensions'><B>Map dimensions</B></A>\n";
print "&nbsp;Length&nbsp;";
print $query->textfield(-name=>'mlen',
			-default=>$default{map_len},
			-size=>5);

print "&nbsp;thickness&nbsp;";
print $query->textfield(-name=>'mapthick',
			-default=>$default{thickness},
			-size=>5);

print "&nbsp;spacing&nbsp;";
print $query->textfield(-name=>'mspacing',
			-default=>$default{spacing},
			-size=>5);

print "</div><BR>\n";
print "<div id='grey_background'>\n";
print "<B><A HREF='help.feature-map.html#palette'>Color palette</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'palette',
			 -Values=>['color','monochrome'],
			 -default=>'color');

print "<B><A HREF='help.feature-map.html#color_file'>&nbsp;&nbsp;";
print "Color File ";
print "</A></B>&nbsp;";
print $query->filefield(-name=>'color_file',
			-size=>10);

print "<BR>\n";
print "<B><A HREF='help.feature-map.html#bgcolor'>Background color (R,G,B)</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->textfield(-name=>'bgcolor',
			-default=>'220,220,220',
			-size=>11);


print "</div><BR>\n";
print "<div id='grey_background'>\n";
print "<B><A HREF='help.feature-map.html#handle'>Feature handle</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'handle',
			 -Values=>['color dot','symbol','none'],
			 -default=>$default{handle});

print "<BR>\n";
print "<B>Feature thickness</B> \n";
print "max ";
print $query->textfield(-name=>'maxfthick',
			-default=>'auto',
			-size=>5);
print "&nbsp;&nbsp;&nbsp;\n";

print "min ";
print $query->textfield(-name=>'minfthick',
			-default=>'auto',
			-size=>5);
print "&nbsp;&nbsp;&nbsp;\n";


print $query->checkbox(-name=>'scorethick',
		       -checked=>'checked',
		       -label=>'');
print "<A HREF='help.feature-map.html#scorethick'>Proportional to score</A>";
print "</B>";



print "<BR>\n";
print "<B>";

print $query->checkbox(-name=>'htmap',
		       -checked=>'checked',
		       -label=>'');
print "<A HREF='help.feature-map.html#htmap'>Dynamic map</A>";
print "&nbsp;&nbsp;&nbsp;";

print "<BR>\n";

print "<B><A HREF='help.feature-map.html#dynamic'>";
print "Label keys";
print "</A></B>\n";
print $query->checkbox(-name=>'label_strand',
		       -label=>' strand ');
print "&nbsp;&nbsp;&nbsp;";

print $query->checkbox(-name=>'label_pos',
		       -label=>' position ');
print "&nbsp;&nbsp;&nbsp;";

print $query->checkbox(-name=>'label_id',
		       -label=>' identifier ');
print "&nbsp;&nbsp;&nbsp;";

print $query->checkbox(-name=>'label_descr',
		       -label=>' description ');
print "&nbsp;&nbsp;&nbsp;";

print $query->checkbox(-name=>'label_score',
		       -label=>' score ');
print "&nbsp;&nbsp;&nbsp;";


## Image format
print "</div><BR>\n";
print "<B><A HREF='help.feature-map.html#img_format'>Image format</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'img_format',
			 -Values=>['jpg','png','ps'],
			 -default=>$default{img_format});


print "<P>";


print "<UL><UL><TABLE class = 'formbutton'><TR>";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print"<TD>", $query->reset, "</TD>\n";
print $query->end_form;


### data for the demo 
print $query->start_multipart_form(-action=>"feature-map_form.cgi");
$demo_data = "PHO5	dnapat	acgtgc|gcacgt	R	438	443	ACGTGC	4.37
PHO8	dnapat	acgtgc|gcacgt	D	268	273	ACGTGC	4.37
PHO11	dnapat	acgtgc|gcacgt	R	586	591	ACGTGC	4.37
PHO81	dnapat	acgtgc|gcacgt	D	458	463	ACGTGC	4.37
PHO81	dnapat	acgtgc|gcacgt	D	794	799	ACGTGC	4.37
PHO81	dnapat	acgtgc|gcacgt	R	456	461	ACGTGC	4.37
PHO84	dnapat	acgtgc|gcacgt	D	242	247	ACGTGC	4.37
PHO84	dnapat	acgtgc|gcacgt	D	540	545	ACGTGC	4.37
PHO84	dnapat	acgtgc|gcacgt	R	213	218	ACGTGC	4.37
PHO84	dnapat	acgtgc|gcacgt	R	386	391	ACGTGC	4.37
PHO5	dnapat	acgtgg|ccacgt	D	549	554	ACGTGG	2.84
PHO8	dnapat	acgtgg|ccacgt	R	266	271	ACGTGG	2.84
PHO11	dnapat	acgtgg|ccacgt	D	519	524	ACGTGG	2.84
PHO11	dnapat	acgtgg|ccacgt	R	384	389	ACGTGG	2.84
PHO81	dnapat	acgtgg|ccacgt	D	86	91	ACGTGG	2.84
PHO84	dnapat	acgtgg|ccacgt	D	366	371	ACGTGG	2.84
PHO84	dnapat	acgtgg|ccacgt	D	388	393	ACGTGG	2.84
PHO84	dnapat	acgtgg|ccacgt	R	364	369	ACGTGG	2.84
PHO5	dnapat	cgcacg|cgtgcg	D	218	223	CGCACG	1.50
PHO8	dnapat	cgcacg|cgtgcg	R	73	78	CGCACG	1.50
PHO11	dnapat	cgcacg|cgtgcg	D	585	590	CGCACG	1.50
PHO81	dnapat	cgcacg|cgtgcg	R	459	464	CGCACG	1.50
PHO84	dnapat	cgcacg|cgtgcg	D	212	217	CGCACG	1.50
PHO84	dnapat	cgcacg|cgtgcg	R	541	546	CGCACG	1.50
PHO5	dnapat	ctgcac|gtgcag	D	29	34	CTGCAC	1.77
PHO5	dnapat	ctgcac|gtgcag	D	411	416	CTGCAC	1.77
PHO5	dnapat	ctgcac|gtgcag	D	461	466	CTGCAC	1.77
PHO8	dnapat	ctgcac|gtgcag	R	270	275	CTGCAC	1.77
PHO8	dnapat	ctgcac|gtgcag	R	390	395	CTGCAC	1.77
PHO81	dnapat	ctgcac|gtgcag	D	671	676	CTGCAC	1.77
PHO84	dnapat	ctgcac|gtgcag	R	244	249	CTGCAC	1.77
PHO84	dnapat	ctgcac|gtgcag	R	457	462	CTGCAC	1.77
PHO5	dnapat	tgccaa|ttggca	D	469	474	TGCCAA	2.57
PHO5	dnapat	tgccaa|ttggca	D	610	615	TGCCAA	2.57
PHO5	dnapat	tgccaa|ttggca	R	36	41	TGCCAA	2.57
PHO5	dnapat	tgccaa|ttggca	R	538	543	TGCCAA	2.57
PHO5	dnapat	tgccaa|ttggca	R	638	643	TGCCAA	2.57
PHO11	dnapat	tgccaa|ttggca	D	597	602	TGCCAA	2.57
PHO11	dnapat	tgccaa|ttggca	D	663	668	TGCCAA	2.57
PHO81	dnapat	tgccaa|ttggca	R	287	292	TGCCAA	2.57
PHO84	dnapat	tgccaa|ttggca	D	184	189	TGCCAA	2.57
PHO84	dnapat	tgccaa|ttggca	D	348	353	TGCCAA	2.57
PHO84	dnapat	tgccaa|ttggca	D	634	639	TGCCAA	2.57
PHO84	dnapat	tgccaa|ttggca	R	598	603	TGCCAA	2.57
PHO5	dnapat	cacgtg|cacgtg	D	548	553	CACGTG	1.79
PHO5	dnapat	cacgtg|cacgtg	R	548	553	CACGTG	1.79
PHO8	dnapat	cacgtg|cacgtg	D	267	272	CACGTG	1.79
PHO8	dnapat	cacgtg|cacgtg	R	267	272	CACGTG	1.79
PHO11	dnapat	cacgtg|cacgtg	D	518	523	CACGTG	1.79
PHO11	dnapat	cacgtg|cacgtg	R	518	523	CACGTG	1.79
PHO81	dnapat	cacgtg|cacgtg	D	457	462	CACGTG	1.79
PHO81	dnapat	cacgtg|cacgtg	R	457	462	CACGTG	1.79
PHO84	dnapat	cacgtg|cacgtg	D	365	370	CACGTG	1.79
PHO84	dnapat	cacgtg|cacgtg	D	387	392	CACGTG	1.79
PHO84	dnapat	cacgtg|cacgtg	R	365	370	CACGTG	1.79
PHO84	dnapat	cacgtg|cacgtg	R	387	392	CACGTG	1.79
PHO5	dnapat	cccacg|cgtggg	D	197	202	CCCACG	0.45
PHO5	dnapat	cccacg|cgtggg	R	550	555	CCCACG	0.45
PHO11	dnapat	cccacg|cgtggg	R	300	305	CCCACG	0.45
PHO11	dnapat	cccacg|cgtggg	R	520	525	CCCACG	0.45
PHO84	dnapat	cccacg|cgtggg	R	389	394	CCCACG	0.45
PHO5	dnapat	aacgtg|cacgtt	D	77	82	AACGTG	0.22
PHO5	dnapat	aacgtg|cacgtt	R	439	444	AACGTG	0.22
PHO8	dnapat	aacgtg|cacgtt	D	421	426	AACGTG	0.22
PHO11	dnapat	aacgtg|cacgtt	R	385	390	AACGTG	0.22
PHO81	dnapat	aacgtg|cacgtt	D	793	798	AACGTG	0.22
PHO84	dnapat	aacgtg|cacgtt	D	539	544	AACGTG	0.22
PHO84	dnapat	aacgtg|cacgtt	R	214	219	AACGTG	0.22
PHO5	dnapat	aaacgt|acgttt	D	76	81	AAACGT	0.79
PHO5	dnapat	aaacgt|acgttt	D	695	700	AAACGT	0.79
PHO5	dnapat	aaacgt|acgttt	R	440	445	AAACGT	0.79
PHO8	dnapat	aaacgt|acgttt	D	420	425	AAACGT	0.79
PHO8	dnapat	aaacgt|acgttt	D	700	705	AAACGT	0.79
PHO8	dnapat	aaacgt|acgttt	R	702	707	AAACGT	0.79
PHO11	dnapat	aaacgt|acgttt	D	672	677	AAACGT	0.79
PHO11	dnapat	aaacgt|acgttt	R	386	391	AAACGT	0.79
PHO81	dnapat	aaacgt|acgttt	D	792	797	AAACGT	0.79
PHO81	dnapat	aaacgt|acgttt	R	54	59	AAACGT	0.79
PHO84	dnapat	aaacgt|acgttt	D	538	543	AAACGT	0.79
";
print "<TD><B>";
print $query->hidden(-name=>'data',-default=>$demo_data);
print $query->hidden(-name=>'title',-default=>'Motifs discovered in upstream sequences of 5 PHO genes');
print $query->hidden(-name=>'scalestep',-default=>'50');
print $query->hidden(-name=>'from',-default=>'0');
print $query->hidden(-name=>'to',-default=>'801');
print $query->hidden(-name=>'origin',-default=>'801');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;
print "<TD><B><A HREF='help.feature-map.html'>MANUAL</A></B></TD>\n";
##print "<TD><B><A HREF='demo.feature-map.html'>DEMO</A></B></TD>\n";



print "<TD><B><A HREF='mailto:jvanheld\@scmbb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>";


print "<HR>\n";


print "</FONT>";

print $query->end_html;

exit(0);






