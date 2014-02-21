#!/usr/bin/perl
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
$demo_data = "MET8		SEQ_START	DR	-463	-463	-	0.00
MET8		SEQ_END	DR	-1	-1	-	0.00
MET8		cacgtg|cacgtg	DR	-209	-204	atttCACGTGttat	9.02
MET8		acgtga|tcacgt	R	-210	-205	taacACGTGAaatt	4.86
MET8		aactgt|acagtt	R	-270	-265	acggAACTGTttga	0.64
MET8		agtcat|atgact	R	-151	-146	ccgtAGTCATcgga	0.60
MET8		tagtca|tgacta	DR	-363	-358	aaccTAGTCAattg	0.47
MET8		tagtca|tgacta	R	-150	-145	cccgTAGTCAtcgg	0.47
MET32		SEQ_START	DR	-547	-547	-	0.00
MET32		SEQ_END	DR	-1	-1	-	0.00
MET32		cacgtg|cacgtg	DR	-378	-373	catgCACGTGacat	9.02
MET32		acgtga|tcacgt	DR	-377	-372	atgcACGTGAcatt	4.86
MET32		ccacag|ctgtgg	DR	-326	-321	tacgCCACAGttta	3.69
MET32		gccaca|tgtggc	DR	-327	-322	gtacGCCACAgttt	3.41
MET32		actgtg|cacagt	R	-325	-320	ataaACTGTGgcgt	2.69
MET32		cgtgca|tgcacg	R	-380	-375	gtcaCGTGCAtgtt	1.43
MET32		aactgt|acagtt	R	-324	-319	tataAACTGTggcg	0.64
MET32		tagtca|tgacta	DR	-527	-522	atacTAGTCAaaga	0.47
MET32		cgtgac|gtcacg	DR	-376	-371	tgcaCGTGACattt	0.30
MET32		acgtgc|gcacgt	R	-379	-374	tgtcACGTGCatgt	0.11
MET18		SEQ_START	DR	-568	-568	-	0.00
MET18		SEQ_END	DR	-1	-1	-	0.00
MET18		acgtga|tcacgt	R	-490	-485	taggACGTGActgg	4.86
MET18		gccaca|tgtggc	R	-392	-387	gacaGCCACAcagg	3.41
MET18		cgtgca|tgcacg	DR	-101	-96	tttaCGTGCAatga	1.43
MET18		aactgt|acagtt	R	-415	-410	atagAACTGTagaa	0.64
MET18		agccac|gtggct	R	-391	-386	cgacAGCCACacag	0.34
MET18		cgtgac|gtcacg	R	-491	-486	aggaCGTGACtggt	0.30
MET18		acgtgc|gcacgt	DR	-102	-97	ttttACGTGCaatg	0.11
MET18		gactca|tgagtc	R	-278	-273	ctttGACTCAttga	0.05
MET30		SEQ_START	DR	-177	-177	-	0.00
MET30		SEQ_END	DR	-1	-1	-	0.00
MET30		cacgtg|cacgtg	DR	-177	-172	CACGTGatcg	9.02
MET30		acgtga|tcacgt	DR	-176	-171	cACGTGAtcgg	4.86
MET30		ccacag|ctgtgg	DR	-162	-157	gaagCCACAGtttg	3.69
MET30		gccaca|tgtggc	DR	-163	-158	ggaaGCCACAgttt	3.41
MET30		actgtg|cacagt	R	-161	-156	gcaaACTGTGgctt	2.69
MET30		aactgt|acagtt	R	-160	-155	cgcaAACTGTggct	0.64
MET30		agccac|gtggct	DR	-164	-159	gggaAGCCACagtt	0.34
MET30		cgcgca|tgcgcg	R	-154	-149	tctcCGCGCAaact	0.26
MET28		SEQ_START	DR	-489	-489	-	0.00
MET28		SEQ_END	DR	-1	-1	-	0.00
MET28		cacgtg|cacgtg	DR	-219	-214	agtgCACGTGactt	9.02
MET28		acgtga|tcacgt	DR	-218	-213	gtgcACGTGActta	4.86
MET28		ccacag|ctgtgg	DR	-284	-279	gatgCCACAGttgt	3.69
MET28		ccacag|ctgtgg	DR	-153	-148	aacaCCACAGtttt	3.69
MET28		gccaca|tgtggc	DR	-285	-280	agatGCCACAgttg	3.41
MET28		actgtg|cacagt	DR	-488	-483	gACTGTGataa	2.69
MET28		actgtg|cacagt	R	-283	-278	cacaACTGTGgcat	2.69
MET28		actgtg|cacagt	R	-152	-147	caaaACTGTGgtgt	2.69
MET28		cgtgca|tgcacg	DR	-389	-384	ctgcCGTGCAaaaa	1.43
MET28		cgtgca|tgcacg	R	-221	-216	gtcaCGTGCActca	1.43
MET28		aactgt|acagtt	R	-282	-277	ccacAACTGTggca	0.64
MET28		aactgt|acagtt	R	-151	-146	ccaaAACTGTggtg	0.64
MET28		agtcat|atgact	R	-82	-77	tcttAGTCATtggc	0.60
MET28		tagtca|tgacta	R	-81	-76	ttctTAGTCAttgg	0.47
MET28		cgtgac|gtcacg	DR	-217	-212	tgcaCGTGACttat	0.30
MET28		acgtgc|gcacgt	R	-220	-215	agtcACGTGCactc	0.11
MET28		gactca|tgagtc	R	-103	-98	ccttGACTCAcgcc	0.05
MET6		SEQ_START	DR	-687	-687	-	0.00
MET6		SEQ_END	DR	-1	-1	-	0.00
MET6		cacgtg|cacgtg	DR	-533	-528	acatCACGTGcaca	9.02
MET6		cacgtg|cacgtg	DR	-495	-490	atttCACGTGactt	9.02
MET6		acgtga|tcacgt	DR	-494	-489	tttcACGTGActta	4.86
MET6		acgtga|tcacgt	R	-534	-529	gtgcACGTGAtgtg	4.86
MET6		acgtga|tcacgt	R	-496	-491	agtcACGTGAaata	4.86
MET6		ccacag|ctgtgg	DR	-650	-645	tttcCCACAGccat	3.69
MET6		ccacag|ctgtgg	DR	-510	-505	aaagCCACAGgaaa	3.69
MET6		ccacag|ctgtgg	R	-304	-299	actaCCACAGtttt	3.69
MET6		gccaca|tgtggc	DR	-540	-535	tcacGCCACAtcac	3.41
MET6		gccaca|tgtggc	DR	-511	-506	aaaaGCCACAggaa	3.41
MET6		actgtg|cacagt	DR	-305	-300	caaaACTGTGgtag	2.69
MET6		actgtg|cacagt	DR	-269	-264	gagaACTGTGaatg	2.69
MET6		cgtgca|tgcacg	DR	-531	-526	atcaCGTGCAcatt	1.43
MET6		cgtgca|tgcacg	DR	-381	-376	cgaaCGTGCAaaag	1.43
MET6		aactgt|acagtt	DR	-306	-301	gcaaAACTGTggta	0.64
MET6		aactgt|acagtt	DR	-270	-265	tgagAACTGTgaat	0.64
MET6		agtcat|atgact	DR	-353	-348	ctagAGTCATcgca	0.60
MET6		agtcat|atgact	DR	-344	-339	tcgcAGTCATggca	0.60
MET6		agtcat|atgact	DR	-297	-292	tggtAGTCATagct	0.60
MET6		tagtca|tgacta	DR	-298	-293	gtggTAGTCAtagc	0.47
MET6		agccac|gtggct	DR	-512	-507	taaaAGCCACagga	0.34
MET6		cgtgac|gtcacg	DR	-493	-488	ttcaCGTGACttac	0.30
MET6		acgtgc|gcacgt	DR	-532	-527	catcACGTGCacat	0.11
MET6		acgtgc|gcacgt	DR	-382	-377	gcgaACGTGCaaaa	0.11
MET10		SEQ_START	DR	-338	-338	-	0.00
MET10		SEQ_END	DR	-1	-1	-	0.00
MET10		cacgtg|cacgtg	DR	-249	-244	acacCACGTGagct	9.02
MET10		cacgtg|cacgtg	DR	-231	-226	gaagCACGTGacca	9.02
MET10		acgtga|tcacgt	DR	-248	-243	caccACGTGAgctt	4.86
MET10		acgtga|tcacgt	DR	-230	-225	aagcACGTGAccac	4.86
MET10		ccacag|ctgtgg	DR	-212	-207	caccCCACAGgtgt	3.69
MET10		gccaca|tgtggc	R	-205	-200	aaaaGCCACAcctg	3.41
MET10		agtcat|atgact	DR	-150	-145	ccatAGTCATcttc	0.60
MET10		tagtca|tgacta	DR	-151	-146	cccaTAGTCAtctt	0.47
MET10		agccac|gtggct	R	-204	-199	aaaaAGCCACacct	0.34
MET10		cgtgac|gtcacg	DR	-229	-224	agcaCGTGACcaca	0.30
MET10		acgtgc|gcacgt	R	-232	-227	ggtcACGTGCttct	0.11
MET10		gactca|tgagtc	DR	-179	-174	aaaaGACTCAttca	0.05
MET13		SEQ_START	DR	-380	-380	-	0.00
MET13		SEQ_END	DR	-1	-1	-	0.00
MET13		ccacag|ctgtgg	DR	-115	-110	aatgCCACAGcttg	3.69
MET13		ccacag|ctgtgg	DR	-30	-25	accaCCACAGttac	3.69
MET13		gccaca|tgtggc	DR	-116	-111	taatGCCACAgctt	3.41
MET13		actgtg|cacagt	DR	-277	-272	atagACTGTGaata	2.69
MET13		actgtg|cacagt	R	-29	-24	agtaACTGTGgtgg	2.69
MET13		aactgt|acagtt	R	-28	-23	tagtAACTGTggtg	0.64
MET13		agtcat|atgact	DR	-302	-297	gtttAGTCATtcta	0.60
MET13		tagtca|tgacta	DR	-303	-298	cgttTAGTCAttct	0.47
MET13		gactca|tgagtc	DR	-251	-246	ctttGACTCAaatt	0.05
MET3		SEQ_START	DR	-800	-800	-	0.00
MET3		SEQ_END	DR	-1	-1	-	0.00
MET3		cacgtg|cacgtg	DR	-377	-372	aggtCACGTGacca	9.02
MET3		cacgtg|cacgtg	DR	-360	-355	aagtCACGTGtaat	9.02
MET3		acgtga|tcacgt	DR	-376	-371	ggtcACGTGAccag	4.86
MET3		acgtga|tcacgt	R	-378	-373	ggtcACGTGAcctt	4.86
MET3		acgtga|tcacgt	R	-361	-356	ttacACGTGActtt	4.86
MET3		ccacag|ctgtgg	DR	-515	-510	gggtCCACAGatat	3.69
MET3		ccacag|ctgtgg	DR	-253	-248	aaagCCACAGtttt	3.69
MET3		gccaca|tgtggc	DR	-254	-249	caaaGCCACAgttt	3.41
MET3		actgtg|cacagt	R	-252	-247	taaaACTGTGgctt	2.69
MET3		cgtgca|tgcacg	DR	-318	-313	ctgtCGTGCAcact	1.43
MET3		aactgt|acagtt	R	-251	-246	gtaaAACTGTggct	0.64
MET3		aactgt|acagtt	R	-78	-73	tcacAACTGTtacg	0.64
MET3		agccac|gtggct	DR	-255	-250	acaaAGCCACagtt	0.34
MET3		cgtgac|gtcacg	DR	-375	-370	gtcaCGTGACcaga	0.30
MET3		cgtgac|gtcacg	R	-379	-374	gtcaCGTGACcttt	0.30
MET3		cgtgac|gtcacg	R	-362	-357	tacaCGTGACtttt	0.30
MET3		cgcgca|tgcgcg	DR	-412	-407	ccaaCGCGCAtcct	0.26
MET14		SEQ_START	DR	-800	-800	-	0.00
MET14		SEQ_END	DR	-1	-1	-	0.00
MET14		cacgtg|cacgtg	DR	-228	-223	atttCACGTGatca	9.02
MET14		acgtga|tcacgt	DR	-227	-222	tttcACGTGAtcaa	4.86
MET14		acgtga|tcacgt	R	-229	-224	gatcACGTGAaatt	4.86
MET14		gccaca|tgtggc	R	-192	-187	cattGCCACAtttt	3.41
MET14		cgtgca|tgcacg	R	-268	-263	aaagCGTGCAgcta	1.43
MET14		aactgt|acagtt	R	-336	-331	aagcAACTGTattc	0.64
MET14		agtcat|atgact	R	-159	-154	ttgaAGTCATattg	0.60
MET14		agccac|gtggct	DR	-207	-202	tacaAGCCACctca	0.34
MET14		gactca|tgagtc	DR	-167	-162	gaacGACTCAatat	0.05
MET1		SEQ_START	DR	-702	-702	-	0.00
MET1		SEQ_END	DR	-1	-1	-	0.00
MET1		cacgtg|cacgtg	DR	-268	-263	tatgCACGTGacat	9.02
MET1		acgtga|tcacgt	DR	-267	-262	atgcACGTGAcatt	4.86
MET1		actgtg|cacagt	DR	-224	-219	ataaACTGTGaacg	2.69
MET1		cgtgca|tgcacg	R	-270	-265	gtcaCGTGCAtaat	1.43
MET1		aactgt|acagtt	DR	-406	-401	tcgcAACTGTagaa	0.64
MET1		aactgt|acagtt	DR	-225	-220	aataAACTGTgaac	0.64
MET1		cgtgac|gtcacg	DR	-266	-261	tgcaCGTGACatta	0.30
MET1		acgtgc|gcacgt	R	-269	-264	tgtcACGTGCataa	0.11
MET1		gactca|tgagtc	DR	-214	-209	aacgGACTCAtaat	0.05
MET1		gactca|tgagtc	R	-428	-423	gttgGACTCAgaga	0.05
MET17		SEQ_START	DR	-800	-800	-	0.00
MET17		SEQ_END	DR	-1	-1	-	0.00
MET17		cacgtg|cacgtg	DR	-300	-295	atggCACGTGaagc	9.02
MET17		acgtga|tcacgt	DR	-299	-294	tggcACGTGAagct	4.86
MET17		ccacag|ctgtgg	R	-274	-269	accaCCACAGttcc	3.69
MET17		actgtg|cacagt	DR	-275	-270	gggaACTGTGgtgg	2.69
MET17		actgtg|cacagt	DR	-219	-214	gaaaACTGTGtaac	2.69
MET17		aactgt|acagtt	DR	-276	-271	ggggAACTGTggtg	0.64
MET17		aactgt|acagtt	DR	-220	-215	tgaaAACTGTgtaa	0.64
MET17		agtcat|atgact	R	-645	-640	aaatAGTCATctaa	0.60
MET17		agtcat|atgact	R	-258	-253	aattAGTCATttgc	0.60
MET17		tagtca|tgacta	DR	-244	-239	aagtTAGTCAaggc	0.47
MET17		tagtca|tgacta	R	-644	-639	gaaaTAGTCAtcta	0.47
MET17		tagtca|tgacta	R	-257	-252	taatTAGTCAtttg	0.47
MET17		acgtgc|gcacgt	R	-301	-296	cttcACGTGCcatt	0.11
MET2		SEQ_START	DR	-481	-481	-	0.00
MET2		SEQ_END	DR	-1	-1	-	0.00
MET2		cacgtg|cacgtg	DR	-353	-348	ttttCACGTGatgc	9.02
MET2		acgtga|tcacgt	DR	-352	-347	tttcACGTGAtgcg	4.86
MET2		acgtga|tcacgt	R	-354	-349	catcACGTGAaaat	4.86
MET2		gccaca|tgtggc	DR	-326	-321	aggcGCCACAcatt	3.41
MET2		agtcat|atgact	DR	-454	-449	cactAGTCATgaaa	0.60
MET2		agtcat|atgact	DR	-274	-269	actaAGTCATgtta	0.60
MET2		tagtca|tgacta	DR	-455	-450	gcacTAGTCAtgaa	0.47
MET2		tagtca|tgacta	DR	-76	-71	ttttTAGTCAcagg	0.47
MET2		agccac|gtggct	R	-378	-373	ggggAGCCACtaac	0.34
MET4		SEQ_START	DR	-800	-800	-	0.00
MET4		SEQ_END	DR	-1	-1	-	0.00
MET4		cacgtg|cacgtg	DR	-330	-325	taatCACGTGcgcg	9.02
MET4		acgtga|tcacgt	R	-331	-326	gcgcACGTGAttaa	4.86
MET4		gccaca|tgtggc	DR	-264	-259	tcggGCCACAcaag	3.41
MET4		actgtg|cacagt	DR	-10	-5	ttccACTGTGaacg	2.69
MET4		cgtgca|tgcacg	R	-695	-690	tgtgCGTGCAtgta	1.43
MET4		cgtgca|tgcacg	R	-681	-676	tgcgCGTGCAtgta	1.43
MET4		aactgt|acagtt	DR	-76	-71	cgtcAACTGTttag	0.64
MET4		agtcat|atgact	R	-566	-561	ggttAGTCATcgct	0.60
MET4		tagtca|tgacta	DR	-111	-106	atatTAGTCAactc	0.47
MET4		tagtca|tgacta	R	-565	-560	aggtTAGTCAtcgc	0.47
MET4		cgcgca|tgcgcg	DR	-677	-672	tgcaCGCGCAtgac	0.26
MET4		cgcgca|tgcgcg	R	-583	-578	ccgaCGCGCAggaa	0.26
MET4		cgcgca|tgcgcg	R	-326	-321	taccCGCGCAcgtg	0.26
MET4		acgtgc|gcacgt	DR	-329	-324	aatcACGTGCgcgg	0.11
MET22		SEQ_START	DR	-215	-215	-	0.00
MET22		SEQ_END	DR	-1	-1	-	0.00
MET22		cacgtg|cacgtg	DR	-135	-130	atatCACGTGttgc	9.02
MET22		acgtga|tcacgt	R	-136	-131	caacACGTGAtatt	4.86
MET22		gccaca|tgtggc	DR	-92	-87	ccatGCCACAtata	3.41
MET22		agtcat|atgact	R	-155	-150	tctgAGTCATtcgc	0.60
MET22		cgcgca|tgcgcg	R	-162	-157	cattCGCGCAtttc	0.26
MET22		gactca|tgagtc	DR	-153	-148	gaatGACTCAgacg	0.05
MET7		SEQ_START	DR	-250	-250	-	0.00
MET7		SEQ_END	DR	-1	-1	-	0.00
MET7		cgcgca|tgcgcg	DR	-190	-185	tcgtCGCGCAtatt	0.26
MET31		SEQ_START	DR	-161	-161	-	0.00
MET31		SEQ_END	DR	-1	-1	-	0.00
MET31		actgtg|cacagt	R	-25	-20	atcgACTGTGtata	2.69
MET12		SEQ_START	DR	-384	-384	-	0.00
MET12		SEQ_END	DR	-1	-1	-	0.00
MET12		ccacag|ctgtgg	R	-214	-209	atagCCACAGtcaa	3.69
MET12		gccaca|tgtggc	R	-213	-208	tataGCCACAgtca	3.41
MET12		actgtg|cacagt	DR	-215	-210	attgACTGTGgcta	2.69
MET12		agccac|gtggct	R	-212	-207	ctatAGCCACagtc	0.34
MET16		SEQ_START	DR	-443	-443	-	0.00
MET16		SEQ_END	DR	-1	-1	-	0.00
MET16		cacgtg|cacgtg	DR	-178	-173	atttCACGTGgcta	9.02
MET16		acgtga|tcacgt	R	-179	-174	agccACGTGAaatg	4.86
MET16		gccaca|tgtggc	DR	-157	-152	aaaaGCCACAacat	3.41
MET16		agtcat|atgact	R	-149	-144	gctgAGTCATgttg	0.60
MET16		agccac|gtggct	DR	-158	-153	gaaaAGCCACaaca	0.34
MET16		agccac|gtggct	R	-175	-170	tactAGCCACgtga	0.34
MET16		cgcgca|tgcgcg	R	-361	-356	tcagCGCGCAttaa	0.26
MET16		gactca|tgagtc	DR	-147	-142	acatGACTCAgcaa	0.05
";
print "<TD><B>";
print $query->hidden(-name=>'data',-default=>$demo_data);
print $query->hidden(-name=>'title',-default=>'Motifs discovered in upstream sequences of 19 MET genes');
print $query->hidden(-name=>'scalestep',-default=>'50');
print $query->hidden(-name=>'from',-default=>'-800');
print $query->hidden(-name=>'to',-default=>'0');
print $query->hidden(-name=>'origin',-default=>'0');
print $query->submit(-label=>"DEMO");
print "</B></TD>\n";
print $query->end_form;
print "<TD><B><A HREF='help.feature-map.html'>MANUAL</A></B></TD>\n";
##print "<TD><B><A HREF='demo.feature-map.html'>DEMO</A></B></TD>\n";



print "<TD><B><A HREF='mailto:jvanheld\@bigre.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>";


print "<HR>\n";


print "</FONT>";

print $query->end_html;

exit(0);






