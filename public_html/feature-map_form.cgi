#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";

### intialization
$feature_map_command = "$SCRIPTS/feature-map";
$tmp_file_name = sprintf "feature-map.%s", &AlphaDate;

$features_from_swissprot_cmd = "$SCRIPTS/features-from-swissprot";
$features_from_msf_cmd = "$SCRIPTS/features-from-msf";
$features_from_gibbs_cmd = "$SCRIPTS/features-from-gibbs";
$features_from_fugue_cmd = "$SCRIPTS/features-from-fugue";
$features_from_dssp_cmd = "$SCRIPTS/features-from-dssp";
$features_from_matins_cmd = "$SCRIPTS/features-from-matins";
$features_from_sigscan_cmd = "$SCRIPTS/features-from-sigscan";
$features_from_dnapat_cmd = "$SCRIPTS/features-from-dnapat";
$features_from_tffact_cmd = "$SCRIPTS/features-from-tffact";
$features_from_tfgene_cmd = "$SCRIPTS/features-from-tfgene";
$features_from_patser_cmd = "$SCRIPTS/features-from-patser";


### Read the CGI query
$query = new CGI;

### read values to fill the form ###
$default{title} = $query->param('title');
$default{title} =~ s/\"//g;


if (-e $query->param('feature_file')) {
    $file = $query->param('feature_file');
    $default{data} = `cat $file`;
} else {
    $default{data} = $query->param('data');
}
$default{data} =~ s/\"//g;


if ($query->param('format')) {
    $default{format} = $query->param('format');
} else {
    $default{format} = 'feature map';
}

if ($query->param('from')) {
    $default{from} = $query->param('from');
} else {
    $default{from} = 'auto';
}

if ($query->param('to')) {
    $default{to} = $query->param('to');
} else {
    $default{to} = 'auto';
}

if ($query->param('handle')) {
    $default{handle} = $query->param('handle');
} else {
    $default{handle} = 'symbol';
}

if ($query->param('origin')) {
    $default{origin} = $query->param('origin');
} else {
    $default{origin} = '0';
}


### print the form ###
&RSA_header("feature map");

print "<CENTER>";
print "Generates a physical map of genetic features for one or several sequences<P>\n";
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
print "<A HREF='help.feature-map.html#legend'>";
print $query->checkbox(-name=>'legend',
		       -checked=>'checked',
		       -label=>' Legend ');
print "</A>";
print "&nbsp;&nbsp;&nbsp;";

print "<A HREF='help.feature-map.html#scalebar'>";
print $query->checkbox(-name=>'scalebar',
		       -checked=>'checked',
		       -label=>' Scalebar ');
print "</A>";

print "(step ";
print $query->textfield(-name=>'scalebarstep',
			-default=>'auto',
			-size=>5);
print ")\n";
print "&nbsp;&nbsp;&nbsp;";

print "<B><A HREF='help.feature-map.html#orientation'>Orientation</B></A>&nbsp;&nbsp;&nbsp;&nbsp";
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
			-default=>'500',
			-size=>5);

print "&nbsp;thickness&nbsp;";
print $query->textfield(-name=>'mapthick',
			-default=>'auto',
			-size=>5);

print "&nbsp;spacing&nbsp;";
print $query->textfield(-name=>'mspacing',
			-default=>'6',
			-size=>5);

print "<BR>\n";
print "<B><A HREF='help.feature-map.html#handle'>Feature handle</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'handle',
			 -Values=>['color dot','symbol','none'],
			 -default=>$default{handle});

print "<B><A HREF='help.feature-map.html#palette'>Color palette</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'palette',
			 -Values=>['color','monochrome'],
			 -default=>'color');


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

print "<A HREF='help.feature-map.html#scorethick'>";
print $query->checkbox(-name=>'scorethick',
		       -checked=>'checked',
		       -label=>' Proportional to score ');
print "</A>";
print "</B>";



print "<BR>\n";
print "<B>";
print "<A HREF='help.feature-map.html#htmap'>";
print $query->checkbox(-name=>'htmap',
		       -checked=>'checked',
		       -label=>' Dynamic ');
print "</A>";
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


print "<P>";


print "<UL><UL><TABLE><TR>";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print"<TD>", $query->reset, "</TD>\n";
print "<TD><B><A HREF='help.feature-map.html'>MANUAL</A></B></TD>\n";
print "<TD><B><A HREF='demo.feature-map.html'>DEMO</A></B></TD>\n";
print "<TD><B><A HREF='mailto:jvanheld\@ucmb.ulb.ac.be'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>";

print "</FONT>";

print $query->end_form;
print "<HR>\n";
print $query->end_html;

exit(0);






