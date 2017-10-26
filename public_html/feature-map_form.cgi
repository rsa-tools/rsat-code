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
print "<A class='iframe' HREF='help.feature-map.html#formats'>";
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


print $query->textarea(-name=>'data',-id=>'data',
		       -default=>$default{data},
		       -rows=>6,
		       -columns=>60);
print "<BR>\n";

print "<B><A class='iframe' HREF='help.feature-map.html#file'>";
print "File ";
print "</A></B>&nbsp;";
print $query->filefield(-name=>'uploaded_file',
			-size=>45);


print "<B><A class='iframe' HREF='help.feature-map.html#title'>";

print "<P>";
print "Title";
print "</A></B>&nbsp;";
print $query->textfield(-name=>'title',-id=>'title',
			-default=>$default{title},
			-size=>50,
			-maxlength=>80);

print "<BR>\n";
print "<div id='grey_background'>\n";
print $query->checkbox(-name=>'legend',-id=>'legend',
		       -checked=>'checked',
		       -label=>'');
print "<A class='iframe' HREF='help.feature-map.html#legend'><b>Legend</b></a>";

print "&nbsp"x5;
print $query->checkbox(-name=>'scalebar',-id=>'scalebar',
		       -checked=>'checked',
		       -label=>'');
print "<A class='iframe' HREF='help.feature-map.html#scalebar'><b>Scalebar</b></a>";

print "&nbsp;"x2, "step ";
print $query->textfield(-name=>'scalestep',-id=>'scalestep',
			-default=>'auto',
			-size=>5);
print "\n";
print "&nbsp;&nbsp;&nbsp;";

print "<br>\n";
print $query->checkbox(-name=>'seq_names',
		       -checked=>'checked',
		       -label=>'');
print "<A class='iframe' HREF='help.feature-map.html#seq_names'><b>Sequence names</b></a>";

print "&nbsp;"x5, "<B><A class='iframe' HREF='help.feature-map.html#orientation'>Orientation</B></A>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'orientation',
			 -Values=>['horizontal','vertical'],
			 -default=>'horizontal');

print "<BR>\n";
print "<B><A class='iframe' HREF='help.feature-map.html#limits'>Display limits</A></B>&nbsp;";
print "&nbsp;From&nbsp;";
print $query->textfield(-name=>'from',-id=>'from',
			-default=>$default{from},
			-size=>5);

print "&nbsp;To&nbsp;";
print $query->textfield(-name=>'to',-id=>'to',
			-default=>$default{to},
			-size=>5);
print "&nbsp;origin&nbsp;";
print $query->textfield(-name=>'origin',-id=>'origin',
			-default=>$default{origin},
			-size=>5);




print "<BR>\n";
print "<A class='iframe' HREF='help.feature-map.html#dimensions'><B>Map dimensions</B></A>\n";
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
print "<B><A class='iframe' HREF='help.feature-map.html#palette'>Color palette</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'palette',
			 -Values=>['color','monochrome'],
			 -default=>'color');

print "<B><A class='iframe' HREF='help.feature-map.html#color_file'>&nbsp;&nbsp;";
print "Color File ";
print "</A></B>&nbsp;";
print $query->filefield(-name=>'color_file',
			-size=>10);

print "<BR>\n";
print "<B><A class='iframe' HREF='help.feature-map.html#bgcolor'>Background color (R,G,B)</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->textfield(-name=>'bgcolor',
			-default=>'220,220,220',
			-size=>11);


print "</div><BR>\n";
print "<div id='grey_background'>\n";
print "<B><A class='iframe' HREF='help.feature-map.html#handle'>Feature handle</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
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
print "<A class='iframe' HREF='help.feature-map.html#scorethick'>Proportional to score</A>";
print "</B>";



print "<BR>\n";
print "<B>";

print $query->checkbox(-name=>'htmap',
		       -checked=>'checked',
		       -label=>'');
print "<A class='iframe' HREF='help.feature-map.html#htmap'>Dynamic map</A>";
print "&nbsp;&nbsp;&nbsp;";

print "<BR>\n";

print "<B><A class='iframe' HREF='help.feature-map.html#dynamic'>";
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
print "<B><A class='iframe' HREF='help.feature-map.html#img_format'>Image format</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
print $query->popup_menu(-name=>'img_format',
			 -Values=>['jpg','png','ps'],
			 -default=>$default{img_format});


print "<P>";


print "<UL><UL><TABLE class = 'formbutton'><TR>";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print"<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;


### data for the demo 
$demo_data = "";
open ($fh, "demo_files/feature-map_demo_data.ft");
while($row = <$fh>){
    chomp $row;
    $demo_data .= $row;
    $demo_data .= "\\n";
}

print '<script>
function setDemo(demo_data){
    $("#reset").trigger("click");
    $("#data").val(demo_data);
    $("#title").val("Motifs discovered in upstream sequences of 19 MET genes");
    scalestep.value = "50";
    from.value = "-800";
    to.value = "0";
    origin.value = "0";
}
</script>';

print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_data'" .')">DEMO</button>';
print "</B></TD>\n";
print "<TD><B><A class='iframe' HREF='help.feature-map.html'>MANUAL</A></B></TD>\n";
##print "<TD><B><A HREF='demo.feature-map.html'>DEMO</A></B></TD>\n";



print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>";


print "<HR>\n";


print "</FONT>";

print $query->end_html;

exit(0);






