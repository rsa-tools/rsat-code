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
print $query->textarea(-name=>'data', -id=>'data',
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
print $query->textfield(-name=>'title',-id=>'title',
			-default=>$default{title},
			-size=>50,
			-maxlength=>80);

## second title
print "<br>Second Title", "&nbsp;"x2;
print $query->textfield(-name=>'title2',-id=>'title2',
			-default=>$default{title2},
			-size=>50,
			-maxlength=>80);

## point size
print "<br>Point size", "&nbsp;"x2;
print $query->textfield(-name=>'pointsize',-id=>'pointsize',
			-default=>$default{pointsize},
			-size=>50,
			-maxlength=>80);

## background color
print "<br>Background color";
print $query->popup_menu(-name=>'bg',-id=>'bg',
			 -Values=>["white","black","blue","gray"],
			 -default=>$default{bg});

## Image format
print "Format";
print $query->popup_menu(-name=>'format',-id=>'format',
			 -Values=>["jpg","png","eps","pdf"],
			 -default=>$default{format});

print "<br>";
## draw lines to join points
print $query->checkbox(-name=>'lines',-id=>'lines',
		       -checked=>$default{lines},
		       -label=>' Lines between points ');

## legend
print $query->checkbox(-name=>'legend',-id=>'legend',
		       -checked=>$default{legend},
		       -label=>' First row contains legend');

print "<br>";
## symbols
print $query->checkbox(-name=>'symbols',-id=>'symbols',
		       -checked=>$default{symbols},
		       -label=>' Use symbols');


## dynamic map
print $query->checkbox(-name=>'htmap',-id=>'htmap',
		       -checked=>$default{htmap},
		       -label=>' Dynamic map (information in the status bar))');



################################################################
## X axis options
print "<hr>\n";
print "<br><B>X axis</B><br>\n";

## xcol
print "<br>Data column for X axis", "&nbsp;"x2;
print $query->textfield(-name=>'xcol',-id=>'xcol',
			-default=>$default{xcol},
			-size=>5,
			);

## xsize
print "&nbsp;"x2, "Size (pixels)", "&nbsp;"x2;
print $query->textfield(-name=>'xsize',-id=>'xsize',
			-default=>$default{xsize},
			-size=>5,
			);

## Logarithmic scale for the X axis
print "&nbsp;"x2, "Log base", "&nbsp;"x2;
print $query->textfield(-name=>'xlog',-id=>'xlog',
			-default=>$default{xlog},
			-size=>4,
			);

## X labels
print "<br>First label", "&nbsp;"x2;
print $query->textfield(-name=>'xleg1',-id=>'xleg1',
			-default=>$default{xleg1},
			-size=>40,
			);
print "<br>Second label", "&nbsp;"x2;
print $query->textfield(-name=>'xleg2',-id=>'xleg2',
			-default=>$default{xleg2},
			-size=>40,
			);

## xmin
print "<br>Min value", "&nbsp;"x2;
print $query->textfield(-name=>'xmin',-id=>'xmin',
			-default=>$default{xmin},
			-size=>5,
			);
## xmax
print "&nbsp;"x2, "Max value", "&nbsp;"x2;
print $query->textfield(-name=>'xmax',-id=>'xmax',
			-default=>$default{xmax},
			-size=>5,
			);
## X grid
print "<br>First grid step", "&nbsp;"x2;
print $query->textfield(-name=>'xgstep1',-id=>'xgstep1',
			-default=>$default{xgstep1},
			-size=>5,
			);
print "&nbsp;"x2,"Second grid step", "&nbsp;"x2;
print $query->textfield(-name=>'xgstep2',-id=>'xgstep2',
			-default=>$default{xgstep2},
			-size=>5,
			);


################################################################
## Y axis options
print "<hr>\n";
print "<br><B>Y axis</B><br>\n";

## ycol
print "<br>Data columns for Y axis (e.g. 4,5,7)", "&nbsp;"x2;
print $query->textfield(-name=>'ycol',-id=>'ycol',
			-default=>$default{ycol},
			-size=>5,
			);

## ysize
print "&nbsp;"x2, "Size (pixels)", "&nbsp;"x2;
print $query->textfield(-name=>'ysize',-id=>'ysize',
			-default=>$default{ysize},
			-size=>5,
			);

## Logarithmic scale for the Y axis
print "&nbsp;"x2, "Log base", "&nbsp;"x2;
print $query->textfield(-name=>'ylog',-id=>'ylog',
			-default=>$default{ylog},
			-size=>4,
			);

## Y labels
print "<br>First label", "&nbsp;"x2;
print $query->textfield(-name=>'yleg1',-id=>'yleg1',
			-default=>$default{yleg1},
			-size=>40,
			);
print "<br>Second label", "&nbsp;"x2;
print $query->textfield(-name=>'yleg2',-id=>'yleg2',
			-default=>$default{yleg2},
			-size=>40,
			);

## ymin
print "<br>Min value", "&nbsp;"x2;
print $query->textfield(-name=>'ymin',-id=>'ymin',
			-default=>$default{ymin},
			-size=>5,
			);
## ymax
print "&nbsp;"x2, "Max value", "&nbsp;"x2;
print $query->textfield(-name=>'ymax',-id=>'ymax',
			-default=>$default{ymax},
			-size=>5,
			);
## Y grid
print "<br>First grid step", "&nbsp;"x2;
print $query->textfield(-name=>'ygstep1',-id=>'ygstep1',
			-default=>$default{ygstep1},
			-size=>5,
			);
print  "&nbsp;"x2,"Second grid step", "&nbsp;"x2;
print $query->textfield(-name=>'ygstep2',-id=>'ygstep2',
			-default=>$default{ygstep2},
			-size=>5,
			);

print "<hr>\n";

print "<UL><UL><TABLE class = 'formbutton'><TR>";
print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
print"<TD>", $query->reset(-id=>"reset"), "</TD>\n";
print $query->end_form;


### data for the demo 
$demo_data = "";
open($fh, "demo_files/XYgraph_demo_data.ft");
while($row = <$fh>){
    chomp $row;
    $demo_data .= $row;
    $demo_data .= "\\n";
}

print '<script>
function setDemo(demo_data){
    $("#reset").trigger("click");
    $("#data").val(demo_data);
    $("#title").val("XYgraph dmeo");
    $("#title2").val("score comparisons");
    pointsize.value = "5";
    $("#bg").val("white");
    $("#htmap").prop("checked", true);
    $("#legend").prop("checked", false);
    $("#xcol").val("6");
    $("#xleg1").val("positional bias");
    $("#xleg2").val("(chi-square value)");
    $("#xmin").val("0");
    $("#xmax").val("auto");
    $("#xsize").val("500");
    $("#xgstep1").val("auto");
    $("#xgstep2").val("none");
    
    $("#ycol").val("2,3,4,5");
    $("#yleg1").val("z-score");
    $("#yleg2").val("(Markov models, different orders)");
    $("#ymin").val("auto");
    $("#ymax").val("auto");
    $("#ysize").val("400");
    $("#ygstep1").val("auto");
    $("#ygstep2").val("none");
}
</script>';

print "<TD><B>";
print '<button type="button" onclick="setDemo('. "'$demo_data'" .')">DEMO</button>';
print "</B></TD>\n";
#print "<TD><B><A HREF='help.XYgraph.html'>MANUAL</A></B></TD>\n";



print "<TD><B><A HREF='mailto:Jacques.van-Helden\@univ-amu.fr'>MAIL</A></B></TD>\n";
print "</TR></TABLE></UL></UL>";


## credits
print "<hr>";
print <<EndCredits;
<H4>Credits</H4>
<UL>

<P> XYgraph has been written by Jacques van Helden (<A
HREF="mailto:Jacques.van-Helden\@univ-amu.fr (Jacques van
Helden)">Jacques.van-Helden\@univ-amu.fr</A>). This program can be used
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

