#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "cgi-lib.pl";
$XYgraph_command = "$SCRIPTS/XYgraph";
$tmp_file_name = sprintf "XYgraph.%s", &AlphaDate();


## dfault parameters
$image_format = $IMG_FORMAT || "png";

### Read the CGI query
$query = new CGI;

#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ECHO >= 2);


#### read parameters ####
$parameters = "";

### general parameters ###
if ($query->param('title1')) {
    $parameters .= " -title1 \"".$query->param('title1')."\" ";
}
if ($query->param('title2')) {
    $parameters .= " -title2 \"".$query->param('title2')."\" ";
}
if (&IsNatural($query->param('pointsize'))) {
    $parameters .= " -pointsize ".$query->param('pointsize');
}
if ($query->param('lines') eq "yes") {
    $parameters .= " -lines ";
}
if ($query->param('symbols') eq "yes") {
    $parameters .= " -symbols ";
}
if ($query->param('legend') eq "yes") {
    $parameters .= " -legend ";
}
if (&IsNatural($query->param('label_col'))) {
    $parameters .= " -lc ".$query->param('label_col');
}

%accepted_bg = (
		"black"=>1,
		"white"=>1,
		"blue"=>1,
		"gray"=>1,
		);
if ($accepted_bg{$query->param('bg')}) {
    $parameters .= " -bg ".$query->param('bg');
}


    
### X axis parameters ###
if ($query->param('xcol')) {
    $parameters .= " -xcol ".$query->param('xcol');
}
if ($query->param('xlog') eq "yes") {
    $parameters .= " -xlog ";
}
if ($query->param('xleg1')) {
    $parameters .= " -xleg1 \"".$query->param('xleg1')."\"";
}
if ($query->param('xleg2')) {
    $parameters .= " -xleg2 \"".$query->param('xleg2')."\"";
}
if (&IsReal($query->param('xmin'))) {
    $parameters .= " -xmin ".$query->param('xmin');
}
if (&IsReal($query->param('xmax'))) {
    $parameters .= " -xmax ".$query->param('xmax');
}
if (&IsReal($query->param('xgstep1'))) {
    $parameters .= " -xgstep1 ".$query->param('xgstep1');
}
if (&IsReal($query->param('xgstep2'))) {
    $parameters .= " -xgstep2 ".$query->param('xgstep2');
}
if ($query->param('xsize')) {
    $parameters .= " -xsize ".$query->param('xsize');
}

### Y axis parameters ###
if ($query->param('ycol')) {
    $parameters .= " -ycol ".$query->param('ycol');
}
if ($query->param('ylog') eq "yes") {
    $parameters .= " -ylog ";
}
if ($query->param('yleg1')) {
    $parameters .= " -yleg1 \"".$query->param('yleg1')."\"";
}
if ($query->param('yleg2')) {
    $parameters .= " -yleg2 \"".$query->param('yleg2')."\"";
}
if (&IsReal($query->param('ymin'))) {
    $parameters .= " -ymin ".$query->param('ymin');
}
if (&IsReal($query->param('ymax'))) {
    $parameters .= " -ymax ".$query->param('ymax');
}
if (&IsReal($query->param('ygstep1'))) {
    $parameters .= " -ygstep1 ".$query->param('ygstep1');
}
if (&IsReal($query->param('ygstep2'))) {
    $parameters .= " -ygstep2 ".$query->param('ygstep2');
}
if ($query->param('ysize')) {
    $parameters .= " -ysize ".$query->param('ysize');
}

# ### Y axis parameters ###
# if ($query->param('ycol'}) {
#     $parameters .= " -ycol $query->param('ycol'} ";
# }
# if ($query->param('ylog'} eq "yes") {
#     $parameters .= " -ylog ";
# }
# if ($query->param('yleg1'}) {
# 	$parameters .= " -yleg1 \"$query->param('yleg1'}\" ";
#     }
# if ($query->param('yleg2'}) {
#     $parameters .= " -yleg2 \"$query->param('yleg2'}\" ";
# }
# if ($query->param('ymin'}) {
#     $parameters .= " -ymin $query->param('ymin'} ";
# }
# if ($query->param('ymax'}) {
#     $parameters .= " -ymax $query->param('ymax'} ";
# }
# if ($query->param('ygstep1'}) {
#     $parameters .= " -ygstep1 $query->param('ygstep1'} ";
# }
# if ($query->param('ygstep2'}) {
#     $parameters .= " -ygstep2 $query->param('ygstep2'} ";
# }
# if ($query->param('ysize'}) {
#     $parameters .= " -ysize $query->param('ysize'} ";
# }

### data file ####
if ($query->param('data_file')) {
    ### internal data file ###
    $parameters .= " -i ".$query->param('data_file');
} else {
    unless ($query->param('data') =~ /\S/) {
	&RSA_header("XYgraph");
	&FatalError("The data box should not be empty.");
    }
    $data_file = "$tmp_file_name.data";
    open DATA, ">$TMP/$data_file";
    print DATA $query->param('data');
    close DATA;
    $parameters .= " -i $TMP/$data_file ";
}

### graph file ###
$graph_file = "$tmp_file_name.${image_format}";
$parameters .= " -o $TMP/$graph_file ";


if ($query->param('htmap') eq "yes") {
    $htmap = 1;
    $htmap_file = "$tmp_file_name.html";
    $parameters .= " -htmap ";
    $parameters .= " -htmap > $TMP/$htmap_file ";
}

### execute the command ###
@data_report = `$XYgraph_command $parameters`;

if ($ECHO >= 2) {
    print &RSA_header("XYgraph");
    print "<PRE>$XYgraph_command $parameters</PRE>";
}

### print the result ###
if ($htmap) {
    print "Location: $WWW_TMP/$htmap_file", "\n\n";
} else {
    ### display the result ###
    print &RSA_header("XYgraph");
    print <<End_Header;
    <HEADER>
<TITLE>RSA-tools - XY graph result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
RSA-tools</A> - XY graph result</H3>
End_Header
    print "<CENTER><IMG SRC=\"$WWW_TMP/$graph_file\"></CENTER><P>\n";
    print "<H4 ALIGN=CENTER>Data report</H4>";
    print "<PRE>";
    print @data_report;
    print "</PRE>";
    print "<HR SIZE = 3>";
    print &HtmlBot;
}    

&DelayedRemoval("$TMP/$graph_file");
&DelayedRemoval("$TMP/$data_file");

exit(0);



