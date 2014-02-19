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
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "cgi-lib.pl";

### Read the CGI query
$query = new CGI;

$ENV{rsat_echo}=1;

## Open result page
&RSA_header("XYgraph result");

&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$XYgraph_command = "$SCRIPTS/XYgraph";
$prefix = "XYgraph";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);

#### read parameters ####
$parameters = " -v 1";

### general parameters ###
if ($query->param('title')) {
    $parameters .= " -title1 \"".$query->param('title')."\" ";
}
if ($query->param('title2')) {
    $parameters .= " -title2 \"".$query->param('title2')."\" ";
}
if (&IsNatural($query->param('pointsize'))) {
    $parameters .= " -pointsize ".$query->param('pointsize');
}
if ($query->param('lines')) {
    $parameters .= " -lines ";
}
if ($query->param('symbols')) {
    $parameters .= " -symbols ";
}
if ($query->param('legend')) {
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
if (&IsReal($query->param('xlog'))) {
    $parameters .= " -xlog ".$query->param('xlog');
}

### Y axis parameters ###
if ($query->param('ycol')) {
    $parameters .= " -ycol ".$query->param('ycol');
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
if (&IsReal($query->param('ylog'))) {
    $parameters .= " -ylog ".$query->param('ylog');
}

################################################################
## data file 
if ($query->param('data_file') =~ /\S/) {
    ### file on the server
    $data_file = $query->param('data_file');

} elsif ($query->param('uploaded_file')) {
    ### upload file from the client
    $data_file = $tmp_file_path."_data.tab";
    open DATA, ">$data_file";
    $upload_data_file = $query->param('uploaded_file');
    $type = $query->uploadInfo($upload_data_file)->{'Content-Type'};
    while (<$upload_data_file>) {
	print DATA;
    }
    close DATA;

### data from the textarea
} else {
    unless ($query->param('data') =~ /\S/) {
 	&RSA_header("XYgraph");
 	&FatalError("The data box should not be empty.");
    }
    $data_file = $tmp_file_path."_data.tab";
    open DATA, ">$data_file";
    print DATA $query->param('data');
    close DATA;
}
push @result_files, "Input data (tab)", $data_file;

$parameters .= " -i ".$data_file;

### graph file ###
$image_format = $query->param('format') || $ENV{rsat_img_format} || "png";
$graph_file = $tmp_file_path.".".$image_format;
push @result_files, "XY graph ($image_format)", $graph_file;
$parameters .= " -format ".$image_format;
$parameters .= " -o ".$graph_file;

if ($query->param('htmap')) {
    $htmap = 1;
    $htmap_file = $tmp_file_path.".html";
    push @result_files, "HTML map", $htmap_file;
#    $parameters .= " -htmap ";
    $parameters .= " -htmap > ".$htmap_file;
}

$XYgraph_command .= " ".$parameters;

### execute the command ###
@data_report = `$XYgraph_command`;

### print the result ###
### display the result ###
my $graph_URL = $ENV{rsat_www}."/tmp/"; $graph_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $graph_file);
my $data_URL = $ENV{rsat_www}."/tmp/"; $data_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $data_file);
if ($htmap) {
  my $short_graph_file = &ShortFileName($graph_file);
  my $short_data_file = &ShortFileName($data_file);
  my $htmap_content = `cat $htmap_file`;
  $htmap_content =~ s/<\/*html>//gi;
  $htmap_content =~ s/<\/*body>//gi;
  $htmap_content =~ s/<\/*head>//gi;
  $htmap_content =~ s/<title>.*<\/title>//gi;
  $htmap_content =~ s/${short_graph_file}/${graph_URL}/g;
  $htmap_content =~ s/${short_data_file}/${data_URL}/g;
#  $htmap_content =~ s/</&lt;/g;
#  $htmap_content =~ s/>/&gt;/g;
  print $htmap_content;
} else {
  print "<center><a href='".$graph_URL."'><img src='".$graph_URL."'></a></center><P>\n";
}

# ## Print out the data report
# ## OBSOLETE
# if (scalar(@data_report) >= 0) {
#   print "<h4 align='center'>Data report</h4>";
#   print "<pre>";
#   print join "\n",  @data_report;
#   print "</pre>";
# }

&PrintURLTable(@result_files);

&ReportWebCommand($XYgraph_command);


print "<hr size='3'>";
print &HtmlBot();

&DelayedRemoval($graph_file);
&DelayedRemoval($data_file);

exit(0);



