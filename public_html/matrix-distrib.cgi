#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
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

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("matrix-distrib result", 'results');

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = $SCRIPTS."/matrix-distrib -v 1 -top 1";
$prefix = "matrix-distrib";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
#$tmp_file_name = sprintf "matrix-distrib.%s", &AlphaDate();
#$result_file = "$TMP/$tmp_file_name.res";
@result_files = ();

#### read parameters ####
my $parameters;

################################################################
#### Matrix specification
$matrix_file = $tmp_file_path.".input";
open MAT, "> $matrix_file";
print MAT $query->param('matrix');
close MAT;
&DelayedRemoval($matrix_file);
push (@result_files, "input matrix",$matrix_file);

$parameters .= " -m $matrix_file";


################################################################
## pseudo-counts and weights are mutually exclusive
if (&IsReal($query->param('pseudo_counts'))) {
  $parameters .= " -pseudo ".$query->param('pseudo_counts');
} else {
  &FatalError("Pseudo-count should be a real number");
}

if ($query->param('pseudo_distribution') eq "equi_pseudo") {
  $parameters .= " -equi_pseudo ";
}
################################################################
## decimals
if (&IsInteger($query->param('decimals'))) {
  $parameters .= " -decimals ".$query->param('decimals');
} else {
  &FatalError("Decimals should be an integer number");
}


################################################################
## Matrix input format
my $input_format = lc($query->param('matrix_format'));
($input_format) = split (/\s+/, $input_format);
#$input_format =~ s/cluster\-buster/cb/i;
#$input_format =~ s/(\S+)/$1/; ## Only retain the first word
$parameters .= " -matrix_format ".$input_format;


################################################################
## Background model method
my $bg_method = $query->param('bg_method');

if ($bg_method eq "bgfile") {
  ## Select pre-computed background file in RSAT genome directory
  my $organism_name = $query->param("organism");
  my $noov = "ovlp";
  my $background_model = $query->param("background");
  #    my $oligo_length = 1;
  my $markov_order = $query->param('markov_order');
  my $oligo_length = $markov_order + 1;
  $bg_file = &ExpectedFreqFile($organism_name,
			       $oligo_length, $background_model,
			       noov=>$noov, str=>"-1str");
  $parameters .= " -bgfile ".$bg_file;

} elsif ($bg_method =~ /upload/i) {
  ## Upload user-specified background file
  my $bgfile = "${TMP}/${tmp_file_name}_bgfile.txt";
  my $upload_bgfile = $query->param('upload_bgfile');
  if ($upload_bgfile) {
    if ($upload_bgfile =~ /\.gz$/) {
      $bgfile .= ".gz";
    }
    my $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
    open BGFILE, ">$bgfile" ||
      &cgiError("Cannot store background file in temp dir.");
    while (<$upload_bgfile>) {
      print BGFILE;
    }
    close BGFILE;
    $parameters .= " -bgfile $bgfile";
    $parameters .= " -bg_format ".$query->param('bg_format');
  } else {
    &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
  }

} else {
  &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
}

################################################################
## bg_pseudo
if (&IsReal($query->param('bg_pseudo'))) {
  $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
}


&ReportWebCommand($command." ".$parameters);

### execute the command ###

my $error_found = 0; # catch an error if occurs, and then prevent from drawing the graphs

## Tab-delimited output file
$distrib_file = $tmp_file_path.".tab";
push (@result_files, "distribution table", $distrib_file);

if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
  &PipingWarning();
  if ($query->param('output') =~ /server/i) {
    &Info("The result will appear below ...");
  }

  ### prepare data for piping
  open RESULT, "$command $parameters |";

  ### open the sequence file on the server
  if (open MIRROR, ">$distrib_file") {
    $mirror = 1;
    &DelayedRemoval($distrib_file);
  }

  print "<PRE>";
  while (<RESULT>) {
    print "$_" unless ($query->param('output') =~ /server/i);
    print MIRROR $_ if ($mirror);
    if ($_ =~ /Error<\/h4><blockquote/ ) {
      $error_found = 1;
    }	
  }
  print "</PRE>";
  close RESULT;
  close MIRROR if ($mirror);

#   if ($query->param('output') =~ /server/i) {
#     $result_URL = "$ENV{rsat_www}/tmp/${tmp_file_name}.res";
#     print ("The result is available at the following URL: ", "\n<br>",
# 	   "<a href=${result_URL}>${result_URL}</a>",
# 	   "<p>\n");
#   }
  print "<hr/>";
  if (($error_found)&&($query->param('output') =~ /server/i)) {
    &RSAT::error::FatalError("Error has occured, check output file.");
  }

  ## prepare figures
  unless($error_found){
    my $plot_format = "png";
    my $XYgraph_command = "$SCRIPTS/XYgraph";
    my $graph_file1 = $tmp_file_path."_1.".${plot_format};
    #    my $graph_file1 = "$tmp_file_name"."_1.".${plot_format};
    my $figure = $graph_file1;
    my $command2 = "$XYgraph_command";
    $command2 .= " -i $distrib_file";
    $command2 .= " -title1 'Distribution of weights'";
    $command2 .= " -title2 'Score probability'";
    $command2 .= " -xcol 1 -ycol 2 -legend -lines -pointsize 1";
    $command2 .= " -xleg1 'weight'";
    $command2 .= " -yleg1 'frequency'";
    $command2 .= " -format ".${plot_format};
    $command2 .= " -o $figure";
    print "<pre>command2: $command2\n</pre>" if ($ENV{rsat_echo} >= 1);
    system($command2);
    $graph_URL1 = $ENV{rsat_www}."/tmp/".&RSAT::util::RelativePath($TMP, $figure);
    print "<center><a href = \"".$graph_URL1."\"><img src=\"".$graph_URL1."\" width='200'></a>";
#    print "<center><a href = \"$WWW_TMP/$graph_file1\"><IMG SRC=\"$WWW_TMP/$graph_file1\" width='200'></a>";
    &DelayedRemoval("$TMP/$graph_file1");
    push (@result_files, "Weight distrib plot", $graph_file1);

    my $graph_file2 = $tmp_file_path."_2.".${plot_format};
    #    my $graph_file2 = "$tmp_file_name"."_2.".${plot_format};
    $figure = $graph_file2;
    my $command3 = "$XYgraph_command";
    $command3 .= " -i $distrib_file";
    $command3 .= " -title1 'Distribution of weights  (log scale)'";
    $command3 .= " -title2 'Score probability and P-value'";
    $command3 .= " -xcol 1";
    $command3 .= " -ycol 2,4";
    $command3 .= " -legend";
    $command3 .= " -lines";
    $command3 .= " -pointsize 1 -ylog -ymax 1 -ymin 0";
    $command3 .= " -xleg1 'weight' -yleg1 'Frequency (log scale)'";
    $command3 .= " -format ${plot_format}";
    $command3 .= " -o $figure";
    print "<pre>command3: $command3\n</pre>" if ($ENV{rsat_echo} >= 1);
    system($command3);
    $graph_URL2 = $ENV{rsat_www}."/tmp/".&RSAT::util::RelativePath($TMP, $figure);
    print "<center><a href = \"".$graph_URL2."\"><img src=\"".$graph_URL2."\" width='200'></a>";
#    print "<a href = \"$WWW_TMP/$graph_file2\"><IMG SRC=\"$WWW_TMP/$graph_file2\" width='200'></a></CENTER><P>\n";
    &DelayedRemoval("$TMP/$graph_file2");
    push (@result_files, "P-value distrib plot", $graph_file2);

    ## Links to the result files
    &PrintURLTable(@result_files);

    ### prepare data for piping
    &PipingForm();
    print "<HR SIZE = 3>";

  }

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $distrib_file);
}
print $query->end_html;

exit(0);

sub PipingForm {

    ### prepare data for piping
    $title = "Distribution of weights";
    print <<End_of_form;
<CENTER>

<TABLE class='nextstep'>
<TR>
  <TD>
    <H3>Next step</H3>
  </TD>
  </tr>
  <tr>
  <TD>
<FORM METHOD="POST" ACTION="XYgraph_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="XYgraph_file" VALUE="$distrib_file">
<INPUT type="hidden" NAME="xleg1" VALUE="weight">
<INPUT type="hidden" NAME="yleg1" VALUE="frequency">
<INPUT type="hidden" NAME="lines" VALUE="on">
<INPUT type="submit" VALUE="XY graph">(Draw graph with user-defined parameters)
</FORM>

</TD>
</TR>
</TABLE>
</CENTER>
End_of_form
}



################################################################
#
# Pipe the result to another command
#
sub PipingFormForDistrib {
    
  print <<End_of_form;
<HR SIZE = 3>
<table class = 'nextstep'>
<tr><td colspan = 5><h3>next step</h3></td></tr>


<tr valign="top" align="center">

    <th align=center>
      	<font size=-1>
	Drawing<br>
	</font>
    </th>

    <td align=center>
	<FORM METHOD="POST" ACTION="XYgraph_form.cgi">
	<INPUT type="hidden" NAME="title1" VALUE="Distribution of weights">
	<INPUT type="hidden" NAME="XYgraph_file" VALUE="$result_URL">
	<INPUT type="hidden" NAME="xcol" VALUE="1">
	<INPUT type="hidden" NAME="ycol" VALUE="2">
	<INPUT type="submit" value="XY Graph">
	</FORM>
	Draws a XY graph from a table of numeric data 
    </td>

</tr>
</table>
End_of_form
}
