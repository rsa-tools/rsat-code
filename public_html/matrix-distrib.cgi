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
$command = "$SCRIPTS/matrix-distrib -v 1";
$tmp_file_name = sprintf "matrix-distrib.%s", &AlphaDate();
$result_file = "$TMP/$tmp_file_name.res";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("matrix-distrib result", 'results');

#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

#### read parameters ####
my $parameters;

################################################################
#### Matrix specification
$matrix_file = "$TMP/$tmp_file_name.input";
open MAT, "> $matrix_file";
print MAT $query->param('matrix');
close MAT;
&DelayedRemoval($matrix_file);

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


print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >= 1);

### execute the command ###

my $error_found = 0; # catch an error if occurs, and then prevent from drawing the graphs

if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
    	
    &PipingWarning();
    if ($query->param('output') =~ /server/i) {
	&Info("The result will appear below ...");
    }
    

    ### prepare data for piping
    open RESULT, "$command $parameters |";
    
    ### open the sequence file on the server
    $sequence_file = "$TMP/$tmp_file_name.res";
    if (open MIRROR, ">$sequence_file") {
	$mirror = 1;
	&DelayedRemoval($sequence_file);
    }

    print "<PRE>";
    while (<RESULT>) {
	print "$_" unless ($query->param('output') =~ /server/i);
	print MIRROR $_ if ($mirror);
	if ($_ =~ /Error<\/h4><blockquote/ ){
		$error_found = 1;
	}	
    }
    print "</PRE>";
    close RESULT;
    close MIRROR if ($mirror);
    
  

    if ($query->param('output') =~ /server/i) {
	$result_URL = "$ENV{rsat_www}/tmp/${tmp_file_name}.res";
	print ("The result is available at the following URL: ", "\n<br>",
	       "<a href=${result_URL}>${result_URL}</a>",
	       "<p>\n");
    }
    print "<hr/>";
    
   if (($error_found)&&($query->param('output') =~ /server/i)){
    	 &RSAT::error::FatalError("Error has occured, check output file.");
    }
    ## prepare figures
    unless($error_found){
    my $XYgraph_command = "$SCRIPTS/XYgraph";

    my $graph_file1 = "$tmp_file_name"."_1.png";
    my $figure = "$TMP/$graph_file1";
    my $command2 = "$XYgraph_command -i $sequence_file -o $figure -title1 'Distribution of weights' -title2 'Score probability' -xcol 1 -ycol 2 -legend -lines -pointsize 1 -xleg1 'weight' -yleg1 'frequency' -format png";
    `$command2`;
    print "<CENTER><a href = \"$WWW_TMP/$graph_file1\"><IMG SRC=\"$WWW_TMP/$graph_file1\" width='200'></a>";
    &DelayedRemoval("$TMP/$graph_file1");

    my $graph_file2 = "$tmp_file_name"."_2.png";
    $figure = "$TMP/$graph_file2";
    my $command3 = "$XYgraph_command -i $sequence_file -o $figure -title1 'Distribution of weights  (log scale)' -title2 'Score probability and P-value' -xcol 1 -ycol 2,4 -legend -lines -pointsize 1 -xleg1 'weight' -yleg1 'Frequency (log scale)' -format png -ylog -ymax 1 -ymin 0";
    `$command3`;
    print "<a href = \"$WWW_TMP/$graph_file2\"><IMG SRC=\"$WWW_TMP/$graph_file2\" width='200'></a></CENTER><P>\n";
    &DelayedRemoval("$TMP/$graph_file2");

    ### prepare data for piping
    &PipingForm();
    print "<HR SIZE = 3>";
    
    }

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
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
<INPUT type="hidden" NAME="XYgraph_file" VALUE="$sequence_file">
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
