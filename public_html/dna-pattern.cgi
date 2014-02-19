#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
@result_files = ();

### Read the CGI query
$query = new CGI;

### print the header of the result page
&RSA_header("dna-pattern result ".$query->param("title"), "results");

&ListParameters if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();


$command = "$SCRIPTS/dna-pattern";
$prefix = "dna-pattern";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);

#### read parameters ####
$parameters = " -v 1";

### pattern file ####
unless ($query->param('patterns') =~ /\S/) {
  &cgiError("The pattern box should not be empty.");
}
$pattern_file = $tmp_file_path.".pat";
push @result_files, "Patterns", $pattern_file;
if (open PAT, ">$pattern_file") {
  print PAT $query->param('patterns');
  close PAT;
  &DelayedRemoval($pattern_file);
}
$parameters .= " -pl $pattern_file";


### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("input sequence",$sequence_file);
$parameters .= " -i $sequence_file -format $sequence_format";


### return matching positions
if ($query->param('match_positions')) {
  $parameters .= " -return sites";

  ### origin ###
  if ($query->param('origin') =~ /end/i) {
    $parameters .= " -origin -0";
  }

  #### flanking residues for the matching sequences
  if ($query->param('flanking') =~ /^\d+$/) {
    $parameters .= " -N ".$query->param('flanking');
  }

  #### match format
  if ($query->param('match_format') eq "fasta") {
    $parameters .= " -match_format fasta";
  }

}

### return sequence limits
if ($query->param('limits')) {
  $parameters .= " -return limits";
}

### return sequence limits
if ($query->param('notacgt')) {
  $parameters .= " -return notacgt";
}


### return match count ###
if ($query->param('counts')) {
  $parameters .= " -return counts";
  if (($query->param('threshold') =~ /^\d+$/) &&
      ($query->param('threshold') > 0)) {
    $parameters .= " -th ".$query->param('threshold');
  }
}

### return match scores
if ($query->param('scores')) {
  $parameters .= " -return scores";
}

### return match rank
if ($query->param('rank')) {
    $parameters .= " -return rank";
}

### Sort matches
if ($query->param('sort')) {
    $parameters .= " -sort";
}

### return match count table
if ($query->param('table')) {
    $parameters .= " -return table";
    ### add a rwo and a column with the totals
    if ($query->param('total')) {
	$parameters .= " -return total";
    }
}

### return match statistics
if ($query->param('stats')) {
  $parameters .= " -return stats";
}

### prevent overlapping matches
if (lc($query->param('noov')) eq "on") {
  $parameters .= " -noov";
}


### strands
if ($query->param('strands') =~ /both/i) {
  $parameters .= " -2str";
} elsif ($query->param('strands') =~ /direct/i) {
  $parameters .= " -1str";
} elsif  ($query->param('strands') =~ /reverse/i) {
  $parameters .= " -R";
}

### substitutions
if ($query->param('subst') =~ /^\d+$/) {
  $parameters .= " -subst ".$query->param('subst');
}

&ReportWebCommand($command." ".$parameters);

### execute the command ###
if ($query->param("output") =~ /display/i) {

  ### execute the command ###
  &PipingWarning() if ($query->param('match_positions'));

  $result_file = $tmp_file_path.".dnapat";
  push @result_files, ('dna-pattern result', $result_file);
  open RESULT, "$command $parameters |";

  ### Print the result on Web page
  print "<h2>Result</h2>";
  &PrintHtmlTable(RESULT, $result_file, 1);
  close RESULT;

  &PrintURLTable(@result_files);
  &PipingForm();
  print "<hr size='3'>";

} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $pattern_file);
}

print $query->end_html;


exit(0);

sub PipingForm {
    ### prepare data for piping
    $title = $query->param("title");
    $title =~ s/\"/\'/g;
print "
    <TABLE class='nextstep'>
<TR>
  <TD>
    <H3>Next step</H3>
  </TD>
  </tr>
";

    if ($query->param('match_positions')) {
      print "
  <tr>
  <td>
<FORM METHOD='POST' ACTION='feature-map_form.cgi'>
<INPUT type='hidden' NAME='title' VALUE='$title'>
<INPUT type='hidden' NAME='feature_file' VALUE='$result_file'>
<INPUT type='hidden' NAME='format' VALUE='dna-pattern'>
<INPUT type='hidden' NAME='fill_form' VALUE='on'>
<INPUT type='submit' VALUE='feature map'>
 </form>
  </td>
</tr>";
    }

    print "
  <tr>
  <td>
<FORM METHOD='POST' ACTION='classfreq_form.cgi'>
<INPUT type='hidden' NAME='transferred_file' VALUE='$result_file'>
<INPUT type='submit' VALUE='Frequency distribution'>
 </form>
  </td>
</tr>";


print "</table>";
}
