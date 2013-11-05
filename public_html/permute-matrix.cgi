#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;

#### redirect error log to a file
#BEGIN {
#    $ERR_LOG = "/dev/null";
##    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
#    use CGI::Carp qw(carpout);
#    open (LOG, ">> $ERR_LOG")
#	|| die "Unable to redirect log\n";
#    carpout(*LOG);
#}

require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("permute-matrix result", 'results');

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = $SCRIPTS."/permute-matrix";
$prefix = "permute-matrix";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "permute-matrix.%s", &AlphaDate();
$ENV{rsat_echo} = 1;
@result_files = ();

#### read parameters ####
local $parameters = " -v 0";



################################################################
## Matrix input format
local $input_format = lc($query->param('matrix_format'));
($input_format) = split (/\s+/, $input_format);
$parameters .= " -in_format ".$input_format;

################################################################
#### Input matrix specification
$matrix_file = $tmp_file_path."_input.".$input_format;
if ($query->param('matrix')) {
  open MAT, "> ".$matrix_file;
  print MAT $query->param('matrix');
  close MAT;
  &DelayedRemoval($matrix_file);
  $parameters .= " -i ".$matrix_file;
} else {
  &RSAT::error::FatalError('You did not enter any data in the matrix box');
}
push @result_files, ("Input file",$matrix_file);


################################################################
## Permutations
if (&IsInteger($query->param('perm'))) {
    $parameters .= " -perm ".$query->param('perm');
}

################################################################
## Matrix output format
local $output_format = lc($query->param('output_format'));
$parameters .= " -out_format ".$output_format;

################################################################
## Output file
local $result_file = $tmp_file_path."_output.".$output_format;
push @result_files, ("Result file",$result_file);
$parameters .= " -o ".$result_file;

&ReportWebCommand($command." ".$parameters);

### execute the command ###
if ($query->param('output') eq "display") {
#    &PipingWarning();

 ## prepare figures
    ### prepare data for piping
 #   open RESULT, "$command $parameters |";
  &doit("$command $parameters"); ## DEBUG test
  open RESULT, "$result_file";  ## DEBUG

  print '<H4>Result</H4>';
  print '<PRE>';
  while (<RESULT>) {
      print $_;
  }
  print '</PRE>';
  close(RESULT);

  ################################################################
  ## Prepare tab-delimited matrices with only the counts f the first
  ## matrix, for piping the result to other programs
  local $tab_matrices = $tmp_file_path."_simple.tab";

  local $command = $SCRIPTS."/convert-matrix -v 0 -i  $matrix_file -from ".$input_format." -to tab  -return counts -o $tab_matrices";
  system $command;
  print "<pre><b>Tab conversion:</b> $command</pre>" if ($ENV{rsat_echo} >= 1);
#  local $matrix_content = `$command`;
  push @result_files, ("Tab matrices", $tab_matrices);

  &PrintURLTable(@result_files);
  &PipingForm();

    print "<HR SIZE = 3>";
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'),$result_file);
}
print $query->end_html;

exit(0);


### prepare data for piping
sub PipingForm {
  local $matrix_content = `cat $tab_matrices`;
  $matrix_content =~ s|//\n||gm;
  $matrix_content =~ s|;.*\n||gm;
#  print "<pre>".$command."</pre>";
#  print "<pre>".$matrix_content."</pre>";


  $title = $query->param('title');
  $title =~ s/\"/\'/g;
    print <<End_of_form;
<hr size="3">
<table class="Nextstep">
<tr>
<td colspan="3">
<h3>Next step</h3>
</td>
</tr>

<tr>
<!--
<td valign="bottom" align="center">
<form method="post" action="patser_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$tab_result_file">
<input type="hidden" name="matrix_format" value="tab">
<input type="submit" value="pattern matching (patser)">
</form>
</td>
-->

<td valign="bottom" align="center">
<form method="POST" action="matrix-scan_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="submit" value="pattern matching (matrix-scan)">
</form>
</td>

<td valign=bottom align=center>
<form method="post" action="convert-matrix_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$result_file">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="hidden" name="logo" value="on" checked="checked">
<input type="submit" value="convert-matrix">
</form>
</td>


</table>
End_of_form

#  print "<pre>", $matrix_content, "</pre>";
  $title = $query->param('title');

}
