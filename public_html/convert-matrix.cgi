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
##    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
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
&RSA_header("convert-matrix result", 'results');


## Check security issues
&CheckWebInput($query);

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = $SCRIPTS."/convert-matrix";
$prefix = "convert-matrix";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
$ENV{rsat_echo} = 1;
@result_files = ();

## Read parameters
local $parameters;


################################################################
## Matrix input format
local $input_format = lc($query->param('matrix_format'));
($input_format) = split (/\s+/, $input_format);
$parameters .= " -from ".$input_format;

################################################################
## Matrix output format
local $output_format = lc($query->param('output_format'));
$parameters .= " -to ".$output_format;

################################################################
#### Matrix specification
$matrix_file = $tmp_file_path."_input.".$input_format;
if ($query->param('matrix')) {
    open MAT, "> $matrix_file";
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
    $parameters .= " -i $matrix_file";
} else {
    &RSAT::error::FatalError('You did not enter any data in the matrix box');
}
push @result_files, ("Input file", $matrix_file);

## Result file
$result_file = $tmp_file_path."_output.".$output_format;
push @result_files, ("Result file",$result_file);

################################################################
## Compute reverse complement
if ($query->param('rc')) {
  $parameters .= " -rc";
}

################################################################
## Pseudo-counts
if (&IsReal($query->param('pseudo_counts'))) {
    $parameters .= " -pseudo ".$query->param('pseudo_counts');
} else {
    &FatalError("Pseudo-count should be a real number");
}
if ($query->param('pseudo_distribution') eq "equi_pseudo") {
    $parameters .= " -equi_pseudo ";
}

################################################################
## Multiply
if (&IsReal($query->param('multiply'))) {
    $parameters .= " -multiply ".$query->param('multiply');
} else {
  &RSAT::error::FatalError("Option 'Multiply counts' should be a Real number");
}

################################################################
## Insert columns on left and/or right flanks
if ($query->param('insert_col_left') > 0) {
  if (&IsNatural($query->param('insert_col_left'))) {
    $parameters .= " -insert_col_left ".$query->param('insert_col_left');
  }else {
    &RSAT::error::FatalError("Option 'Insert columns' should be a Natural number");
  }
}
if ($query->param('insert_col_right') > 0) {
  if (&IsNatural($query->param('insert_col_right'))) {
    $parameters .= " -insert_col_right ".$query->param('insert_col_right');
  }else {
    &RSAT::error::FatalError("Option 'Insert columns' should be a Natural number");
  }
}

################################################################
## decimals
if (&IsInteger($query->param('decimals'))) {
    $parameters .= " -decimals ".$query->param('decimals');
} else {
    &FatalError("Decimals should be an integer number");
}

################################################################
## permutations
if (&IsInteger($query->param('perm'))) {
    $parameters .= " -perm ".$query->param('perm');
}


################################################################
## Background model method
&SetBackgroundModel();

################################################################
## bg_pseudo
if (&IsReal($query->param('bg_pseudo'))) {
    $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
}


## Return fields
local @return_fields = ();
foreach my $stat qw (counts frequencies weights info consensus parameters profile header margins logo links) {
  if ($query->param($stat)) {
    push @return_fields, $stat;
    if ($stat eq "logo"){
      $parameters .= " -logo_format png,pdf ";
      $parameters .= " -logo_file ".$result_file."_logo";

      # seqlogo options
      if ($query->param("error_bar")){
	$parameters .= " -logo_opt '-e' ";
      }
      if ($query->param("small_correc")){
	$parameters .= " -logo_opt '-M' ";
      }
      if ($query->param("stretch")){
	$parameters .= " -logo_opt '-S' ";
      }
      $parameters .= " -logo_dir $ENV{RSAT}/public_html/tmp ";
    }
  }
}

if ($output_format eq 'tab') {
  ## verbosity
  if ($query->param("comments")) {
    $parameters .= " -v 1";
  }
  $parameters .= " -return ";
  $parameters .= join ",", @return_fields;
}else {
    $parameters .= " -to ".$output_format;
    $parameters .= " -return counts";
}
$parameters .= " -o ".$result_file;

&ReportWebCommand($command." ".$parameters);

## Update log file
&UpdateLogFile();

### execute the command ###
if ($query->param('output') eq "display") {
#    &PipingWarning();

 ## prepare figures
    ### prepare data for piping
#  open RESULT, "$command $parameters |";

  &doit("$command $parameters"); ## DEBUG test
  open RESULT, "$result_file";  ## DEBUG

  my $public_temp_dir = &RSAT::util::get_pub_temp();
  print '<H4>Result</H4>';
  print '<PRE>';
  while (<RESULT>) {
    next if ($_ =~ /logo file:(.*)\.pdf$/);
    if ($_ =~ /logo file:(.*)\.png$/){
      (local $logo = $1 )=~ s|${public_temp_dir}| ${WWW_TMP}|g;
      $logo =~ s/\.png//;
      print "<a href = '".$logo.".pdf'><IMG SRC='".$logo.".png' height='120'></a> ";
    } else {
      print $_;
    }
  }
  print '</PRE>';
  close(RESULT);

  ################################################################
  ## Prepare tab-delimited matrices with only the counts f the first
  ## matrix, for piping the result to other programs
  local $tab_matrices = $tmp_file_path."_simple.tab";
#  local $tab_matrices = $result_file.".tab";
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
<input type="hidden" name="matrix_file" value="$tab_matrices">
<input type="hidden" name="matrix_format" value="$output_format">
<input type="submit" value="pattern matching (matrix-scan)">
</form>
</td>

<td valign=bottom align=center>
<form method="post" action="convert-matrix_form.cgi">
<input type="hidden" name="title" value="$title">
<input type="hidden" name="matrix_file" value="$tab_matrices">
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
