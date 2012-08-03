#!/usr/bin/perl
#### redirect error log to a file
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

#### TEMPORARY
#$ENV{rsat_echo}=2;
@result_files = ();

$position_analysis_command = "$SCRIPTS/position-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$prefix = "position-analysis";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "position-analysis.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("position-analysis result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

#### read parameters ####
$parameters = " -v ";
$parameters .= " -sort" if ($query->param('sort'));
$parameters .= " -nofilter" unless ($query->param('filter'));
$parameters .= " -nocheck" unless ($query->param('check'));

#### purge sequence option
$purge = $query->param('purge');

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("Input sequence", $sequence_file);
if ($purge) {
  $purged_seq_file = ${sequence_file}.".purged";
  push @result_files, ("Purged sequence", $purged_seq_file);

  $command = "$purge_sequence_command -i $sequence_file -format $sequence_format -o $purged_seq_file";
  $command .= "; $position_analysis_command -i ".$purged_seq_file;
#    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format -o ${sequence_file}.purged;  $position_analysis_command -i ${sequence_file}.purged  ";
} else {
    $command= $position_analysis_command." -i ".$sequence_file;
}


### fields to return
@return_fields = ();

### threshold on occurrences
$lth_occ = $query->param('lth_occ');
&FatalError("$lth_occ Invalid threshold on occurrences") unless &IsNatural($lth_occ);
$parameters .= " -lth_occ ".$lth_occ;

#### chi2
if ($query->param('return_chi')) {
  push @return_fields, "chi";

  ### threshold on chi-square
  $lth = $query->param('lth');
  &FatalError("$lth Invalid threshold on chi2") unless &IsNatural($lth);
  $parameters .= " -lth ".$lth;
}

### graphs
if ($query->param('return_graph')) {
  push @return_fields, "graph";
} 

### rank
if ($query->param('return_rank')) {
  push @return_fields, "rank";
}

### distrib
if ($query->param('return_distrib')) {
  push @return_fields, "distrib";
} 

### expected occurrences
if ($query->param('return_exp')) {
  push @return_fields, "exp";
} 

$return_fields = join ",", @return_fields;


$parameters .= " -return $return_fields";


### single or both strands
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}


### group patterns by pairs of reverse complements
unless ($query->param('grouprc')) {
  $parameters .= " -nogrouprc";
} 

### prevent overlapping matches of the same pattern
if ($query->param('noov')) {
  $parameters .= " -noov";
} 

### verbose
$parameters .= " -v";


#### oligo length
$oligo_length = $query->param('oligo_length') ;
&FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
$parameters .= " -l $oligo_length";

#### class interval
$class_interval = $query->param('class_interval') ;
&FatalError("$class_interval Invalid class interval") unless &IsNatural($class_interval);
$parameters .= " -ci $class_interval";

################################################################
#### origin
$origin = $query->param('origin');
$parameters .= " -origin ".$origin;

################################################################
#### Offset
my $offset = $query->param('offset');
if ((&IsInteger($offset)) && ($offset != 0)) {
  $parameters .= " -offset ".$offset;
}

print "<PRE>command: ", &RSAT::util::hide_RSAT_path($command." ".$parameters), "<P>\n</PRE>" if ($ENV{rsat_echo} >=1);


if ($query->param('output') =~ /display/i) {

  &PipingWarning();

  ### execute the command ###
  $result_file = $tmp_file_path.".tab";
  push @result_files, ("Position analysis result", $result_file);
  open RESULT, "$command $parameters |";


  ### Print result on the web page
  print '<H2>Result</H2>';
  &PrintHtmlTable(RESULT, $result_file, true);
  close(RESULT);

  #### oligonucleotide assembly ####
  if (&IsReal($query->param('occ_significance_threshold'))) {

    ## Assemble the significant patterns with pattern-assembly
      $assembly_file = $tmp_file_path.".asmb";
#      $assembly_file = "$TMP/$tmp_file_name.asmb";
      push @result_files, ('assembly', $assembly_file);
    $pattern_assembly_command = $SCRIPTS."/pattern-assembly -v 1 -subst 0 -top 50";
    if ($query->param('strand') =~ /single/) {
      $pattern_assembly_command .= " -1str";
    } else {
      $pattern_assembly_command .= " -2str";
    }

    #     if (&IsNatural($query->param('max_asmb_nb'))) {
    #       $pattern_assembly_command .= " -max_asmb_nb ".$query->param('max_asmb_nb');
    #     }
    #     $pattern_assembly_command .= " -i ".$result_file;
    #     $pattern_assembly_command .= " -o ".$assembly_file;
    # 	$pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 1";
    # 	if ($query->param('strand') =~ /single/) {
    # 	    $pattern_assembly_command .= " -1str";
    # 	} else {
    # 	    $pattern_assembly_command .= " -2str";
    # 	}

    unless ($ENV{RSA_ERROR}) {
      print "<H2>Pattern assembly</H2>\n";
      open CLUSTERS, "$pattern_assembly_command  |";
      print "<PRE>\n";
      while (<CLUSTERS>) {
	s|$ENV{RSAT}/||g;
	print;
      }
      print "</PRE>\n";
      close(CLUSTERS);
    }
  }

  &PrintURLTable(@result_files);
  &PipingForm();
  print '<HR SIZE=3>';

} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $tmp_file_path);
}

print $query->end_html;

exit(0);


sub PipingForm {
    ### prepare data for piping
    
    #### title
    $title = $query->param('title');
    $title =~ s/\"/\'/g;

    #### strand for pattern-assembly
    if ($query->param('strand') =~ /single/) {
	$strand_opt .= " sensitive";
    } else {
	$strand_opt .= " insensitive";
    }
  print <<End_of_form;
<HR SIZE = 3>
<TABLE class = 'nextstep'>
<TR>

<TD colspan = 2>
<H3>Next step</H3>
</TD>
</tr><tr>
<TD>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="pattern matching (dna-pattern)">
</FORM>
</TD>

</TD>
<TD>
<FORM METHOD="POST" ACTION="pattern-assembly_form.cgi">
<INPUT type="hidden" NAME="local_pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="subst" VALUE=1>
<INPUT type="hidden" NAME="maxfl" VALUE=1>
<INPUT type="hidden" NAME="sc" VALUE="auto">
<INPUT type="hidden" NAME="strand" VALUE=$strand_opt>
<INPUT type="submit" value="pattern assembly">
</FORM>
</TD>


</TR>
</TABLE>
End_of_form

}



