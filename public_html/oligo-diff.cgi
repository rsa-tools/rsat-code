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
require "RSA.disco.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("oligo-diff result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();

#### TEMPORARY
$command = "$SCRIPTS/oligo-diff";
$prefix = "oligo-diff";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "oligo-diff.%s", &AlphaDate();

#### read parameters ####
$parameters = " -v 1";

#### Test sequence set
$upload_test_seq = $query->param('upload_test_seq');
if ($upload_test_seq) {
  $tmp_test_seq = $tmp_file_path."_upload_test_seq.fasta";
  $upload_test_seq = $query->param('upload_test_seq');
  if ($upload_test_seq) {
    if ($upload_file =~ /\.gz$/) {
      $tmp_test_seq .= ".gz";
    }
    $type = $query->uploadInfo($upload_test_seq)->{'Content-Type'};
    open TEST_SEQ, ">$tmp_test_seq" ||
      &cgiError("Cannot store query feature file in temp dir.");
    while (<$upload_test_seq>) {
      print TEST_SEQ;
    }
    close TEST_SEQ;
  }
} elsif ($query->param('test_seq') =~/\S/) {
  $tmp_test_seq = $tmp_file_path."_pasted_test_seq.fasta";
  open TEST_SEQ, "> $tmp_test_seq";
  print TEST_SEQ $query->param('test_seq');
  close TEST_SEQ;
  &DelayedRemoval($tmp_test_seq);
} else {
  &FatalError ("Please select the test sequence file on your hard drive with the Browse button or paste sequences in the text area");
}
$parameters .= " -test $tmp_test_seq";
push @result_files, ("test sequences",$tmp_test_seq);

################################################################
#### Second sequence set

$upload_ctrl_seq = $query->param('upload_ctrl_seq');
if ($upload_ctrl_seq) {
    $tmp_ctrl_seq = $tmp_file_path."_upload_ctrl_seq.fasta";
    $upload_ctrl_seq = $query->param('upload_ctrl_seq');
    if ($upload_ctrl_seq) {
	if ($upload_file =~ /\.gz$/) {
	    $tmp_ctrl_seq .= ".gz";
	}
	$type = $query->uploadInfo($upload_ctrl_seq)->{'Content-Type'};
	open CTRL_SEQ, ">$tmp_ctrl_seq" ||
	  &cgiError("Cannot store expected frequency file in temp dir.");
	while (<$upload_ctrl_seq>) {
	    print CTRL_SEQ;
	}
	close CTRL_SEQ;
    }
} elsif ($query->param('ctrl_seq') =~/\S/) {
    $tmp_ctrl_seq = $tmp_file_path."_pasted_ctrl_seq.fasta";
    open CTRL_SEQ, "> $tmp_ctrl_seq";
    print CTRL_SEQ $query->param('ctrl_seq');
    close CTRL_SEQ;
    &DelayedRemoval($tmp_ctrl_seq);
} else {
    &FatalError ("Please select the second sequence file on your hard drive with the Browse button or paste sequences in the text area");
}
$parameters .= " -ctrl $tmp_ctrl_seq";
push @result_files, ("control sequences",$tmp_ctrl_seq);

### purge or not
if ($query->param('purge')) {
  $parameters .= " -purge";
} else {
  $parameters .= " -nopurge";
}

## oligo size
$oligo_len = $query->param('oligo_len') ;
&FatalError("$oligo_len Invalid oligonucleotide length") unless &IsNatural($oligo_len);
$parameters .= " -l $oligo_len";


### single or both strands
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}

### prevent overlapping matches of the same pattern
if ($query->param('noov')) {
  $parameters .= " -noov";
}

### Thresholds
foreach my $field qw(occ occ_sig occ_Pval occ_Eval) {
  if (&IsReal($query->param('lth_'.$field))) {
    $parameters .= " -lth ".$field;
    $parameters .= " ".$query->param('lth_'.$field);
  }
  if (&IsReal($query->param('uth_'.$field))) {
    $parameters .= " -uth ".$field;
    $parameters .= " ".$query->param('uth_'.$field);
  }
}

$command .= $parameters;

## Output file
$result_file = $tmp_file_path.".tab";
push @result_files, ("oligos", $result_file);

&ReportWebCommand($command);

if ($query->param('output') =~ /display/i) {
    &PipingWarning();

    ### Execute the command
    open RESULT, "$command |";

    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);

     #### oligonucleotide assembly ####
    if (($query->param('return') ne "table") &&
	($query->param('return') ne "distrib") &&
	(&IsReal($query->param('lth_occ_sig')))) {

      ## Pattern-assembly
      $assembly_file = $tmp_file_path.".asmb";
      $top_patterns = 50;
      $pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 1 -top ".$top_patterns;
      if ($query->param('strand') =~ /single/) {
	$pattern_assembly_command .= " -1str";
      } else {
	$pattern_assembly_command .= " -2str";
      }
      $pattern_assembly_command .= "  -i $result_file";
      $pattern_assembly_command .= "  -o $assembly_file";

      unless ($ENV{RSA_ERROR}) {

	## Assemble the significant patterns
	print "<H2>Pattern assembly</H2>\n";
	print "<PRE>pattern-assembly command: $pattern_assembly_command<P>\n</PRE>" if ($ENV{rsat_echo} >=1);
	system "$pattern_assembly_command";
	open ASSEMBLY, $assembly_file;
	print "<PRE>\n";
	while (<ASSEMBLY>) {
	  s|$ENV{RSAT}/||g;
	  print;
	}
	print "</PRE>\n";
	close(ASSEMBLY);
	push @result_files, ('assembly', $assembly_file);


	## Convert pattern-assembly result into PSSM
	if ($query->param('to_matrix')) {
	  local $sequence_file = $tmp_test_seq;
	  local $sequence_format = "fasta";
	  &MatrixFromPatterns_run();
	}
      }
    }

    &PrintURLTable(@result_files);
    &OligoDyadPipingForm();
    print '<HR SIZE=3>';

} else {
    &EmailTheResult("$command", $query->param('user_email'), $result_file);
}

print $query->end_html;

exit(0);


