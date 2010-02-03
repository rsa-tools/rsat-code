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

$command = "$SCRIPTS/oligo-diff";
$tmp_file_name = sprintf "oligo-diff.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("oligo-diff result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();

#### read parameters ####
$parameters = " -v 1";

#### First sequence set
$upload_seq1 = $query->param('upload_seq1');
if ($upload_seq1) {
    $tmp_seq1 = "${TMP}/${tmp_file_name}_upload_seq1.tab";
    $upload_seq1 = $query->param('upload_seq1');
    if ($upload_seq1) {
	if ($upload_file =~ /\.gz$/) {
	    $tmp_seq1 .= ".gz";
	}
	$type = $query->uploadInfo($upload_seq1)->{'Content-Type'};
	open SEQ1, ">$tmp_seq1" ||
	  &cgiError("Cannot store query feature file in temp dir.");
	while (<$upload_seq1>) {
	    print SEQ1;
	}
	close SEQ1;
    }
}elsif ($query->param('seq1') =~/\S/) {
    $tmp_seq1 = "${TMP}/${tmp_file_name}_pasted_seq1.tab";
    open SEQ1, "> $tmp_seq1";
    print SEQ1 $query->param('seq1');
    close SEQ1;
    &DelayedRemoval($tmp_seq1);
}else {
    &FatalError ("Please select the first sequence file on your hard drive with the Browse button or paste sequences in the text area");
}
$parameters .= " -file1 $tmp_seq1";

################################################################
#### Second sequence set

$upload_seq2 = $query->param('upload_seq2');
if ($upload_seq2) {
    $tmp_seq2 = "${TMP}/${tmp_file_name}_upload_seq2.tab";
    $upload_seq2 = $query->param('upload_seq2');
    if ($upload_seq2) {
	if ($upload_file =~ /\.gz$/) {
	    $tmp_seq2 .= ".gz";
	}
	$type = $query->uploadInfo($upload_seq2)->{'Content-Type'};
	open SEQ2, ">$tmp_seq2" ||
	  &cgiError("Cannot store expected frequency file in temp dir.");
	while (<$upload_seq2>) {
	    print SEQ2;
	}
	close SEQ2;
    }
} elsif ($query->param('seq2') =~/\S/) {
    $tmp_seq2 = "${TMP}/${tmp_file_name}_pasted_seq2.tab";
    open SEQ2, "> $tmp_seq2";
    print SEQ2 $query->param('seq2');
    close SEQ2;
    &DelayedRemoval($tmp_seq2);
} else {
    &FatalError ("Please select the second sequence file on your hard drive with the Browse button or paste sequences in the text area");
}
$parameters .= " -file2 $tmp_seq2";

### purge or not
if ($query->param('purge')) {
  $parameters .= " -purge";
} else {
  $parameters .= " -nopurge";
}

## oligo size
$oligo_length = $query->param('oligo_length') ;
&FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
$parameters .= " -l $oligo_length";


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




print "<PRE>command: $command $return_fields $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >=1);

if ($query->param('output') =~ /display/i) {

#    &PipingWarning();

    ### execute the command ###
    $result_file = "$TMP/${tmp_file_name}.res";
    open RESULT, "$command $parameters $return_fields |";

    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

#    &PipingForm();
    print '<HR SIZE=3>';

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters $return_fields", $query->param('user_email'), $tmp_file_name);
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $tmp_file_name);
}

print $query->end_html;

exit(0);


