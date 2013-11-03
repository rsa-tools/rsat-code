#!/usr/bin/perl
############################################################
#
# $Id: consensus.cgi,v 1.25 2013/11/03 19:33:31 jvanheld Exp $
#
# Time-stamp: <2003-07-03 10:06:42 jvanheld>
#
############################################################
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


### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("Consensus result file", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$command = $BIN."/consensus";
#$convert_matrix_command = "$SCRIPTS/matrix-from-consensus -v 1";
$convert_matrix_command = "$SCRIPTS/convert-matrix -from consensus -return counts";
$convert_seq_command = "$SCRIPTS/convert-seq";
#$tmp_file_name = sprintf "consensus.%s", &AlphaDate();
$prefix = "consensus";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
@result_files = ();

################################################################
#
# Read parameters
#

### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", no_format=>1, add_rc=>0);
push @result_files, "Input sequences ($sequence_format)",$sequence_file;

$parameters .= " -f $sequence_file ";


## Number of matrices to save
if (&IsNatural($query->param('matrices_to_save'))) {
    $parameters .= " -pt ".$query->param('matrices_to_save');
    $parameters .= " -pf ".$query->param('matrices_to_save');
}

### strands and pattern symmetry 
if ($query->param('symmetrical')) {
    $parameters .= " -c3 ";
} elsif ($query->param('strands') =~ /ignore/i) {
    $parameters .= " -c0 ";
} elsif ($query->param('strands') =~ /separate/i) {
    $parameters .= " -c1 ";
} elsif ($query->param('strands') =~ /single/i) {
    $parameters .= " -c2 ";
}

### pattern length ###
if (&IsNatural($query->param('length'))) {
    $parameters .= " -L ".$query->param('length');
}

### alphabet ###
$parameters .= " -A ".$query->param('alphabet');

### use designated prior frequencies ###
if ($query->param('prior_freq') eq "on") {
    $parameters .= " -d ";
}

$seed_ok = 1; #### for option compatibility

### expected matches ###
if (&IsNatural($query->param('cycles'))) {
    if ($query->param('one_per_seq')) {
	$parameters .= " -N ".$query->param('cycles');
    } else {
	$parameters .= " -n ".$query->param('cycles');
	$seed_ok = 0;
    }
}

### seed with first sequence and proceed linearly ###
if ($query->param('seed') eq "on") {
    if ($seed_ok) {
	$parameters .= " -l ";
    } else {
	&cgiWarning("Seed option will be ignored because it requires at least one match within each sequence. ");
    }
}

#### Output and input files
$result_file = $tmp_file_path.".txt";
push @result_files, ('Consensus result', $result_file);

$matrix_file = $tmp_file_path.".tab";
push @result_files, "Output matrix (tab)", $matrix_file;

#### Matrix conversion command
$convert_matrix_command.= " -i ".$result_file." -o ".$matrix_file;

#### error file
#$error_file  = "$TMP/$tmp_file_name.err";

#### execute the command ###
if ($query->param('output') eq "display") {
    #### echo the commands
    if ($ENV{rsat_echo} >= 1) {
      &ReportWebCommand($command." ".$parameters);
      &ReportWebCommand($convert_matrix_command);
    }

    ### prepare data for piping
    &PipingWarning();

    open RESULT, "$command $parameters | ";
    open RES_FILE, ">$result_file";
  #    system "$command $parameters >$result_file ";

    ### Print result on the web page
    print '<H4>Result</H4>';
    print "<PRE>";
    while (<RESULT>) {
	print $_;
	print RES_FILE $_;
    }
    print "</PRE>";
    close(RESULT);
    close(RES_FILE);


    ## Display matrices with logos
    my ($out_matrix_file) = &display_matrices_web($result_file, "consensus");
#    push @result_files, ("Converted matrix", $out_matrix_file);

    &PrintURLTable(@result_files);
    &PipingForm();

    system "$convert_matrix_command -i $result_file -o $matrix_file";

    &DelayedRemoval($result_file);
    &DelayedRemoval($matrix_file);

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}

print $query->end_html;

exit(0);




################################################################
## Piping form: to send the result as input for another tool

sub PipingForm {
    print <<End_of_form;
<TABLE class = 'nextstep'  CELLSPACING=0 CELLPADDING=10 BORDER=0 NOWRAP>

<TR VALIGN="top" ALIGN="center">
    <Th VALIGN=BOTTOM ALIGN=CENTER COLSPAN=3>
	Next step
    </Th>

</TR>

<TD VALIGN=BOTTOM ALIGN=CENTER>
<FORM METHOD="POST" ACTION="patser_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$matrix_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="consensus">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="wconsensus">
<INPUT type="submit" value="pattern matching (patser)">
</FORM>
</TD>

<TD VALIGN=BOTTOM ALIGN=CENTER>
<b><font color=red>New</a></b>
<FORM METHOD="POST" ACTION="matrix-scan_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$matrix_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="consensus">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="wconsensus">
<INPUT type="submit" value="pattern matching (matrix-scan)">
</FORM>
</TD>

<TD VALIGN=BOTTOM ALIGN=CENTER>
<FORM METHOD="POST" ACTION="convert-matrix_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$result_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="consensus">
<INPUT type="submit" value="convert-matrix">
</FORM>
</TD>

</TR>

</TABLE>
End_of_form
  
}
