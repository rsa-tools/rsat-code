#!/usr/bin/perl
############################################################
#
# $Id: consensus.cgi,v 1.10 2004/06/10 04:42:45 jvanheld Exp $
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
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = "$BIN/consensus";
#$convert_matrix_command = "$SCRIPTS/matrix-from-consensus -v 1";
$convert_matrix_command = "$SCRIPTS/convert-matrix -format consensus -return alignment";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "consensus.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("consensus result");
&ListParameters() if ($ECHO >= 2);

#### update log file ####
&UpdateLogFile();

################################################################
#
# read parameters
#

### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1, 0);
$parameters .= " -f $sequence_file ";

## Number of matrices to save
if (&IsNatural($query->param('matrices_to_save'))) {
    $parameters .= " -q ".$query->param('matrices_to_save');
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

#### Output files    
$result_file = "$TMP/$tmp_file_name.res";
$matrix_file = "$TMP/$tmp_file_name.matrix";

#### Matrix conversion command
$convert_matrix_command.= " -i".$result_file." -o ".$matrix_file;

#### error file
#$error_file  = "$TMP/$tmp_file_name.err";

#### execute the command ###
if ($query->param('output') eq "display") {
    #### echo the commands
    if ($ECHO >= 1) {
	print "<PRE><B>Command:</B> $command $parameters </PRE>";
	print "<PRE><b>Conversion:</b> $convert_matrix_command</PRE>";
    }

    open RESULT, "$command $parameters | ";
    open RES_FILE, ">$result_file";
    
    ### prepare data for piping
    &PipingWarning();
    
    ### Print result on the web page
    print '<H3>Result</H3>';
    print "<PRE>";
    while (<RESULT>) {
	print $_;
	print RES_FILE $_;
    }
    print "</PRE>";
    close(RESULT);
    close(RES_FILE);
    
    &PipingForm();
    
    system "$convert_matrix_command -i $result_file -o $matrix_file";

    print "<HR SIZE = 3>";
    
    &DelayedRemoval($result_file);
    &DelayedRemoval($matrix_file);

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'), $tmp_file_name);
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $tmp_file_name);
}

print $query->end_html;

exit(0);




################################################################
## Piping form: to send the result as input for another tool

sub PipingForm {
    print <<End_of_form;
<TABLE>

<TR>
<TD>
<H3>Next step</H3>
</TD>
<TD>
<FORM METHOD="POST" ACTION="patser_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$matrix_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="consensus">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="wconsensus">
<INPUT type="submit" value="pattern matching (patser)">
</FORM>
</TD>
<TD>
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
