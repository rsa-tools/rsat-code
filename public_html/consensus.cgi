#!/usr/bin/perl
############################################################
#
# $Id: consensus.cgi,v 1.7 2003/10/29 09:04:36 jvanheld Exp $
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
$matrix_from_consensus_command = "$SCRIPTS/matrix-from-consensus";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "consensus.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("consensus result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

################################################################
#
# read parameters
#

### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1, 0);
$parameters .= " -f $sequence_file ";




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
if (&IsNatural($query->param('repeats'))) {
    if ($query->param('one_per_seq')) {
	$parameters .= " -N ".$query->param('repeats');
    } else {
	$parameters .= " -n ".$query->param('repeats');
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


#### error file
#$error_file  = "$TMP/$tmp_file_name.err";

#### execute the command ###
if ($query->param('output') eq "display") {

    
    $result_file = "$TMP/$tmp_file_name.res";
    $matrix_file = "$TMP/$tmp_file_name.matrix";

    #### echo the command
    print "<PRE><B>Command:</B> $command $parameters </PRE>" if ($ECHO >= 1);
    print "<PRE>$matrix_from_consensus_command -i $result_file -o $matrix_file</PRE>" if ($ECHO >= 1);

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
    
    system "$matrix_from_consensus_command -i $result_file -o $matrix_file";

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
</TR>
</TABLE>
End_of_form
  
}
