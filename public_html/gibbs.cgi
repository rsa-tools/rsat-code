#!/usr/bin/perl
############################################################
#
# $Id: gibbs.cgi,v 1.12 2004/06/10 04:42:45 jvanheld Exp $
#
# Time-stamp: <2003-05-13 11:30:48 jvanheld>
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

$command = "$BIN/gibbs";
#$convert_matrix_command = "$SCRIPTS/matrix-from-gibbs";
$convert_matrix_command = "$SCRIPTS/convert-matrix -format gibbs -return alignment";
$convert_seq_command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "gibbs.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("gibbs result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

################################################################
#
# read parameters
#

#### add reverse complement
if (lc($query->param("add_rc")) eq "on") {
    $add_rc = 1;
    $convert_seq_options .= "-addrc ";
} else {
    $add_rc = 0;
}

### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("fasta", 1, $add_rc);
$parameters = " $sequence_file ";

### matrix length
if (&IsNatural($query->param('length'))) {
    $parameters .= $query->param('length')." ";
}

### expected number of matches
if (&IsNatural($query->param('expected'))) {
    $parameters .= $query->param('expected')." ";
}

### sequence type
if (lc($query->param('seq_type')) eq "dna") {
    $parameters .= "-n ";
}

### inactivate frqgmentation
unless (lc($query->param('fragmentation')) eq "on") {
    $parameters .= "-d ";
}

if ($query->param('output') eq "display") {  
    &PipingWarning();

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    $matrix_file = "$TMP/$tmp_file_name.matrix";
    
    system "$command $parameters > $result_file";

    $convert_matrix_command .= " -i ".$result_file." -o ".$matrix_file;
    system "$convert_matrix_command";
    if ($ECHO >= 1) {
	print "<PRE><B>Command:</B> $command $parameters </PRE>";
	print "<PRE><B>Conversion:</B> $convert_matrix_command </PRE>";
    }
    
    ### Print result on the web page
    print '<H4>Result</H4>';
    print "<PRE>";
    print `cat $result_file`;
    print "</PRE>";
    
    &PipingForm();
  
    &DelayedRemoval($result_file);
    &DelayedRemoval($matrix_file);

    print "<HR SIZE = 3>";
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

print $query->end_html;
exit(0);



### prepare data for piping
sub PipingForm {
  $title = $query->param('title');
  $title =~ s/\"/\'/g;
    print <<End_of_form;
<HR SIZE = 3>
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
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="pattern matching (patser)">
</FORM>
</TD>
<TD>
<FORM METHOD="POST" ACTION="convert-matrix_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$result_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="gibbs">
<INPUT type="submit" value="convert-matrix">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form
  
}
