#!/usr/bin/perl
############################################################
#
# $Id: purge-sequence.cgi,v 1.4 2003/10/29 09:04:37 jvanheld Exp $
#
# Time-stamp: <2003-10-01 00:38:45 jvanheld>
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

$command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "purge-sequence.%s", &AlphaDate;
$out_format = "fasta";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("purge-sequence result");
&ListParameters() if ($ECHO >= 2);

#### update log file ####
&UpdateLogFile;

################################################################
#
# read parameters
#

### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("fasta", 1, $add_rc);
$parameters = " -i $sequence_file ";

$result_file = "$TMP/$tmp_file_name.res";
#$parameters .= " -o $result_file ";
#&DelayedRemoval($result_file);
    

#### add reverse complement
if ($query->param("both_strands")) {
    $parameters .= "-2str ";
} else {
    $parameters .= "-1str ";
}

### treatment
if ($query->param('treatment') eq 'delete') {
    $parameters .= " -del ";
}

### match length
if (&IsNatural($query->param('match_len'))) {
    $parameters .= " -ml ".$query->param('match_len');
}

### mismatches
if (&IsNatural($query->param('mismatches'))) {
    $parameters .= " -mis ".$query->param('mismatches');
}

#$parameters .= " -o $TMP/bounpurged";

print "<PRE><B>Command:</B> $command $parameters </PRE>" if ($ECHO >= 1);

if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

    ### execute the command ###
    open RESULT, "$command $parameters |";

    ### print the result ### 
    &PipingWarning();
    if ($query->param('output') =~ /server/i) {
	&Info("The server is currently retrieving your sequences...");
    }

    print '<H2>Result</H2>';
    
    ### open the mirror file ###
    $mirror_file = "$TMP/$tmp_file_name.res";
    if (open MIRROR, ">$mirror_file") {
	$mirror = 1;
	&DelayedRemoval($mirror_file);
    }
    
    ### Print result on the web page
    print "<PRE>";
    while (<RESULT>) {
	print "$_" unless ($query->param('output') =~ /server/i);
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close(RESULT);
    close MIRROR if ($mirror);
    
    if ($query->param('output') =~ /server/i) {
	$result_URL = "${WWW_RSA}/tmp/${tmp_file_name}.res";
	print ("Result is available in the file ",
	       "<a href=${result_URL}>${result_URL}</a>",
	       "\n");
    }
    
    ### prepare data for piping
    &PipingForm();
    
    print "<HR SIZE = 3>";

#} elsif ($query->param('output') =~ /server/i) {
#    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

print $query->end_html;
exit(0);


################################################################
#
# Pipe the result to other programs
#
sub PipingForm {
  print <<End_of_form;
<HR SIZE = 3>
<TABLE CELLSPACING=0 CELLPADDING=10 BORDER=0 NOWRAP BGCOLOR= #FFEEDD>

<TR VALIGN="top" ALIGN="center">
    <TD COLSPAN=5 BGCOLOR=	#FFEEDD>
	<H3>Next step</H3>
    </TD>

</TR>

<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=		#FFEEDD>
	<B>Pattern discovery</B><BR>
	(unknown patterns)
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="oligo-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$result_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="oligonucleotide analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dyad-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$result_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dyad analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="consensus_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$result_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="consensus">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="gibbs_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$result_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="gibbs sampler">
	</FORM>
    </TD>

</TR>

<!--
<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=		#FFEEDD>
	<B>Pattern matching</B><BR>
	(known patterns)
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$result_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dna-pattern (IUPAC)">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="patser_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$result_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="patser (matrices)";
	</FORM>
    </TD>

    <TD>
    &nbsp;
    </TD>

    <TD>
    &nbsp;
    </TD>

</TR>
-->

</TABLE>
End_of_form
}
