#!/usr/bin/perl
############################################################
#
# $Id: purge-sequence.cgi,v 1.17 2013/11/03 19:33:31 jvanheld Exp $
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
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("purge-sequence result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$command = "$SCRIPTS/purge-sequence";
$prefix = "purge-sequence";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
#$tmp_file_name = sprintf "purge-sequence.%s", &AlphaDate;
$out_format = "fasta";
@result_files = ();


################################################################
#
# read parameters
#

### Input sequence file ####
($in_sequence_file,$sequence_format) = &GetSequenceFile("fasta", no_format=>1, add_rc=>0);
$parameters = " -i $in_sequence_file ";
push @result_files, ("Input sequences (fasta)",$in_sequence_file);

## Output file (purged sequences)
$sequence_file = $tmp_file_path.".fasta";
push @result_files, ("Purged sequences (fasta)",$sequence_file);

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

&ReportWebCommand($command." ".$parameters);

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
    if (open MIRROR, ">$sequence_file") {
	$mirror = 1;
	&DelayedRemoval($sequence_file);
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

    ### prepare data for piping
    &PrintURLTable(@result_files);
    &PipingFormForSequence();
    print "<HR SIZE = 3>";
#    ### prepare data for piping
#    &PipingForm();

    print "<hr size = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $sequence_file);
}

print $query->end_html;
exit(0);


# ################################################################
# #
# # Pipe the result to other programs
# #
# sub PipingForm {
#   print <<End_of_form;
# <HR SIZE = 3>
# <TABLE class='nextstep'>

# <TR VALIGN="top" ALIGN="center">
#     <TD COLSPAN=5>
# 	<H3>Next step</H3>
#     </TD>

# </TR>

# <TR VALIGN="top" ALIGN="center">

#     <TD>
# 	<B>motif discovery</B><BR>
# 	(unknown patterns)
#     </TD>

#     <TD>
# 	<FORM METHOD="POST" ACTION="oligo-analysis_form.cgi">
# 	<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
# 	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
# 	<INPUT type="submit" value="oligonucleotide analysis">
# 	</FORM>
#     </TD>

#     <TD>
# 	<FORM METHOD="POST" ACTION="dyad-analysis_form.cgi">
# 	<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
# 	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
# 	<INPUT type="submit" value="dyad analysis">
# 	</FORM>
#     </TD>

#     <TD>
# 	<FORM METHOD="POST" ACTION="consensus_form.cgi">
# 	<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
# 	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
# 	<INPUT type="submit" value="consensus">
# 	</FORM>
#     </TD>

#     <TD>
# 	<FORM METHOD="POST" ACTION="gibbs_form.cgi">
# 	<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
# 	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
# 	<INPUT type="submit" value="gibbs sampler">
# 	</FORM>
#     </TD>

# </TR>

# <!--
# <TR VALIGN="top" ALIGN="center">

#     <TD BGCOLOR=		#FFEEDD>
# 	<B>Pattern matching</B><BR>
# 	(known patterns)
#     </TD>

#     <TD>
# 	<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
# 	<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
# 	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
# 	<INPUT type="submit" value="dna-pattern (IUPAC)">
# 	</FORM>
#     </TD>

#     <TD VALIGN=BOTTOM ALIGN=CENTER>
#         <b><font color=red>New</a></b>
#         <FORM METHOD="POST" ACTION="matrix-scan_form.cgi">
#         <INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
#         <INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
#         <INPUT type="submit" value="matrix-scan (matrices)">
#         </FORM>
#     </TD>


#     <TD>
# 	<FORM METHOD="POST" ACTION="patser_form.cgi">
# 	<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
# 	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
# 	<INPUT type="submit" value="patser (matrices)";
# 	</FORM>
#     </TD>

#     <TD>
#     &nbsp;
#     </TD>

#     <TD>
#     &nbsp;
#     </TD>

# </TR>
# -->

# </TABLE>
# End_of_form
# }
