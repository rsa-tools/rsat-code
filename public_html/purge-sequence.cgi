#!/usr/bin/perl
############################################################
#
# $Id: purge-sequence.cgi,v 1.2 2002/09/05 10:06:38 jvanheld Exp $
#
# Time-stamp: <2002-09-05 05:05:56 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "purge-sequence.%s", &AlphaDate;
$out_format = "fasta";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("purge-sequence result");
#&ListParameters;

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


### match length
if (&IsNatural($query->param('match_len'))) {
    $parameters .= " -ml ".$query->param('match_len');
}

### mismatches
if (&IsNatural($query->param('mismatches'))) {
    $parameters .= " -mis ".$query->param('mismatches');
}

print "<PRE><B>Command:</B> $purge_sequence_command $parameters </PRE>" if ($ECHO);

if ($query->param('output') eq "display") {  
    &PipingWarning();
    
    ### execute the command ###
#    system "$purge_sequence_command $parameters";
    open RESULT, "$purge_sequence_command $parameters |";
    
    ### open the mirror file ###
    $mirror_file = "$TMP/$tmp_file_name.res";
    if (open MIRROR, ">$mirror_file") {
	$mirror = 1;
	&DelayedRemoval($mirror_file);
    }
    
    ### Print result on the web page
    print '<H2>Result</H2>';
    print "<PRE>";
    while (<RESULT>) {
  	print;
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close(RESULT);
    close MIRROR if ($mirror);
    
#    ### Print result on the web page
#    print '<H4>Result</H4>';
#    print "<PRE>";
#    print `cat $result_file`;
#    print "</PRE>";
    
    &PipingForm();
    
    print "<HR SIZE = 3>";

} else {

    &EmailTheResult("$purge_sequence_command $parameters", $query->param('user_email'));

#   #### send e-mail with the result
#   if ($query->param('user_email') =~ /(\S+\@\S+)/) {
#     $address = $1;
#     print "<B>Result will be sent to your e-mail address: <P>";
#     print "$address</B><P>";
#     system "$purge_sequence_command $parameters | $mail_command $address &"; 
#   } else {
#     if ($query->param('user_email') eq "") {
#       print "<B>ERROR: you did not enter your e-mail address<P>";
#     } else {
#       print "<B>ERROR: the e-mail address you entered is not valid<P>";
#       print "$query->param('user_email')</B><P>";      
#     }
#     }
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
