#!/usr/bin/perl
############################################################
#
# $Id: consensus.cgi,v 1.3 2000/12/26 22:51:48 jvanheld Exp $
#
# Time-stamp: <2000-12-26 23:51:42 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

$consensus_command = "$BIN/consensus";
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

### strands ###
if ($query->param('strands') =~ /ignore/i) {
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

  ### seed with first sequence and proceed linearly ###
if ($query->param('seed') eq "on") {
    $parameters .= " -l ";
} elsif (&IsNatural($query->param('repeats'))) {
    ### expected matches ###
    if ($query->param('one_per_seq') eq "on") {
	  $parameters .= " -n ".$query->param('repeats');
      } else {
	  $parameters .= " -N ".$query->param('repeats');
      }
}

### result file
#  $result_file = "$TMP/$tmp_file_name.res";
#  $matrix_file = "$TMP/$tmp_file_name.matrix";
$error_file  = "$TMP/$tmp_file_name.err";


### print the header
print <<End_Header;
<HEADER>
<TITLE>CONSENSUS result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER>Matrix extraction (consensus) result $query->param('set_name')</H3>
End_Header



    ### execute the command ###
if ($query->param('output') eq "display") {
    
    $result_file = "$TMP/$tmp_file_name.res";
    $matrix_file = "$TMP/$tmp_file_name.matrix";
#    print "<PRE>";
#    print "$consensus_command $parameters | ", "\n";
#    print "$matrix_from_consensus_command -i $result_file -o $matrix_file";
#    print "</PRE>";
    open RESULT, "$consensus_command $parameters | ";
    open RES_FILE, ">$result_file";
    
    #print "<PRE><B>Command:</B> $consensus_command $parameters </PRE>";
    
    ### prepare data for piping
#	$title = $query->param('title');
#	$title =~ s/\"/\'/g;
    print <<End_of_form;
<TABLE>
<TR>
<TD>
<H4>Next step</H4>
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
</TR>
</TABLE>
End_of_form
  
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
    
    
    system "$matrix_from_consensus_command -i $result_file -o $matrix_file";
    
} else {
    ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your e-mail address: <P>";
	print "$address</B><P>";
	system "$consensus_command $parameters | $mail_command $address &";
    } else {
	if ($query->param('user_email') eq "") {
	    print "<B>ERROR: you did not enter your e-mail address<P>";
	} else {
	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
	    print $query->param('user_email')."</B><P>";      
	}
    } 
}

print "<HR SIZE = 3>";
print $query->end_html;
exit(0);






