#!/usr/bin/perl
############################################################
#
# $Id: purge-sequence.cgi,v 1.1 2002/01/07 02:03:34 jvanheld Exp $
#
# Time-stamp: <2001-10-08 00:51:11 jvanheld>
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

$gibbs_command = "nice -n 30 $BIN/gibbs";
$matrix_from_gibbs_command = "$SCRIPTS/matrix-from-gibbs";
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
    &DelayedRemoval($result_file);
    &DelayedRemoval($matrix_file);
    
    system "$gibbs_command $parameters > $result_file";
    system "$matrix_from_gibbs_command -i $result_file -o $matrix_file";
    print "<PRE><B>Command:</B> $gibbs_command $parameters </PRE>";
    
    ### Print result on the web page
    print '<H4>Result</H4>';
    print "<PRE>";
    print `cat $result_file`;
    print "</PRE>";
    
    &PipingForm();
  
} else {
  #### send e-mail with the result
  if ($query->param('user_email') =~ /(\S+\@\S+)/) {
    $address = $1;
    print "<B>Result will be sent to your e-mail address: <P>";
    print "$address</B><P>";
    system "$gibbs_command $parameters | $mail_command $address &"; 
  } else {
    if ($query->param('user_email') eq "") {
      print "<B>ERROR: you did not enter your e-mail address<P>";
    } else {
      print "<B>ERROR: the e-mail address you entered is not valid<P>";
      print "$query->param('user_email')</B><P>";      
    }
    }
  print '<HR SIZE=3>';
}

print "<HR SIZE = 3>";
print $query->end_html;
exit(0);



sub PipingForm {
  ### prepare data for piping
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
</TR>
</TABLE>
End_of_form
  
}
