#!/usr/bin/perl
############################################################
#
# $Id: patser.cgi,v 1.7 2001/07/18 11:56:47 jvanheld Exp $
#
# Time-stamp: <2001-07-18 13:55:14 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";

$patser_command = "$BIN/patser";
$matrix_from_transfac_command = "$SCRIPTS/matrix-from-transfac";
$matrix_from_gibbs_command = "$SCRIPTS/matrix-from-gibbs";
$convert_seq_command = "$SCRIPTS/convert-seq";
$features_from_patser_cmd = "$SCRIPTS/features-from-patser";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
$add_yeast_link_command = "$SCRIPTS/add-yeast-link";
$tmp_file_name = sprintf "patser.%s", &AlphaDate;


### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("patser result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

################################################################
#
# read parameters
#

#### read parameters ####
$parameters = " -A a:t c:g ";

### matrix ####
unless ($query->param('matrix') =~ /\S/) { ### empty matrix
    &cgiError("You did not enter the matrix");
}

$matrix_file = "$TMP/$tmp_file_name.matrix";

$matrix_format = lc($query->param('matrix_format'));
if ($matrix_format =~ /transfac/i) {
    open MAT, "| $matrix_from_transfac_command > $matrix_file";
} elsif ($matrix_format =~ /gibbs/i) {
    open MAT, "| $matrix_from_gibbs_command > $matrix_file";
} elsif ($matrix_format =~ /consensus/i) {
    open MAT, "> $matrix_file";
} else {
    &cgiError("Invalid matrix format.");
}
print MAT $query->param('matrix');
close MAT;
&DelayedRemoval($matrix_file);
$parameters .= " -m $matrix_file";

### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1);
$parameters .= " -f $sequence_file";

### strands ###
if ($query->param('strands') =~ /both/i) {
    $parameters .= " -c";
}

### top value only ###
if ($query->param('return') =~ /top/i) {
    $parameters .= " -t";
}

### thresholds ###
if (&IsReal($query->param('lthreshold'))) {
    $parameters .= " -ls ".$query->param('lthreshold');
    $parameters .= " -M ".$query->param('lthreshold');
} 

if (&IsReal($query->param('uthreshold'))) {
    $parameters .= " -u ".$query->param('uthreshold');
}

### execute the command ###
if ($query->param('output') eq "display") {
    ### parameters for the piping to the feature map ###
    $feature_file =  "$TMP/$tmp_file_name.ft";
    $features_from_patser_cmd .= " -seq $sequence_file";
    $features_from_patser_cmd .= " -o $feature_file";

    print &PipingWarning();

    ### Print the result on Web page
    open RESULT, "$patser_command $parameters & |";
    open FEATURES, "| $features_from_patser_cmd";
    

    print "<PRE>";
    while (<RESULT>) {
	print;
	print FEATURES;
    }
    close FEATURES;
    close RESULT;
    print "</PRE>";
    print "<HR SIZE=3>\n";

    &PipingForm();
    
} else {
    ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your account: <P>";
	print "$address</B><P>";
	system "$patser_command $parameters | $mail_command $address &";
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


sub PipingForm {
    ### prepare data for piping
    print <<End_of_form;
<CENTER>
<TABLE>
<TR>
  <TD>
    <H4>Next step</H4>
  </TD>
  <TD>
    <FORM METHOD="POST" ACTION="feature-map_form.cgi">
    <INPUT type="hidden" NAME="feature_file" VALUE="$feature_file">
    <INPUT type="hidden" NAME="format" VALUE="feature-map">
    <INPUT type="hidden" NAME="handle" VALUE="none">
    <INPUT type="hidden" NAME="fill_form" VALUE="on">
    <INPUT type="submit" value="feature map">
    </FORM>
  </TD>
</TR>
</TABLE>
</CENTER>
End_of_form
}
