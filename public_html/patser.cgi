#!/usr/bin/perl
############################################################
#
# $Id: patser.cgi,v 1.20 2003/10/29 09:04:37 jvanheld Exp $
#
# Time-stamp: <2003-06-16 00:59:07 jvanheld>
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
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = "$BIN/patser";
$matrix_from_transfac_command = "$SCRIPTS/matrix-from-transfac";
$matrix_from_gibbs_command = "$SCRIPTS/matrix-from-gibbs";
$convert_seq_command = "$SCRIPTS/convert-seq";
$features_from_patser_cmd = "$SCRIPTS/features-from-patser";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
$add_yeast_link_command = "$SCRIPTS/add-yeast-link";
$tmp_file_name = sprintf "patser.%s", &AlphaDate;


### Read the CGI query
$query = new CGI;

#$ECHO=2;

### print the result page
&RSA_header("patser result");
&ListParameters() if ($ECHO >=2);

#### update log file ####
&UpdateLogFile();


################################################################
#### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1);

#### patser parameters
&ReadPatserParameters();
$patser_parameters .= " -f $sequence_file";

#### feature-from-patser parameters
&ReadFeaturesFromPatserParams(); 


################################################################
### vertically print the matrix
if ($query->param('vertically_print')){
    $patser_parameters .= " -p";
}

$command .= " $patser_parameters";
$command .= "| $features_from_patser_cmd";

################################################################
#### echo the command (for debugging)
print "<pre>$command</pre>" if ($ECHO >= 1);

################################################################
### execute the command ###
if ($query->param('output') eq "display") {

    unless ($query->param('table')) {
	&PipingWarning();
    }

    ### Print the result on Web page
    $result_file =  "$TMP/$tmp_file_name.ft";
    open RESULT, "$command |";
#    open FEATURES, "| $features_from_patser_cmd";
    



    print "<PRE>";
    &PrintHtmlTable(RESULT, $result_file, true);
    
#    while (<RESULT>) {
#	s|$RSA/||g;
#	print;
#	print FEATURES;
#    }
#    close FEATURES;
#    close RESULT;
    print "</PRE>";

    if ($query->param('output_format') eq 'feature list') {
	&PipingForm();
    }

    print "<HR SIZE = 3>";
    
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $patser_parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $patser_parameters", $query->param('user_email'));
}


print $query->end_html;

exit(0);


sub PipingForm {
    ### prepare data for piping
    print <<End_of_form;
<CENTER>
<TABLE>
<TR>
  <TD>
    <H3>Next step</H3>
  </TD>
  <TD>
    <FORM METHOD="POST" ACTION="feature-map_form.cgi">
    <INPUT type="hidden" NAME="feature_file" VALUE="$result_file">
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
