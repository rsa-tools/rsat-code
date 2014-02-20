#!/usr/bin/perl
############################################################
#
# $Id: patser.cgi,v 1.33 2013/11/03 19:33:31 jvanheld Exp $
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
require "RSA2.cgi.lib";
require "patser.lib.pl";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
### Read the CGI query
$query = new CGI;

#$ENV{rsat_echo}=2;

### print the result page
&RSA_header("patser result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();

$command = $BIN."/patser";
#$convert_seq_command = $SCRIPTS."/convert-seq";
$features_from_patser_cmd = $SCRIPTS."/features-from-patser -v 1";
$add_yeast_link_command = $SCRIPTS."/add-yeast-link";
$prefix = "patser";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

################################################################
#### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", no_format=>1,add_rc=>0, skip_short=>30);

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

## Output file
$result_file =  $tmp_file_path.".ft";
push @result_files, "Patser result (ft)", $result_file;

################################################################
#### echo the command (for debugging)
&ReportWebCommand($command);

################################################################
### execute the command ###
if ($query->param('output') eq "display") {

    unless ($query->param('table')) {
	&PipingWarning();
    }

    ### Print the result on Web page
    open RESULT, "$command |";
#    open FEATURES, "| $features_from_patser_cmd";
    



    print "<PRE>";
    &PrintHtmlTable(RESULT, $result_file, true);
    
#    while (<RESULT>) {
#	s|$ENV{RSAT}/||g;
#	print;
#	print FEATURES;
#    }
#    close FEATURES;
#    close RESULT;
    print "</PRE>";

    &PrintURLTable(@result_files);
    if ($query->param('output_format') eq 'feature list') {
	&PipingForm();
    }
    
    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $patser_parameters", $query->param('user_email'), $result_file);
}


print $query->end_html;

exit(0);


## Prepare data for piping
sub PipingForm {
    print <<End_of_form;
<CENTER>
<TABLE class = 'nextstep'>
<TR>
  <TD>
    <H3>Next step</H3>
  </TD>
  </tr>
  <tr>
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
