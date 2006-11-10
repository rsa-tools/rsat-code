#!/usr/bin/perl
############################################################
#
# $Id: matrix-scan.cgi,v 1.3 2006/11/10 07:26:43 jvanheld Exp $
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
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$command = "$SCRIPTS/matrix-scan -v 1 ";
$tmp_file_name = sprintf "matrix-scan.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

$ECHO=1;

### print the result page
&RSA_header("matrix-scan result");
&ListParameters() if ($ECHO >=2);

#### update log file ####
&UpdateLogFile();

################################################################
## sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();

#### matrix-scan parameters
&ReadMatrixScanParameters();
$parameters .= " -i $sequence_file";

################################################################
### vertically print the matrix
if ($query->param('vertically_print')){
    $parameters .= " -p";
}

$command .= " $parameters";

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

    print "<PRE>";
    &PrintHtmlTable(RESULT, $result_file, true);
    print "</PRE>";

    if ($query->param('return_sites') eq 'on') {
	&PipingForm();
    }

    print "<HR SIZE = 3>";
    
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
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


################################################################
#
# read patser parameters
#
sub ReadMatrixScanParameters {

    ################################################################
    ## matrix specification
    unless ($query->param('matrix') =~ /\S/) { ### empty matrix
	&RSAT::error::FatalError("You did not enter the matrix");
    }
    $matrix_file = "$TMP/$tmp_file_name.matrix";

    $matrix_format = lc($query->param('matrix_format'));

    open MAT, "> ".$matrix_file;
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
    $parameters .= " -m $matrix_file";

    ################################################################
    ## pseudo-counts and weights are mutually exclusive
    if (&IsReal($query->param('pseudo_counts'))) {
      $parameters .= " -pseudo ".$query->param('pseudo_counts');
    }

    ################################################################
    ## strands 
    my $str = "-1str";
    if ($query->param('strands') =~ /both/i) {
      $str = "-2str";
      $parameters .= " -2str";
    } else {
      $str = "-1str";
      $parameters .= " -1str";
    }

    ################################################################
    #### origin
    if ($query->param('origin') =~ /end/i) {
	$parameters .= " -origin -0";
    }

    ################################################################
    ## Markov order
    my $markov_order = $query->param('markov_order');
    &RSAT::error::FatalError("Markov model should be a Natural number") unless &IsNatural($markov_order);

    ################################################################
    ## Background model method
    my $bg_method = $query->param('bg_method');
    if ($bg_method eq "bginput") {
      $parameters .= " -bginput";
      $parameters .= " -markov ".$markov_order;

    } elsif ($bg_method eq "window") {
      my $window_size = $query->param('window_size');
      &RSAT::message::FatalError(join("\t",$window_size, "Invalid value for the window size. Should be a Natural number." )) unless (&IsNatural($window_size));

      $parameters .= " -window ".$window_size;
      $parameters .= " -markov ".$markov_order;

    } elsif ($bg_method eq "bgfile") {

      my $organism_name = $query->param("organism");
      my $noov = "ovlp";
      my $background_model = $query->param("background");
      my $oligo_length = $markov_order + 1;
      $bg_file = &ExpectedFreqFile($organism_name,
				   $oligo_length, $background_model,
				   noov=>$noov, str=>$str);
      $parameters .= " -bgfile ".$bg_file;

    } else {
      &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
    }


    ################################################################
    ## Return fields
    @return_fields = qw(sites rank normw limits bg_model matrix freq_matrix weight_matrix);
    foreach my $field (@return_fields) {
      if ($query->param("return_".$field) eq "on") {
	$parameters .= " -return ".$field;
      }
    }

    ################################################################
    ## thresholds 
    @threshold_fields = qw(score rank proba_M proba_B normw);
    foreach my $field (@threshold_fields) {
      if ($query->param("lth_".$field) ne "none") {
	my $lth = $query->param("lth_".$field);
	&RSAT::error::FatalError($lth." is not a valid value for the lower $field threshold. Should be a number. ") unless (&IsReal($lth));
	$parameters .= " -lth $field $lth ";
      }

      if ($query->param("uth_".$field) ne "none") {
	my $uth = $query->param("uth_".$field);
	&RSAT::error::FatalError($uth." is not a valid value for the upper $field threshold. Should be a number. ") unless (&IsReal($uth));
	$parameters .= " -uth $field $uth ";

      }
    }

}


## Prepare data for piping
sub PipingForm {
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
