#!/usr/bin/perl
############################################################
#
# $Id: patser.cgi,v 1.15 2003/04/17 20:45:42 jvanheld Exp $
#
# Time-stamp: <2003-04-17 22:45:12 jvanheld>
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

### print the result page
&RSA_header("patser result");
&ListParameters if ($ECHO >=2);

#### update log file ####
&UpdateLogFile;

################################################################
#
# read parameters
#

#### read parameters ####

################################################################
### alphabet ###
$parameters .= " -A ".$query->param('alphabet');


################################################################
### matrix specification
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

#### [-w <Matrix is a weight matrix>]
if ($query->param('matrix_is_weight')) {
    $parameters .= " -w";
} else {
    #### pseudo-counts and weights are mutually exclusive
    if (&IsReal($query->param('pseudo_counts'))) {
	$parameters .= " -b ".$query->param('pseudo_counts');
    }

}

#### [-v <Vertical matrix---rows correspond to positions>]
if ($query->param('matrix_is_vertical')) {
    $parameters .= " -v";
}


################################################################
### sequence file ####
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1);
$parameters .= " -f $sequence_file";

################################################################
### strands ###
if ($query->param('strands') =~ /both/i) {
    $parameters .= " -c";
}

################################################################
### return top values
if ($query->param('return') =~ /top/i) {
    $parameters .= " -t";
    $top_scores = $query->param('top_scores');
    if (&IsNatural($query->param('top_scores'))) {
	if ($top_scores == 0) {
	    &FatalError("number of top scores must be >= 1");
	} else {
	    $parameters .= " $top_scores";
	}
    } else {
	&FatalError("Number of top scores must be a strictly positive integer");
    }
}

################################################################
#### case sensitivity
if ($query->param('case') eq "sensitive") {
    $parameters .= ' -CS'; #### [-CS <Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)>]
} elsif ($query->param('case') =~ /mark/) {
    $parameters .= ' -CM'; #### [-CM <Ascii alphabet is case insensitive, but mark the location of lowercase letters>]
}

################################################################
#### unrecognized characters
if ($query->param('unrecognized') eq "errors") {
    $parameters .= " -d0";
} elsif ($query->param('unrecognized') =~ /no warning/) {
    $parameters .= " -d2";
} else {
    $parameters .= " -d1";
}

################################################################
### thresholds ###
if ($query->param('lthreshold_method') =~ /weight/) {
    if (&IsReal($query->param('lthreshold'))) {
	$parameters .= " -ls ".$query->param('lthreshold');

	
    } else {
	&FatalError ("Lower threshold value must ne a real number");
    }
} elsif  ($query->param('lthreshold_method') =~ /p\-value/) {
   ### [-lp <Determine lower-threshold score from a maximum ln(p-value)>]
    if (&IsReal($query->param('lthreshold'))) {
	$parameters .= " -lp ".$query->param('lthreshold');
    } else {
	&FatalError ("Lower threshold value must ne a real number");
    }
} elsif  ($query->param('lthreshold_method') =~ /adjusted information content/) {
   ### [-li <Determine lower-threshold score from adjusted information content>]
    $parameters .= ' -li';
} else {
    &FatalError("Unknown method for estimating lower threshold");
}

if (&IsReal($query->param('uthreshold'))) {
    $parameters .= " -u ".$query->param('uthreshold');
}


#### TEMPORARILY DISACTIVATED, BECAUSE INTERFERES WITH ADJUSTED INFO THRESHOLD
################################################################
#### minimum score for calculating the P-value
# if (&IsReal($query->param('min_calc_P'))) {
#     $parameters .= " -M ".$query->param('min_calc_P');
# } else {
#     &FatalError("Minimum score for calculating P-value must e a real number");
# }

################################################################
### vertically print the matrix
if ($query->param('vertically_print')){
    $parameters .= " -p";
}


################################################################
#### echo the command (for debugging)
print "<pre>$command $parameters</pre>" if ($ECHO);

################################################################
### execute the command ###
if ($query->param('output') eq "display") {
    ### parameters for the piping to the feature map ###
    $feature_file =  "$TMP/$tmp_file_name.ft";
    $features_from_patser_cmd .= " -seq $sequence_file";
    $features_from_patser_cmd .= " -o $feature_file";

    &PipingWarning();

    ### Print the result on Web page
    open RESULT, "$command $parameters & |";
    open FEATURES, "| $features_from_patser_cmd";
    

    print "<PRE>";
    while (<RESULT>) {
	s|$RSA/||g;
	print;
	print FEATURES;
    }
    close FEATURES;
    close RESULT;
    print "</PRE>";

    &PipingForm();

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
<HR SIZE = 3>
<CENTER>
<TABLE>
<TR>
  <TD>
    <H3>Next step</H3>
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
