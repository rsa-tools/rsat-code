#!/usr/bin/perl
############################################################
#
# $Id: patser.cgi,v 1.19 2003/06/03 22:18:33 jvanheld Exp $
#
# Time-stamp: <2003-06-04 00:18:03 jvanheld>
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
$features_from_patser_cmd = "$SCRIPTS/features-from-patser -v 1";
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

#### alphabet
my $alphabet = $query->param('alphabet') || " a:t c:g ";
$patser_parameters = " -A $alphabet";


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
$patser_parameters .= " -m $matrix_file";

#### [-w <Matrix is a weight matrix>]
if ($query->param('matrix_is_weight')) {
    $patser_parameters .= " -w";
} else {
    #### pseudo-counts and weights are mutually exclusive
    if (&IsReal($query->param('pseudo_counts'))) {
	$patser_parameters .= " -b ".$query->param('pseudo_counts');
    }

}

#### [-v <Vertical matrix---rows correspond to positions>]
if ($query->param('matrix_is_vertical')) {
    $patser_parameters .= " -v";
}


################################################################
#### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1);
$patser_parameters .= " -f $sequence_file";

################################################################
#### strands 
if ($query->param('strands') =~ /both/i) {
    $patser_parameters .= " -c";
}

################################################################
### return top values
if ($query->param('return') =~ /top/i) {
    $patser_parameters .= " -t";
    $top_scores = $query->param('top_scores');
    if (&IsNatural($query->param('top_scores'))) {
	if ($top_scores == 0) {
	    &FatalError("number of top scores must be >= 1");
	} else {
	    $patser_parameters .= " $top_scores";
	}
    } else {
	&FatalError("Number of top scores must be a strictly positive integer");
    }
}

################################################################
#### case sensitivity
if ($query->param('case') eq "sensitive") {
    $patser_parameters .= ' -CS'; #### [-CS <Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)>]
} elsif ($query->param('case') =~ /mark/) {
    $patser_parameters .= ' -CM'; #### [-CM <Ascii alphabet is case insensitive, but mark the location of lowercase letters>]
}

################################################################
#### unrecognized characters
if ($query->param('unrecognized') eq "errors") {
    $patser_parameters .= " -d0";
} elsif ($query->param('unrecognized') =~ /no warning/) {
    $patser_parameters .= " -d2";
} else {
    $patser_parameters .= " -d1";
}

################################################################
### thresholds ###

#### lower threshold on the weight
if ($query->param('lthreshold_method') =~ /weight/) {
    if (&IsReal($query->param('lthreshold'))) {
	$patser_parameters .= " -ls ".$query->param('lthreshold');
    } elsif ($query->param('lthreshold') eq 'none') {
	### no lower threshold
    } else {
	&Warning("Lower threshold ignored (not a real value)");
    }
    
#### lower threshold on P-value
} elsif  ($query->param('lthreshold_method') =~ /p\-value/) {
   ### [-lp <Determine lower-threshold score from a maximum ln(p-value)>]
    if (&IsReal($query->param('lthreshold'))) {
	$patser_parameters .= " -lp ".$query->param('lthreshold');
    } elsif ($query->param('lthreshold') eq 'none') {
	### no lower threshold
    } else {
	&FatalError ("Lower threshold value must be a real number");
    }

#### automatic threshold on the basis of adjusted information content
} elsif  ($query->param('lthreshold_method') =~ /adjusted information content/) {
    ### [-li <Determine lower-threshold score from adjusted information content>]
    $patser_parameters .= ' -li';
} else {
    &FatalError("Unknown method for estimating lower threshold");
}

#### upper threshold
if (&IsReal($query->param('uthreshold'))) {
    $patser_parameters .= " -u ".$query->param('uthreshold');
}


#### TEMPORARILY DISACTIVATED, BECAUSE INTERFERES WITH ADJUSTED INFO THRESHOLD
################################################################
#### minimum score for calculating the P-value
# if (&IsReal($query->param('min_calc_P'))) {
#     $patser_parameters .= " -M ".$query->param('min_calc_P');
# } else {
#     &FatalError("Minimum score for calculating P-value must e a real number");
# }

################################################################
### vertically print the matrix
if ($query->param('vertically_print')){
    $patser_parameters .= " -p";
}

$command .= " $patser_parameters";

################################################################
### parameters for the piping to the feature map ###
#$feature_file =  "$TMP/$tmp_file_name.ft";
$features_from_patser_cmd .= " -seq $sequence_file ";
#### origin 
if ($query->param('origin') =~ /end/i) {
    $features_from_patser_cmd .= " -origin -0";
}

#### flanking residues for the matching sequences
if ($query->param('flanking') =~ /^\d+$/) {
    $features_from_patser_cmd .= " -N ".$query->param('flanking');
}
$features_from_patser_cmd .= " -origin -0";
#$features_from_patser_cmd .= " -o $feature_file";

### return matching positions
if ($query->param('positions')) {
    $features_from_patser_cmd .= " -return matches";
}

### return score table
if ($query->param('table')) {
    $features_from_patser_cmd .= " -return table";
}


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

    unless ($query->param('table')) {
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
