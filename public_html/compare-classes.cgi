#!/usr/bin/perl
#### redirect error log to a file
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

#### TEMPORARY

$command = "$SCRIPTS/compare-classes";
$tmp_file_name = sprintf "compare-classes.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("compare-classes result");
&ListParameters() if ($ECHO >=2);

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = " -v 1";

#### return a confusion table
if ($query->param('return') eq "matrix") {
    $parameters .= " -matrix"; 
} 

### fields to return
$return_fields = "";

### occurrences
if ($query->param('occ')) {
    $return_fields .= "occ,";

    ### lower threshold on occurrences
    if ($query->param('lth_occ') =~ /^\d+$/) {
	$parameters .= " -lth occ ".$query->param('lth_occ');
    }

    ### upper threshold on occurrences
    if ($query->param('uth_occ') =~ /^\d+$/) {
	$parameters .= " -uth occ ".$query->param('uth_occ');
    }
} 


### percentages
if ($query->param('percent')) {
    $return_fields .= "percent,";

    ### lower threshold on percentages
    if ($query->param('lth_percent') =~ /^\d+$/) {
	$lth_percent = $query->param('lth_percent');
	if ((&IsReal($lth_percent)) &&
	     ($lth_percent <= 100) &&
	     ($lth_percent >= 0)) {
	    $parameters .= " -lth percent ".$lth_percent;
	} else {
	    &FatalError("Lower threshold on percentages: $lth_percent invalid value.");
	}
    }

    ### upper threshold on percentages
    if ($query->param('uth_percent') =~ /^\d+$/) {
	$uth_percent = $query->param('uth_percent');
	if ((&IsReal($uth_percent)) &&
	     ($uth_percent <= 100) &&
	     ($uth_percent >= 0)) {
	    $parameters .= " -uth percent ".$uth_percent;
	} else {
	    &FatalError("Upper threshold on percentages: $uth_percent invalid value.");
	}
    }
} 


### probabilities
if ($query->param('proba')) {
    $return_fields .= "proba,";

    ### lower threshold on probabilities
    if ($query->param('lth_proba') =~ /^\d+$/) {
	$lth_proba = $query->param('lth_proba');
	if ((&IsReal($lth_proba)) &&
	     ($lth_proba <= 100) &&
	     ($lth_proba >= 0)) {
	    $parameters .= " -lth proba ".$lth_proba;
	} else {
	    &FatalError("Lower threshold on probabilities: $lth_proba invalid value.");
	}
    }

    ### upper threshold on probabilities
    if ($query->param('uth_proba') =~ /^\d+$/) {
	$uth_proba = $uth_proba;
	if ((&IsReal($uth_proba)) &&
	     ($uth_proba <= 100) &&
	     ($uth_proba >= 0)) {
	    $parameters .= " -uth proba ".$uth_proba;
	} else {
	    &FatalError("Upper threshold on probabilities: $uth_proba invalid value.");
	}
    }


    ### lower threshold on significance
    if ($query->param('lth_sig') =~ /^\d+$/) {
	$lth_sig = $query->param('lth_sig');
	if (&IsReal($lth_sig)) {
	    $parameters .= " -lth sig ".$lth_sig;
	} else {
	    &FatalError("Lower threshold on significance: $lth_sig invalid value.");
	}
    }

    ### upper threshold on significance
    if ($query->param('uth_sig') =~ /^\d+$/) {
	$uth_sig = $query->param('uth_sig');
	if (&IsReal($uth_sig)) {
	    $parameters .= " -uth sig ".$uth_sig;
	} else {
	    &FatalError("Upper threshold on significance: $uth_sig invalid value.");
	}
    }

} 

### members
if ($query->param('members')) {
    $return_fields .= "members,";
}

$return_fields =~ s/,$//;


unless ($return_fields) {
    &FatalError('You must select at least one return field')
}
$parameters .= " -return $return_fields";

################################################################
## sort keys
if ($query->param('sort_key')) {
    $parameters .= " -sort ".$query->param('sort_key');
}

################################################################
## sort keys
if (&IsNatural($query->param('pop_size'))) {
    $parameters .= " -pop ".$query->param('pop_size');
}

#### load the query classification file
$tmp_query_classes = "${TMP}/${tmp_file_name}_query_classes.tab";
$upload_query_classes = $query->param('upload_query_classes');
if ($upload_query_classes) {
    if ($upload_file =~ /\.gz$/) {
	$tmp_query_classes .= ".gz";
    }
    $type = $query->uploadInfo($upload_query_classes)->{'Content-Type'};
    open CLASS, ">$tmp_query_classes" ||
	&cgiError("Cannot store expected frequency file in temp dir.");
    while (<$upload_query_classes>) {
	print CLASS;
    }
    close CLASS;
} else {
    &FatalError ("Please select the query classification file on your hard drive with the Browse button");
}
$parameters .= " -q $tmp_query_classes";

#### load the reference classification file
$tmp_ref_classes = "${TMP}/${tmp_file_name}_ref_classes.tab";
$upload_ref_classes = $query->param('upload_ref_classes');
if ($upload_ref_classes) {
    if ($upload_file =~ /\.gz$/) {
	$tmp_ref_classes .= ".gz";
    }
    $type = $query->uploadInfo($upload_ref_classes)->{'Content-Type'};
    open CLASS, ">$tmp_ref_classes" ||
	&cgiError("Cannot store expected frequency file in temp dir.");
    while (<$upload_ref_classes>) {
	print CLASS;
    }
    close CLASS;
} else {
    &FatalError ("Please select the reference classification file on your hard drive with the Browse button");
}
$parameters .= " -r $tmp_ref_classes";


print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ECHO >=1);

if ($query->param('output') =~ /display/i) {

#    &PipingWarning();
    
    ### execute the command ###
    $result_file = "$TMP/${tmp_file_name}.res";
    open RESULT, "$command $parameters |";
    
    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);
    
#    &PipingForm();
    print '<HR SIZE=3>';
  
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'), $tmp_file_name);
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $tmp_file_name);
}

print $query->end_html;

exit(0);


