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
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("compare-classes result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

## Commands
$command = $SCRIPTS."/compare-classes -quick";
#$tmp_file_name = sprintf "compare-classes.%s", &AlphaDate();
$prefix = "compare-classes";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

#### read parameters ####
$parameters = " -v 1";

#### return a confusion table
# if ($query->param('return') eq "matrix") {
#     $parameters .= " -matrix"; 
# } 

### fields to return
$return_fields = "";

### occurrences
if ($query->param('occ')) {
    $return_fields .= "occ,";
} 

### frequencies
if ($query->param('freq')) {
    $return_fields .= "freq,";
}

### probabilities
if ($query->param('proba')) {
	$return_fields .= "proba,";
}

### Members
if ($query->param('members')) {
    $return_fields .= "members,";
}

### Jaccard similarity
if ($query->param('jac')) {
    $return_fields .= "jac_sim,";
}

### Entropy
if ($query->param('entropy')) {
    $return_fields .= "entropy";
}
# Thresolds
### lower threshold on occurrences in the query class
if ($query->param('lth_q') =~ /^\d+$/) {
  $parameters .= " -lth Q ".$query->param('lth_q');
}
### upper threshold on occurrences in the query class
if ($query->param('uth_q') =~ /^\d+$/) {
  $parameters .= " -uth Q ".$query->param('lth_q');
}
### lower threshold on occurrences in the reference class
if ($query->param('lth_r') =~ /^\d+$/) {
  $parameters .= " -lth R ".$query->param('lth_r');
}
### upper threshold on occurrences in the reference class
if ($query->param('uth_r') =~ /^\d+$/) {
  $parameters .= " -uth R ".$query->param('lth_r');
}
### lower threshold on occurrences in the intersection
if ($query->param('lth_qr') =~ /^\d+$/) {
  $parameters .= " -lth QR ".$query->param('lth_qr');
}
### upper threshold on occurrences in the intersection
if ($query->param('uth_qr') =~ /^\d+$/) {
  $parameters .= " -uth QR ".$query->param('lth_qr');
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
### lowe threshold on p-value
if ($query->param('lth_pval') =~ /^\d+$/) {
  $lth_pval = $query->param('lth_pval');
  if (&IsReal($lth_pval)) {
    $parameters .= " -lth P_val ".$lth_pval;
  } else {
    &FatalError("Lower threshold on P-value: $lth_pval invalid value.");
  }
}
### upper threshold on p-value
if ($query->param('uth_pval') =~ /^\d+$/) {
  $uth_pval = $query->param('uth_pval');
  if (&IsReal($uth_pval)) {
    $parameters .= " -uth P_val ".$uth_pval;
  } else {
    &FatalError("Upper threshold on P-value: $uth_pval invalid value.");
  }
}
### lowe threshold on e-value
if ($query->param('lth_eval') =~ /^\d+$/) {
  $lth_eval = $query->param('lth_eval');
  if (&IsReal($lth_eval)) {
    $parameters .= " -lth E_val ".$lth_eval;
  } else {
    &FatalError("Lower threshold on e-value: $lth_eval invalid value.");
  }
}
### upper threshold on e-value
if ($query->param('uth_eval') =~ /^\d+$/) {
  $uth_eval = $query->param('uth_eval');
  if (&IsReal($uth_eval)) {
    $parameters .= " -uth E_val ".$uth_eval;
  } else {
    &FatalError("Upper threshold on e-value: $uth_eval invalid value.");
  }
}

### lowe threshold on Jaccard Index
if ($query->param('lth_jac') =~ /^\d+$/) {
  $lth_jac = $query->param('lth_jac');
  if ((&IsReal($lth_jac)) && ($lth_jac < 1)) {
    $parameters .= " -lth jac ".$lth_jac;
  } else {
    &FatalError("Lower threshold on Jaccard Index: $lth_jac invalid value.<br> Must be comprised between 0 and 1");
  }
}
### upper threshold on Jaccard Index
if ($query->param('uth_jac') =~ /^\d+$/) {
  $uth_jac = $query->param('uth_jac');
  if ((&IsReal($lth_jac)) && ($lth_jac < 1)) {
    $parameters .= " -uth jac ".$uth_jac;
  } else {
    &FatalError("Upper threshold on Jaccard Index : $uth_jac invalid value.<br> Must be comprised between 0 and 1");
  }
}
### lower threshold on Mutual information Index
if ($query->param('lth_mi') =~ /^\d+$/) {
  $lth_mi = $query->param('lth_mi');
  if (&IsReal($lth_mi)) {
    $parameters .= " -lth mi ".$lth_mi;
  } else {
    &FatalError("Lower threshold on Mutual information: $lth_mi invalid value.");
  }
}
### upper threshold on Mutual information Index
if ($query->param('uth_mi') =~ /^\d+$/) {
  $uth_mi = $query->param('uth_mi');
  if (&IsReal($uth_mi)) {
    $parameters .= " -uth mi ".$uth_mi;
  } else {
    &FatalError("Upper threshold on Mutual information : $uth_mi invalid value.");
  }
}
$return_fields =~ s/,$//;


unless ($return_fields) {
    &FatalError('You must select at least one return field')
}
$parameters .= " -return $return_fields";

################################################################
## sort keys
if ($query->param('sort_key')) {
    $sort_key = $query->param("sort_key");
    if ($sort_key eq "Mutual information") {
       $sort_key = "'I(Q,R)'";
    } elsif ($sort_key eq "Jaccard index") {
       $sort_key = 'jac_sim';
    }
    $parameters .= " -sort ".$sort_key;
}

################################################################
## Population size
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


&ReportWebCommand($command." ".$parameters);

if ($query->param('output') =~ /display/i) {

#    &PipingWarning();
    
    ### execute the command ###
    $result_file = &RSAT::util::get_pub_temp()."/${tmp_file_name}.res";
    open RESULT, "$command $parameters |";
    
    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);
    
#    &PipingForm();
    print '<HR SIZE=3>';

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $tmp_file_name);
}

print $query->end_html;

exit(0);


