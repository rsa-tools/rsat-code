#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
require "$RSA/public_html/genome-scale.lib.pl";

$patser_command = "$BIN/patser";
$matrix_from_transfac_command = "$SCRIPTS/matrix-from-transfac";
$matrix_from_gibbs_command = "$SCRIPTS/matrix-from-gibbs";
$convert_seq_command = "$SCRIPTS/convert-seq";
#$features_from_patser_cmd = "$SCRIPTS/features-from-patser";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
#$add_yeast_link_command = "$SCRIPTS/add-yeast-link";
$add_orf_function_command = "$SCRIPTS/add-orf-function";
$link_command = "$SCRIPTS/add-yeast-link";
$tmp_file_name = sprintf "genome-scale-patser.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the header of the result page
&RSA_header("patser result ".$query->param("title"));

#&ListParameters;


#### update log file ####
&UpdateLogFile;

&ReadRetrieveSeqParams();

################################################################
#
# patser parameters
#

#### read parameters ####
$patser_parameters = " -A a:t c:g ";

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
$patser_parameters .= " -m $matrix_file";

### sequence file ####
#($sequence_file,$sequence_format) = &GetSequenceFile("wconsensus", 1);
#$patser_parameters .= " -f $sequence_file";

### strands ###
if ($query->param('strands') =~ /both/i) {
    $patser_parameters .= " -c";
}

### top value only ###
if ($query->param('return') =~ /top/i) {
    $patser_parameters .= " -t";
}

### thresholds ###
if (&IsReal($query->param('lthreshold'))) {
    $patser_parameters .= " -ls ".$query->param('lthreshold');
    $patser_parameters .= " -M ".$query->param('lthreshold');
} 
if (&IsReal($query->param('uthreshold'))) {
    $patser_parameters .= " -u ".$query->param('uthreshold');
}


#### pseudo-counts
if (&IsReal($query->param('pseudo_counts'))) {
    $parameters .= " -b ".$query->param('pseudo_counts');
}


### parameters for the piping to the feature map ###
$feature_file =  "$TMP/$tmp_file_name.ft";
#$features_from_patser_cmd .= " -seq $sequence_file";
#$features_from_patser_cmd .= " -o $feature_file";

$command = "$retrieve_seq_command $retrieve_seq_parameters ";
$command .= "| $patser_command $patser_parameters ";

if ($org eq "yeast") { #### not yet supported for other organisms
    $command .= "| $add_orf_function_command  ";
    $command .= "| $link_command  ";
}

### execute the command ###
if ($query->param("output") =~ /display/i) {

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
#    print "<PRE>$command</PRE>";
#    print "<PRE>$features_from_patser_cmd</PRE>";
    open RESULT, "$command & |";
#    open FEATURES, "| $features_from_patser_cmd";

#    &PipingWarning();

    ### Print the result on Web page
    print "<PRE>";
    while (<RESULT>) {
	s|$RSA/||g;
	print;
#	print FEATURES;
    }
#    close FEATURES;
    close RESULT;
    print "</PRE>";
    
#    &PipingForm();

#    ### Print the result on Web page
#    &PrintHtmlTable(RESULT, $result_file);
#    close RESULT;
    
    
} else {
    ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your e-mail address: <P>";
	print "$address</B><P>";
	system "$command | $mail_command $address &";
    } else {
	if ($query->param('user_email') eq "") {
	    print "<B>ERROR: you did not enter your e-mail address<P>";
	} else {
	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
	    print "$query->param('user_email')</B><P>";      
	}
    } 
}

print "<HR SIZE = 3>";

print $query->end_html;


exit(0);

sub PipingForm {
  ### prepare data for piping
  $title = $query->param("title");
  $title =~ s/\"/'/g;
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
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="feature_file" VALUE="$result_file">
<INPUT type="hidden" NAME="format" VALUE="patser">
<INPUT type="hidden" NAME="fill_form" VALUE="on">
<INPUT type="submit" VALUE="feature map">
</FORM>
</TD>
</TR>
</TABLE>
</CENTER>
End_of_form
}
