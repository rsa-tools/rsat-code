#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis.cgi,v 1.10 2001/09/21 06:30:03 jvanheld Exp $
#
# Time-stamp: <2001-09-21 08:27:10 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$output_context = "cgi";

$dyad_analysis_command = "$SCRIPTS/dyad-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "dyad-analysis.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("dyad-analysis result");
&ListParameters if $ECHO;

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "-v -sort -return proba,rank -timeout 3600 ";

#### purge sequence option
$purge = $query->param('purge');

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
if ($purge) {
    $command= "$convert_seq_command -i $sequence_file -from $sequence_format -to fasta | $purge_sequence_command | $dyad_analysis_command -format fasta ";
} else {
    $command= "$dyad_analysis_command -i $sequence_file -from $sequence_format ";
}

### sequence file
#($sequence_file,$sequence_format) = &GetSequenceFile();
#$parameters .= " -i $sequence_file -format $sequence_format";

### dyad type
if ($query->param('dyad_type') =~ /invert/) {
  $parameters .= " -type ir";
} elsif ($query->param('dyad_type') =~ /direct/) {
  $parameters .= " -type dr";
} elsif ($query->param('dyad_type') =~ /repeat/) {
  $parameters .= " -type rep";
} else {
  $parameters .= " -type any";
}

### single strand search
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}

### overlapping matches
if ($query->param('noov')) {
  $parameters .= " -noov";
} 

### organism
if ($organism = $query->param('organism')) {
  $parameters .= " -org $organism";
} 


### threshold ###
if ($query->param('occ_significance_threshold') =~ /^-{0,1}[\d\.\-+e]+$/i) {
  $parameters .= "  -thosig ".$query->param('occ_significance_threshold');
}

#### oligo size ####
if (&IsNatural($query->param('oligo_size'))) {
  $oligo_length = $query->param('oligo_size') ;
} 
$parameters .= " -l $oligo_length";

#### spacing ####
unless ((&IsNatural($query->param('spacing_from'))) && (&IsNatural($query->param('spacing_from'))))  {
  &cgiError("Invalid spacing specification.");    
}
if ($query->param('spacing_from') >$query->param('spacing_to')) {
  $spacing = $query->param('spacing_to')."-".$query->param('spacing_from') ;   
} else {
  $spacing = $query->param('spacing_from')."-".$query->param('spacing_to') ;   
}
$parameters .= " -spacing $spacing";

#### expected frequency estimation ####
if ($query->param('exp_freq') =~ /non\-coding/i) {
  $parameters .= " -ncf";
} 

print "<PRE><B>Command:</B> $command $parameters </PRE>" if $ECHO;

if ($query->param('output') eq "display") {  

    &PipingWarning();

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$command $parameters | ";


    ### Print result on the web page
    print '<H4>Result</H4>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

    #### pattern assembly ####
    if ((&IsReal($query->param('occ_significance_threshold'))) && ($query->param('occ_significance_threshold')>= -1)) {
	$fragment_assembly_command = "$SCRIPTS/pattern-assembly -v";
	if ($query->param('strand') =~ /single/) {
	    $fragment_assembly_command .= " -1str";
	} else {
	    $fragment_assembly_command .= " -2str";
	}
	$fragment_assembly_command .= "-maxfl 2 -subst 0 ";
	
	print "<H3>Pattern assembly</H3>\n";
	open CLUSTERS, "$fragment_assembly_command -i $result_file |";
	print "<PRE>\n";
	while (<CLUSTERS>) {
	    print;
	}
	print "</PRE>\n";
	close(CLUSTERS);
    }

    &PipingForm();
    
} else {
    #### send e-mail with the result
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your e-mail address: <P>";
	print "$address</B><P>";
	system "$dyad_analysis_command $parameters | $mail_command $address &"; 
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

unless ($graph_request) {
    print "<HR SIZE = 3>";
    print $query->end_html;
}


exit(0);

sub PipingForm {
    ### prepare data for piping
    $title = $query->param('title');
    $title =~ s/\"/\'/g;
    print <<End_of_form;
<TABLE>
<TR>
<TD>
<H3>Next step</H3>
</TD>
<TD>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="pattern matching (dna-pattern)">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form
}
