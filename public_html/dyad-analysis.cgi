#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis.cgi,v 1.17 2002/11/16 14:11:37 jvanheld Exp $
#
# Time-stamp: <2002-11-16 07:46:24 jvanheld>
#
############################################################
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$dyad_analysis_command = "$SCRIPTS/dyad-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "dyad-analysis.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("dyad-analysis result");
&ListParameters if ($ECHO >=2);

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = " -v 1 -sort -return proba,rank -timeout 3600 ";

#### purge sequence option
$purge = $query->param('purge');

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
if ($purge) {
    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format | $dyad_analysis_command ";
} else {
    $command= "$dyad_analysis_command -i $sequence_file  ";
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
if ($query->param('freq_estimate') eq 'background') {
    $parameters .= " -bg ".$query->param('background');
} else {
    unless ($query->param('freq_estimate') eq 'monads') {
	&FatalError("Invalid expected frequency calibration");
    }
}

print "<PRE><B>Command:</B> $command $parameters </PRE>" if ($ECHO);

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
	$pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 0";
	if ($query->param('strand') =~ /single/) {
	    $pattern_assembly_command .= " -1str";
	} else {
	    $pattern_assembly_command .= " -2str";
	}
	$pattern_assembly_command .= " -maxfl 2 -subst 0 ";
	
	print "<H4>Pattern assembly</H4>\n";
	open CLUSTERS, "$pattern_assembly_command -i $result_file |";
	print "<PRE>\n";
	while (<CLUSTERS>) {
	    s|$RSA/||g;
	    print;
	}
	print "</PRE>\n";
	close(CLUSTERS);
    }

    &PipingForm();
    
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}



print $query->end_html;


exit(0);

sub PipingForm {
    ### prepare data for piping

    #### title
    $title = $query->param('title');
    $title =~ s/\"/\'/g;

    #### strand for pattern-assembly
    if ($query->param('strand') =~ /single/) {
	$strand_opt .= " sensitive";
    } else {
	$strand_opt .= " insensitive";
    }
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



</TD>
<TD>
<FORM METHOD="POST" ACTION="pattern-assembly_form.cgi">
<INPUT type="hidden" NAME="local_pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="subst" VALUE=0>
<INPUT type="hidden" NAME="maxfl" VALUE=1>
<INPUT type="hidden" NAME="sc" VALUE="auto">
<INPUT type="hidden" NAME="strand" VALUE=$strand_opt>
<INPUT type="submit" value="pattern assembly">
</FORM>
</TD>


</TR>
</TABLE>
End_of_form
}
