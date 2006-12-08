#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis.cgi,v 1.22 2006/12/08 15:57:57 jvanheld Exp $
#
# Time-stamp: <2003-10-11 00:30:17 jvanheld>
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
require "RSA.disco.lib";
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
$parameters = " -v 1 -sort -timeout 3600 ";
#$parameters = " -v 1 -sort -return proba,rank -timeout 3600 ";


### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
$purge = $query->param('purge');
if ($purge) {
    #### purge sequence option
#    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format | $dyad_analysis_command ";
    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format -o $sequence_file.purged; $dyad_analysis_command -i $sequence_file.purged ";
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
#if ($query->param('lth_occ_sig') =~ /^-{0,1}[\d\.\-+e]+$/i) {
#  $parameters .= "  -thosig ".$query->param('lth_occ_sig');
#}

&CGI_return_fields();

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
} elsif ($query->param('freq_estimate') =~ /upload/i) {
    $exp_freq_file = "${TMP}/$tmp_file_name.expfreq";
    $upload_freq_file = $query->param('upload_freq_file');
    if ($upload_freq_file) {
	if ($upload_file =~ /\.gz$/) {
	    $exp_freq_file .= ".gz";
	}
	$type = $query->uploadInfo($upload_freq_file)->{'Content-Type'};
	open FREQ, ">$exp_freq_file" ||
	    &cgiError("Cannot store expected frequency file in temp dir.");
	while (<$upload_freq_file>) {
	    print FREQ;
	}
	close FREQ;
	$parameters .= " -expfreq $exp_freq_file";
    } else {
	&FatalError ("If you want to upload an expected frequency file, you should specify the location of this file on your hard drive with the Browse button");
    }

} else {
    unless ($query->param('freq_estimate') eq 'monads') {
	&FatalError("Invalid expected frequency calibration");
    }
}

$command .= $parameters;

print "<PRE><B>Command:</B> $command </PRE>" if ($ECHO);

&SaveCommand($command, "$TMP/$tmp_file_name");

if ($query->param('output') eq "display") {  

    &PipingWarning();

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$command | ";


    ### Print result on the web page
    print '<H4>Result</H4>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

    #### pattern assembly ####
    if ((&IsReal($query->param('lth_occ_sig'))) && ($query->param('lth_occ_sig')>= -1)) {

# 	$pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 0";
# 	if ($query->param('strand') =~ /single/) {
# 	    $pattern_assembly_command .= " -1str";
# 	} else {
# 	    $pattern_assembly_command .= " -2str";
# 	}
# 	$pattern_assembly_command .= " -maxfl 1 -subst 0 ";
# 	print "<PRE>pattern-assembly command: $pattern_assembly_command<P>\n</PRE>" if ($ECHO >=1);
	
# 	print "<H4>Pattern assembly</H4>\n";
# 	open CLUSTERS, "$pattern_assembly_command -i $result_file |";
# 	print "<PRE>\n";
# 	while (<CLUSTERS>) {
# 	    s|$RSA/||g;
# 	    print;
# 	}
# 	print "</PRE>\n";
# 	close(CLUSTERS);


	## Assemble the significant patterns with pattern-assembly
	$assembly_file = "$TMP/$tmp_file_name.asmb";
	$pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 1";
	if ($query->param('strand') =~ /single/) {
	  $pattern_assembly_command .= " -1str";
	} else {
	  $pattern_assembly_command .= " -2str";
	}
	$pattern_assembly_command .= "  -i $result_file";
	$pattern_assembly_command .= "  -o $assembly_file";

	print "<H2>Pattern assembly</H2>\n";
	print "<PRE>pattern-assembly command: $pattern_assembly_command<P>\n</PRE>" if ($ECHO >=1);
	system "$pattern_assembly_command";
	open ASSEMBLY, $assembly_file;
	print "<PRE>\n";
	while (<ASSEMBLY>) {
	  s|$RSA/||g;
	  print;
	}
	print "</PRE>\n";
	close(ASSEMBLY);


	## Convert pattern-assembly result into matrix profiles to be displayed on the screen
	$profile_file = "$TMP/$tmp_file_name.profile";
	$profile_command = "$SCRIPTS/convert-matrix -v 1 ";
	$profile_command .= " -in_format assembly -out_format patser";
	$profile_command .= " -return profile,parameters";
	$profile_command .= " -i $assembly_file";
	$profile_command .= " -o $profile_file";
	print "<PRE>command to generate profiles: $profile_command<P>\n</PRE>" if ($ECHO >=1);
	system "$profile_command";
	print "<H2>Position-specific scoring matrices (PSSM)</H2>\n";
	open PROFILE, $profile_file;
	print "<PRE>\n";
	while (<PROFILE>) {
	  s|$RSA/||g;
	  print;
	}
	print "</PRE>\n";
	close(PROFILE);

	## Convert pattern-assembly result into PSSM for piping to other tools
	$pssm_file = "$TMP/$tmp_file_name.pssm";
	$pssm_command = "$SCRIPTS/convert-matrix -v 0 ";
	$pssm_command .= " -in_format assembly -out_format patser";
	$pssm_command .= " -return counts";
	$pssm_command .= " -i $assembly_file";
	$pssm_command .= " -o $pssm_file";
	print "<PRE>command to generate matrices: $pssm_command<P>\n</PRE>" if ($ECHO >=1);
	system "$pssm_command";

    }

    &PipingForm();
    
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command", $query->param('user_email'));
} else {
    &EmailTheResult("$command", $query->param('user_email'));
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

<TD VALIGN=BOTTOM ALIGN=CENTER>
<H3>Next step</H3>
</TD>

<TD VALIGN=BOTTOM ALIGN=CENTER>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="string-based pattern matching (dna-pattern)">
</FORM>
</TD>

<td valign=bottom align=center>
<b><font color="red">New !</font></b>
<FORM METHOD="POST" ACTION="matrix-scan_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$pssm_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="tab">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="matrix-based pattern matching (matrix-scan)">
</FORM>
</TD>

<TD VALIGN=BOTTOM ALIGN=CENTER>
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

