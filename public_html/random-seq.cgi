#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
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
$command = "$SCRIPTS/random-seq";

$size_limit = 5e+6;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Random sequence result");

#### update log file ####
&UpdateLogFile;

&ListParameters() if ($ECHO >= 2);

################################################################
#### read parameters ####
$parameters = "";

#### sequence length
$length = $query->param('length');
if (&IsNatural($length)) {
    $parameters .= " -l $length ";
} else {
    &FatalError("Sequence length must be a natural number");
}

#### number of repetitions
$repet = $query->param('repet');
if (&IsNatural($repet)) {
    $parameters .= " -r $repet";
} else {
    &FatalError("Repetitions must be a natural number");
}

#### check lengths and repetitions
if ($repet * $length == 0) {
    &FatalError("Sequence length and repeitions must be non-null");
}
if ($repet*$length > $size_limit) {
    &FatalError("The web interface does not support queries of this size. Maximum size per query (length * repetitions) = ".$size_limit);
}

#### line width
if (&IsNatural($query->param('lw'))) {
    $parameters .= " -lw ".$query->param('lw');
}

#### output format
$out_format = $query->param('format');
&CheckOutputSeqFormat($out_format);
$parameters .= " -format $out_format ";

#### alphabet
if ($query->param('proba') eq "alphabet") {
    $at_freq = $query->param('ATfreq');
    $cg_freq = $query->param('CGfreq');
    unless ((&IsReal($at_freq))  && (&IsReal($cg_freq))) {
	&FatalError("Invalid residue frequencies");
    }
    $parameters .= " -a a:t $at_freq c:g $cg_freq ";

#### expected frequency estimation ####
} elsif ($query->param('proba') =~ /ncf/i) {
    ### check organism
    unless ($organism = $query->param('organism')) {
	&cgiError("You should specify an organism to use intergenic frequency calibration");
    }
    unless (defined(%{$supported_organism{$organism}})) {
	&cgiError("Organism $organism is not supported on this site");
    }
    $oligo_size = $query->param("oligo_size");
    unless (&IsNatural($oligo_size)) {
	&cgiError("Invalid oligonucleotide length $oligo_size");
    }
    $parameters .= " -ncf -org $organism -ol $oligo_size";
}

print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ECHO >= 1);


### execute the command ###
if ($query->param('output') eq "display") {
    open RESULT, "$command $parameters |";
    
    ### open the mirror file ###
    $mirror_file = "$TMP/$tmp_file_name.res";
    if (open MIRROR, ">$mirror_file") {
	$mirror = 1;
	&DelayedRemoval($mirror_file);
    }
    
    
    ### print the result ### 
    &PipingWarning();
    print '<H3>Result</H3>';
    print "<PRE>";
    while (<RESULT>) {
	print "$_";
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close RESULT;
    
    ### prepare data for piping
    &PipingForm();

#      &PipingWarning();
#      ### Print the result on Web page
#      open RESULT, "$command $parameters  & |";
#      print "<PRE>";
#      while (<RESULT>) {
#  	print "$_";
#      }
#      print "</PRE>";
#      close RESULT;
#      ### prepare data for piping
#      &PipingForm();
    print "<HR SIZE = 3>";

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}
print $query->end_html;

exit(0);


################################################################
#
# Pipe the result to other commands
#
sub PipingForm {
  print <<End_of_form;
<HR SIZE = 3>
<TABLE CELLSPACING=0 CELLPADDING=10 BORDER=0 NOWRAP BGCOLOR= #FFEEDD>

<TR VALIGN="top" ALIGN="center">
    <TD COLSPAN=5 BGCOLOR=	#FFEEDD>
	<H3>Next step</H3>
    </TD>

</TR>

<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=		#FFEEDD>
	<B>Pattern discovery</B><BR>
	(unknown patterns)
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="oligo-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="oligonucleotide analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dyad-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dyad analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="position-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="position analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="consensus_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="consensus">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="gibbs_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="gibbs sampler">
	</FORM>
    </TD>

</TR>

<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=#FFEEDD>
	<B>Pattern matching</B><BR>
	(known patterns)
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dna-pattern (IUPAC)">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="patser_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism">
	<INPUT type="hidden" NAME="from" VALUE="$from">
	<INPUT type="hidden" NAME="to" VALUE="$to">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="patser (matrices)";
	</FORM>
    </TD>

    <TD>
    &nbsp;
    </TD>

    <TD>
    &nbsp;
    </TD>

</TR>


<!--
<TR VALIGN="top" ALIGN="center">

    <TD BGCOLOR=		#FFEEDD>
	<B>Utilities</B>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="purge-sequence_form.cgi">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="purge sequence">
	</FORM>
    </TD>

    <TD>
    &nbsp;
    </TD>

    <TD>
    &nbsp;
    </TD>

    <TD>
    &nbsp;
    </TD>

</TR>
-->

</TABLE>
End_of_form
}

