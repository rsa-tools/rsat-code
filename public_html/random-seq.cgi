#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$random_seq_command = "$SCRIPTS/random-seq";

$size_limit = 1000000;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Random sequence result");

#### update log file ####
&UpdateLogFile;

&ListParameters() if ($ECHO);

#### read parameters ####
$parameters = "";
$length = $query->param('length');
if (&IsNatural($length)) {
    $parameters .= " -l $length ";
} else {
    &FatalError("Sequence length must be a natural number");
}
$repet = $query->param('repet');
if (&IsNatural($repet)) {
    $parameters .= " -r $repet";
} else {
    &FatalError("Repetitions must be a natural number");
}
if ($repet * $length == 0) {
    &FatalError("Sequence length and repeitions must be non-null");
}
if ($repet*$length > $size_limit) {
    &FatalError("The web interface does not support queries of this size. Maximum size per query (length * repetitions) = ".$size_limit);
}

if (&IsNatural($query->param('lw'))) {
    $parameters .= " -lw ".$query->param('lw');
}
&CheckInputSeqFormat($query->param('format'));
$parameters .= " -format ".$query->param('format');

if ($query->param('proba') eq "alphabet") {
    $at_freq = $query->param('ATfreq');
    $cg_freq = $query->param('CGfreq');
    unless ((&IsReal($at_freq))  && (&IsReal($cg_freq))) {
	&FatalError("Invalid residue frequencies");
    }
    $parameters .= " -a a:t $at_freq c:g $cg_freq ";
} elsif ($query->param('proba') eq "expfreq") {
    $oligo_length = $query->param('oligo_size');
    $freq_file = "${oligo_length}nt";
    
    if ($query->param('seq_type') eq "genomic") {
	$freq_file .= ".genomic.freq";
    } elsif ($query->param('seq_type') eq "coding") {
	$freq_file .= ".coding.freq";
    } elsif ($query->param('seq_type') eq "non coding") {
	$freq_file .= ".non-coding.freq";
    }
    
    $parameters .= " -expfreq $RSA/data/yeast/oligo-frequencies/$freq_file ";
} 

print "<PRE>command: $random_seq_command $parameters<P>\n</PRE>" if $ECHO;

### execute the command ###
if ($query->param('output') eq "display") {
    ### Print the result on Web page
    open RESULT, "$random_seq_command $parameters  & |";
    print "<PRE>";
    while (<RESULT>) {
	print "$_";
    }
    print "</PRE>";
    close RESULT;

} else {
    ### send an e-mail with the result ###
    if ($query->param('user_email') =~ /(\S+\@\S+)/) {
	$address = $1;
	print "<B>Result will be sent to your account: <P>";
	print "$address</B><P>";
	system "$random_seq_command $parameters | $mail_command $address &";
    } else {
	if ($query->param('user_email') eq "") {
	    print "<B>ERROR: you did not enter your e-mail address<P>";
	} else {
	    print "<B>ERROR: the e-mail address you entered is not valid<P>";
	    print $query->param('user_email')."</B><P>";      
	}
    } 
}
print "<HR SIZE = 3>";
print $query->end_html;

exit(0);

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

