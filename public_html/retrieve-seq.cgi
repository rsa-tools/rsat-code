#!/usr/bin/perl
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

$command = "$SCRIPTS/retrieve-seq";
$tmp_file_name = sprintf "retrieve-seq.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("retrieve-seq result");


#### update log file ####
&UpdateLogFile;

&ListParameters() if ($ECHO >= 2);

#### organism
if (defined($supported_organism{$query->param('organism')})) {
    $organism = $supported_organism{$query->param('organism')}->{'name'};
    $parameters .= " -org ".$query->param('organism');
} else {
    &cgiError("Organism '",
	      $query->param('organism'),
	      "' is not supported on this web site.");
}

### feature type
if ($query->param('feattype')) {
    my ($feattype) = split " ", $query->param('feattype'); ### take the first word
    $parameters .= " -feattype ".$feattype;
}

### sequence type
if ($query->param('sequence_type')) {
    ($seq_type) = split " ", $query->param('sequence_type'); ### take the first word
    $parameters .= " -type ".$seq_type;
}

### output format ###
$out_format = lc($query->param('format'));
if ($accepted_output_seq{$out_format}) {
    $parameters .= " -format $out_format";
}

### sequence label
$seq_label = lc($query->param('seq_label'));
if (($seq_label =~ /gene/) && 
    ($seq_label =~ /orf/)) {
    $parameters .= " -label orf_gene";
} elsif ($seq_label =~ /gene/) {
    $parameters .= " -label gene";
} elsif ($seq_label =~ /orf/) {
    $parameters .= " -label orf";
} elsif ($seq_label =~ /full/) {
    $parameters .= " -label full";
} else {
    &cgiError("Invalid option for sequence label '$seq_label'");
}

### limits ###
if (&IsInteger($query->param('from'))) {
    $parameters .= " -from ".$query->param('from');
}  
if (&IsInteger($query->param('to'))) {
    $parameters .= " -to ".$query->param('to');
}


### prevent orf overlap ###
if (lc($query->param('noorf')) eq "on") {
    $noorf = 1;
    $parameters .= " -noorf ";
}

#### queries ####
if ($query->param('genes') eq "all") {
    ### take all genes as query
    $parameters .= " -all ";
} elsif ($query->param('uploaded_file')) {
    $upload_file = $query->param('uploaded_file');
    $gene_list_file = "${TMP}/${tmp_file_name}.genes";
    if ($upload_file =~ /\.gz$/) {
	$gene_list_file .= ".gz";
    }
    $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open SEQ, ">$gene_list_file" ||
	&cgiError("Cannot store gene list file in temporary directory");
    while (<$upload_file>) {
	print SEQ;
    }
    close SEQ;
    $parameters .= " -i $gene_list_file ";

} else {
    $gene_selection = $query->param('gene_selection');
    unless ($gene_selection =~ /\S/) {
	&cgiError("You should enter at least one gene identifier in the query box. <P>Read the on-line manual for more information.");
    }
    $gene_selection .= "\n";			### make sure there is a carriage return at the end
    
    
    @query_lines = split("\n",$gene_selection);
    foreach $line (@query_lines) {
	next if ($line =~ /^;/);
	next if ($line =~ /^--/);
	next if ($line =~ /^\#/);
	if ($line =~ /(\S+)/) {
	    $parameters .= " -q '$1' "; 
	}
    }
}

print  "<PRE><B>Command :</B> $command $parameters</PRE><P>" if ($ECHO >= 1);

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
    open RESULT, "$command $parameters |";

    ### print the result ### 
    &PipingWarning();
    if ($query->param('output') =~ /server/i) {
	&Info("The result will appear below ...");
    }

    print '<H3>Result</H3>';

    ### open the mirror file ###
    $mirror_file = "$TMP/$tmp_file_name.res";
    if (open MIRROR, ">$mirror_file") {
	$mirror = 1;
	&DelayedRemoval($mirror_file);
    }
    
    print "<PRE>";
    while (<RESULT>) {
	print "$_" unless ($query->param('output') =~ /server/i);
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close RESULT;
    close MIRROR if ($mirror);
    
    if ($query->param('output') =~ /server/i) {
	$result_URL = "${WWW_RSA}/tmp/${tmp_file_name}.res";
	print ("The result is available at the following URL: ", "\n<br>",
	       "<a href=${result_URL}>${result_URL}</a>",
	       "<p>\n");
    }

    ### prepare data for piping
    &PipingForm();
    
    print "<HR SIZE = 3>";

#} elsif 
#    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

print $query->end_html;

exit(0);

################################################################
#
# Pipe the result to another command
#
sub PipingForm {
    #### choose background model for oligo-analysis
    if ($seq_type =~ /upstream/) {
	if ($noorf) {
	    $background = "upstream-noorf";
	} else {
	    $background = "upstream";
	}
    } else {
	$background = "intergenic";
    }

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
	<INPUT type="hidden" NAME="background" VALUE="$background">
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
	<INPUT type="hidden" NAME="background" VALUE="$background">
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
	<INPUT type="hidden" NAME="background" VALUE="$background">
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
	<INPUT type="submit" value="patser (matrices)">
	</FORM>
    </TD>


    <TD>
    &nbsp;
    </TD>

    <TD>
    &nbsp;
    </TD>

</TR>


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


</TABLE>
End_of_form
}
