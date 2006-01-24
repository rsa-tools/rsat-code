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

$tmp_file_name = sprintf "retrieve-seq.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;


#$ECHO=2;


### print the header
&RSA_header("retrieve-seq result");


#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ECHO >= 2);

$parameters = "";

## ##############################################################
## The query contains only identifiers (no need to load synonyms)
if ($query->param('ids_only')) {
    $parameters .= " -ids_only";
}

################################################################
## Single or multi-genome query
if ($query->param('single_multi_org') eq 'multi') {
    $command = "$SCRIPTS/retrieve-seq-multigenome";

    &cgiMessage(join("<P>",
		     "The computation can take a more or less important time depending on the taxon size.",
		     "If the answer does not appear in due time, use the option <i>output email</i>"));


#     my $gene_col = $query->param('gene_col');
#     if (&IsNatural($gene_col) && ($gene_col > 0)) {
# 	$parameters .= " -gene_col ".$gene_col;
#     }
#     my $org_col = $query->param('org_col');
#     if (&IsNatural($org_col) && ($org_col > 0)) {
# 	$parameters .= " -org_col ".$org_col;
#     }
    
} else {
    $command = "$SCRIPTS/retrieve-seq";

}



#### organism
$organism = $query->param('organism');
if (defined($supported_organism{$organism})) {
    $organism_name = $supported_organism{$organism}->{'name'};

    $parameters .= " -org ".$organism unless ($query->param('single_multi_org') eq 'multi'); ## For multi-genome retrieval, the query organism name is passed to pattern discovery programs, but it is not necessary
} else {
    &cgiError("Organism '",
	      $organism,
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
if ($seq_label eq 'gene identifier') {
    $parameters .= " -label id";
} elsif ($seq_label eq 'gene name') {
    $parameters .= " -label name";
} elsif ($seq_label eq 'gene identifier + name') {
    $parameters .= " -label id,name";
} elsif ($seq_label eq 'gene identifier + organism + gene name') {
    $parameters .= " -label id,organism_name,name";
} elsif ($seq_label eq 'full identifier') {
    $parameters .= " -label id,name,organism_name,sequence_type,current_from,current_to,ctg,orf_strand,reg_left,reg_right";
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
    my $gene_selection = $query->param('gene_selection');
    $gene_selection =~ s/\r/\n/g;
    my @gene_selection = split ("\n", $gene_selection);
    if ($gene_selection =~ /\S/) {
	open QUERY, ">$TMP/$tmp_file_name";
	foreach my $row (@gene_selection) {
	    $row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option. 
	    print QUERY $row, "\n";
	}
	close QUERY;
	&DelayedRemoval("$TMP/$tmp_file_name");
	$parameters .= " -i $TMP/$tmp_file_name";
    } else {
	&cgiError("You should enter at least one gene identifier in the query box..");
    }
}

print  "<PRE><B>Command :</B> $command $parameters</PRE><P>" if ($ECHO >= 1);

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

#    $ENV{RSA_OUTPUT_CONTEXT} = "text";
    open RESULT, "$command $parameters |";

    ### print the result ### 
    &PipingWarning();
    if ($query->param('output') =~ /server/i) {
	&Info("The result will appear below ...");
    }

    print '<H3>Result</H3>';

    ### open the mirror file
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
    ## Choose organism-specific or alternative background
    if ($query->param('single_multi_org') eq 'multi') {
	$oligo_background_model = "<INPUT type=\"hidden\" NAME=\"freq_estimate\" VALUE=\"Residue frequencies from input sequence\">";
	$dyad_background_model = "<INPUT type=\"hidden\" NAME=\"freq_estimate\" VALUE=\"monads\">";
    } else {
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
	$oligo_background_model = "\n<INPUT type=\"hidden\" NAME=\"background\" VALUE=\"$background\">";
	$dyad_background_model = $oligo_background_model;
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
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	$oligo_background_model
	<INPUT type="submit" value="oligonucleotide analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="dyad-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	$dyad_background_model
	<INPUT type="submit" value="dyad analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="position-analysis_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="hidden" NAME="background" VALUE="$background">
	<INPUT type="submit" value="position analysis">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="consensus_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="consensus">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="gibbs_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
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
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
	<INPUT type="hidden" NAME="sequence_file" VALUE="$mirror_file">
	<INPUT type="hidden" NAME="sequence_format" VALUE="$out_format">
	<INPUT type="submit" value="dna-pattern (IUPAC)">
	</FORM>
    </TD>

    <TD>
	<FORM METHOD="POST" ACTION="patser_form.cgi">
	<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
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
