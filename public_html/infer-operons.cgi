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
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("infer-operon result", 'results');

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);


$prefix = "infer-operons";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "infer-operon.%s", &AlphaDate();

@result_files = ();


$parameters = "";

################################################################
## Single or multi-genome query
# if ($query->param('single_multi_org') eq 'multi') {
#     $command = "$SCRIPTS/infer-operon-multigenome";

#     &cgiMessage(join("<P>",
# 		     "The computation can take a more or less important time depending on the taxon size.",
# 		     "If the answer does not appear in due time, use the option <i>output email</i>"));
# } else {
     $command = "$SCRIPTS/infer-operon";

# }



#### organism
$organism = $query->param('organism');
if (defined($supported_organism{$organism})) {
    $organism_name = $supported_organism{$organism}->{'name'};

    $parameters .= " -org ".$organism unless ($query->param('single_multi_org') eq 'multi'); ## For multi-genome retrieval, the query organism name is passed to motif discovery programs, but it is not necessary
} else {
    &cgiError("Organism '",
	      $organism,
	      "' is not supported on this web site.");
}



#### distance threshold
my $dist_thr = $query->param('dist_thr');
&RSAT::error::FatalError($dist_thr, "Invalid value for distance threshold. Should be an Integer value.") 
  unless (&IsInteger($dist_thr));
$parameters .= " -dist ".$dist_thr;

#### Min gene number
my $min_gene_nb = $query->param('min_gene_nb');
&RSAT::error::FatalError($min_gene_nb, "Invalid value for min gene number. Should be a strictly positive Natural number.") 
  unless ((&IsNatural($min_gene_nb)) && ($min_gene_nb > 0));
$parameters .= " -min_gene_nb ".$min_gene_nb;



### return fields
my $i=0;
foreach my $field ("query", "name", "leader","trailer","operon", "upstr_dist", "q_info","up_info","down_info", "gene_nb") {
    my $return_field = "return_".$field;
#    my $return_field = $field;
    if ($query->param($return_field) eq "on"){
	$parameters .= " -return ".$field;
	$i++;
    }
}
&cgiError("Invalid output fields, please check at least one output field.") if ($i==0);

#### queries ####
if ($query->param('genes') eq "all") {
    ### take all genes as query
    $parameters .= " -all ";
} elsif ($query->param('uploaded_file')) {
    $upload_file = $query->param('uploaded_file');
    $gene_list_file = "${TMP}/${tmp_file_name}.genes";
    push (@result_files, 'input genes', $gene_list_file);
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

## Output file
$result_file = $tmp_file_path.".tab";
push (@result_files, 'operons', $result_file);

&ReportWebCommand($command." ".$parameters);


################################################################
#### run the command
if ($query->param('output') eq "display") {
    &PipingWarning();

    print '<H2>Result</H2>';
    open RESULT, "$command $parameters |";
    &PrintHtmlTable(RESULT, $result_file, 1, 5000);
    close(RESULT);

    &PrintURLTable(@result_files);

    &PipingForm();

    print "<HR SIZE = 3>";

} else { 
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}
print $query->end_html();


exit(0);


################################################################
#
# Pipe the result to other commands
#
sub PipingForm {
    my $genes = `cat $result_file`;

    ### prepare data for piping
    print <<End_of_form;
<HR SIZE = 3>
<TABLE class = 'nextstep'>
<tr>
<td>
<H3>Next step</H3>
</td>
</tr>
<tr>
<TD>
<FORM METHOD="POST" ACTION="retrieve-seq_form.cgi">
<INPUT type="hidden" NAME="organism" VALUE="$organism">
<INPUT type="hidden" NAME="genes" VALUE="selection">
<INPUT type="hidden" NAME="gene_selection" VALUE="$genes">
<INPUT type="hidden" NAME="feattype" VALUE="$feattype">
<INPUT type="submit" value="retrieve sequences">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form
}
