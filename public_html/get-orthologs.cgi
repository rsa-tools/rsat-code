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
$command = "$SCRIPTS/get-orthologs";
$tmp_file_name = sprintf "get-orthologs.%s", &AlphaDate();
$result_file = "$TMP/$tmp_file_name.res";

$ECHO=2;

### Read the CGI query
$query = new CGI;

#### update log file ####
&UpdateLogFile();

### Print the header
&RSA_header("get-orthologs result");

#### read parameters ####
$parameters .= "";

################################################################
#### queries
if ( $query->param('queries') =~ /\S/) {
    open QUERY, ">$TMP/$tmp_file_name";
    print QUERY $query->param('queries');
    close QUERY;
    &DelayedRemoval("$TMP/$tmp_file_name");
    $parameters .= " -i $TMP/$tmp_file_name";
} else {
    &cgiError("You should enter at least one query in the box\n");
}  

################################################################
#### organism
my $organism = "";
unless ($organism = $query->param('organism')) {
    &cgiError("You should specify a query organism");
}
unless (defined(%{$supported_organism{$organism}})) {
    &cgiError("Organism $org is not supported on this site");
}
$parameters .= " -org $organism";


################################################################
#### Taxon
my $taxon = "";
unless ($taxon = $query->param('taxon')) {
    &cgiError("You should specify a taxon");
}
$parameters .= " -taxon $taxon";


## ##############################################################
## Thresholds
my @parameters = $query->param;
foreach my $param (@parameters) {
    if ($param =~ /^lth_(.+)/) {
	my $field = $1 ;
	my $value = $query->param($param);
	next unless (&IsReal($value));
	$parameters .= " -lth ".$field." ".$value;
    } elsif ($param =~ /^uth_(.+)/) {
	my $field = $1 ;
	my $value = $query->param($param);
	next unless (&IsReal($value));
	$parameters .= " -uth ".$field." ".$value;
    }
}


## Report the command
print "<PRE>$command $parameters </PRE>" if ($ECHO);

################################################################
#### run the command
if ($query->param('output') eq "display") {
    &PipingWarning();

    open RESULT, "$command $parameters |";

    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);

    &PipingForm();

    print "<HR SIZE = 3>";
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else { 
    &EmailTheResult("$command $parameters", $query->param('user_email'));
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
<TABLE>
<TR>
<TD>
<H3>Next step</H3>
</TD>
<TD>
<FORM METHOD="POST" ACTION="retrieve-seq_form.cgi">
<INPUT type="hidden" NAME="organism" VALUE="$organism">
<INPUT type="hidden" NAME="multigenome" VALUE="OK">
<INPUT type="hidden" NAME="genes" VALUE="selection">
<INPUT type="hidden" NAME="gene_selection" VALUE="$genes">
<INPUT type="submit" value="retrieve sequences">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form

}

