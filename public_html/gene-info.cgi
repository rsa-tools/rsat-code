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
$command = "$SCRIPTS/gene-info";
$tmp_file_name = sprintf "gene-info.%s", &AlphaDate;
$result_file = "$TMP/$tmp_file_name.res";

### Read the CGI query
$query = new CGI;

#### update log file ####
&UpdateLogFile;

### Print the header
&RSA_header("gene-info result");

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
### feature type
if ($query->param('feattype')) {
    my ($feattype) = split " ", $query->param('feattype'); ### take the first word
    $parameters .= " -feattype ".$feattype;
}

################################################################
#### full match
$parameters .= " -full" if ($query->param('full'));

################################################################
#### match queries against description
$parameters .= " -desc" if ($query->param('match_description'));

################################################################
#### organism
my $organism = "";
unless ($organism = $query->param('organism')) {
    &cgiError("You should specify an organism to use intergenic frequency calibration");
}
unless (defined(%{$supported_organism{$organism}})) {
    &cgiError("Organism $org is not supported on this site");
}
$parameters .= " -org $organism";


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
<INPUT type="hidden" NAME="genes" VALUE="selection">
<INPUT type="hidden" NAME="gene_selection" VALUE="$genes">
<INPUT type="submit" value="retrieve sequences">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form

}

