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
$random_genes_command = "$SCRIPTS/random-genes";
$tmp_file_name = sprintf "random-genes.%s", &AlphaDate;
$result_file = "$TMP/$tmp_file_name.res";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("Random gene selection result");

#### update log file ####
&UpdateLogFile;

&ListParameters() if ($ECHO >= 2);

#### read parameters ####
$parameters = "";
$gene_nb = $query->param('gene_nb');
if (&IsNatural($gene_nb)) {
    $parameters .= " -n $gene_nb ";
} else {
    &FatalError("Gene number must be a natural number");
}

unless ($organism = $query->param('organism')) {
    &FatalError("You should specify an organism to use non-coding frequency calibration");
}
if (defined(%{$supported_organism{$organism}})) {
    $parameters .= " -org $organism ";
} else {
    &FatalError("Organism $organism is not supported on this site");
}

if ($query->param('replacement')) {
    $parameters .= " -rep ";
}


print "<PRE>command: $random_genes_command $parameters<P>\n</PRE>" if ($ECHO);

### execute the command ###
if ($query->param('output') eq "display") {
    &PipingWarning();

    ### prepare data for piping
    open RESULT, "$random_genes_command $parameters |";

    print '<H2>Result</H2>';
    print '<PRE>';
    while (<RESULT>) {
	print $_;
	$genes .= $_;
    }
    print '</PRE>';
    close(RESULT);

    &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$random_genes_command $parameters", $query->param('user_email'));
}
print $query->end_html;

exit(0);


################################################################
#
# Pipe the result to other commands
#
sub PipingForm {
#    my $genes = `cat $result_file`;

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
