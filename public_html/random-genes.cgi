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
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
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
&RSA_header("Random gene selection result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/random-genes";
$prefix = "random-genes";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();

#### read parameters ####
$parameters = "";

#### number of genes
$gene_nb = $query->param('gene_nb');
if (&IsNatural($gene_nb)) {
    $parameters .= " -n $gene_nb ";
} else {
    &FatalError("Gene number must be a natural number");
}

#### number of groups
$group_nb = $query->param('group_nb');
if (&IsNatural($group_nb)) {
    $parameters .= " -g $group_nb ";
} else {
    &FatalError("Number of groups must be a natural number");
}

#### organism
unless ($organism = $query->param('organism')) {
    &FatalError("You should specify an organism to use upstream frequency calibration");
}
if (%{$supported_organism{$organism}}) {
    $parameters .= " -org $organism ";
} else {
    &FatalError("Organism $organism is not supported on this site");
}

### feature type
my ($feattype) = split " ", $query->param('feattype'); ### take the first word
if ($feattype) {
    $parameters .= " -feattype ".$feattype;
}


### Replacement
if ($query->param('replacement')) {
    $parameters .= " -rep ";
}


## Output file
$result_file = $tmp_file_path.".txt";
push @result_files, ("random genes",$result_file);

&ReportWebCommand($command." ".$parameters);

### execute the command ###
if ($query->param('output') eq "display") {

    &PipingWarning();

    ### prepare data for piping
    open RESULT, "$command $parameters |";

    print '<H2>Result</H2>';
    print '<PRE>';
    while (<RESULT>) {
	print $_;
	$genes .= $_;
    }
    print '</PRE>';
    close(RESULT);

    open RES_FILE, ">".$result_file;
    print RES_FILE $genes;
    close RES_FILE;

## STILL TO BE DONE: WRITE A COPY OF GENE LIST IN A FILE
    &PrintURLTable(@result_files);
    &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
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
