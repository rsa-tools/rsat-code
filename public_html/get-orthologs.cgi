#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push @INC, "$`lib/";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

#$ENV{rsat_echo}=2;

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("get-orthologs result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/get-orthologs";
$prefix = "get-orthologs";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();


&RSAT::message::Warning("The computation can take a more or less important time depending on the taxon size.",
			"If the answer does not appear in due time, use the option <i>output email</i>");


#### read parameters ####
$parameters = " -v 1";

################################################################
#### queries
if ( $query->param('queries') =~ /\S/) {
  $query_file = $tmp_file_path."_query.txt";
  push @result_files, ("query",$query_file);

  open QUERY, ">".$query_file;
  print QUERY $query->param('queries');
  close QUERY;
  &DelayedRemoval($query_file);
  $parameters .= " -i ".$query_file;
} else {
    &cgiError("You should enter at least one query in the box\n");
}  

################################################################
#### organism
my $organism = "";
unless ($organism = $query->param('organism')) {
    &cgiError("You should specify a query organism");
}
unless (%{$supported_organism{$organism}}) {
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
my @parameters = $query->param();
foreach my $param (@parameters) {
    if ($param =~ /^return_(.+)/) {
	my $field = $1;
	$parameters .= " -return ".$field;
    } elsif ($param =~ /^ortho_lth_(.+)/) {
	my $field = $1 ;
	my $value = $query->param($param);
	next unless (&IsReal($value));
	$parameters .= " -lth ".$field." ".$value;
    } elsif ($param =~ /^ortho_uth_(.+)/) {
	my $field = $1 ;
	my $value = $query->param($param);
	next unless (&IsReal($value));
	$parameters .= " -uth ".$field." ".$value;
    }
}
## Show only resulst for organisms for which we have a blast file
$parameters .= " -nowarn";

## Output file
$result_file = $tmp_file_path.".tab";
push @result_files, ("Orthology table",$result_file);

## Report the command
&ReportWebCommand($command." ".$parameters);

################################################################
#### run the command
if ($query->param('output') eq "display") {
    &PipingWarning();

    open RESULT, "$command $parameters |";
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);

    &PrintURLTable(@result_files);
    &PipingForm();

    print "<HR SIZE = 3>";

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'), $result_file);

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
<TR>

<TD>
<H3>Next step</H3>
</TD>

</tr>
<tr>
<TD>
<FORM METHOD="POST" ACTION="retrieve-seq_form.cgi">
<INPUT type="hidden" NAME="organism" VALUE="$organism">
<INPUT type="hidden" NAME="taxon" VALUE="$taxon">
<INPUT type="hidden" NAME="single_multi_org" VALUE="multi">
<INPUT type="hidden" NAME="seq_label" VALUE="gene identifier + organism + gene name">
<INPUT type="hidden" NAME="genes" VALUE="selection">
<INPUT type="hidden" NAME="gene_selection" VALUE="$genes">
<INPUT type="hidden" NAME="ids_only" VALUE="checked">
<INPUT type="submit" value="retrieve sequences">
</FORM>
</TD>
</TR>
</TABLE>
End_of_form

}

