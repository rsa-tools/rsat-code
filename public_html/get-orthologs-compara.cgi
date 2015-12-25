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
&RSA_header("get-orthologs-compara result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/get-orthologs-compara";
$prefix = "get-orthologs-compara";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); 
$tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();


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
#generic organism names are not listed in organisms.tab
#unless (%{$supported_organism{$organism}}) {
#    &cgiError("Organism $org is not supported on this site");
#}
$parameters .= " -ref_org $organism";



## ##############################################################
## Thresholds and type
my $type = "";
unless ($type = $query->param('type')) {
    &cgiError("You should select an homology type");
}
$parameters .= " -type $type";

my $ident_target =  $query->param('ident_target');
$ident_target = 0 if (lc($ident_target) eq "none");
unless (&IsReal($ident_target) && ($ident_target >= 0) && ($ident_target <= 100)) {
    &cgiError("Invalid value for target identity: should be a Real number between 0 and 100.");
}
$parameters .= " -ident_target $ident_target";

my $ident_query =  $query->param('ident_query');
$ident_query = 0 if (lc($ident_query) eq "none");
unless (&IsReal($ident_query) && ($ident_query >= 0) && ($ident_query <= 100)) {
    &cgiError("Invalid value for query identity: should be a Real number between 0 and 100.");
}
$parameters .= " -ident_query $ident_query";

# my $ident_query = "";
# unless ($ident_query = $query->param('ident_query')) {
#     &cgiError("Please set ident_query");
# }
# &cgiError("Please choose a [0,100] %value for ident_query") if($ident_query <0 || $ident_query > 100);
# $parameters .= " -ident_query $ident_query";


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

