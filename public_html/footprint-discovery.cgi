#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push @INC, "$`lib/";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/footprint-discovery";

#$ENV{rsat_echo}=2;

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("footprint-discovery result", "results");

#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

#### read parameters ####
$parameters = " -v 1 -index ";

################################################################
#### Compute the query prefix
my $query_prefix = "footprints";
if ($query->param('queries') =~ /\S/) {
  my @query_lines = split "\n", $query->param('queries');
  my $l = 0;
  foreach my $line (@query_lines) {
    $l++;
    $line =~ s/^\s+//;
    my @fields = split /\s+/, $line;
    push @query_genes,  $fields[0];
  }

  if (scalar(@query_genes) == 1) {
    $query_prefix = $query_genes[0];
  } elsif (scalar(@query_genes) <= 10) {
    $query_prefix = join "_", @query_genes;
  } else {
    $query_prefix = scalar(@query_genes)." genes";
  }
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

################################################################
## File prefix
$tmp_file_name = join( "_", "footprint-discovery", $taxon, $organism, $query_prefix, &AlphaDate());
$result_subdir = $tmp_file_name;
$result_dir = $TMP."/".$result_subdir;
$result_dir =~ s|\/\/|\/|g;
`mkdir -p $result_dir`;
$file_prefix = $result_dir."/".$query_prefix;
$query_file = $file_prefix."_genes";

################################################################
## Prepare a file on the server with the query genes
if ($query->param('queries') =~ /\S/) {
  open QUERY, ">".$query_file;
  print QUERY $query->param('queries');
  close QUERY;
  &DelayedRemoval($query_file);
  $parameters .= " -i ".$query_file;
} else {
  &cgiError("You should enter at least one query in the box\n");
}

## Return fields and threshold values for dyad-analysis
&CGI_return_fields();

## Infer operon leader genes
if ($query->param('leaders')) {
  $parameters .= " -infer_operons";
}

## Dyad filter
if (!$query->param('dyads_filter')) {
  $parameters .= " -no_filter";
}

## Background model
$parameters .= " -bg_model ".$query->param('bg_model');

## Output prefix
$parameters .= " -o ".$file_prefix;

## Report the command
print "<PRE>$command $parameters </PRE>" if ($ENV{rsat_echo} >= 1);

$index_file = $result_subdir."/".$query_prefix."_index.html";
my $mail_title = join (";", "RSAT", "footprint-discovery", &AlphaDate(), $taxon, $organism, $query_prefix);
&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file,
		title=>$mail_title,
	       );

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

