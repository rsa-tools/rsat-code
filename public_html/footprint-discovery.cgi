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

## Limit the analysis to only the 100 first genes
#$parameters .= " -max_genes 2 ";
my $max_genes = 100;
my $max_genename_size = 12;
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
    if ($fields[0] =~ /^>[actg]+$/i){ # Avoid fasta sequences
      &cgiError("Fasta sequences are not valid as input. Please use gene identifiers (eg: YFL021W) or gene names (eg: GAL4, NIL1).<P>\n");
    }elsif ($fields[0] =~ /^[actg]+$/i){
      &cgiError("Sequences format are not valid as input. Please use gene identifiers (eg: YFL021W) or gene names (eg: GAL4, NIL1).<P>\n");      
    }elsif(length($fields[0])>= $max_genename_size){ # put a threshold on the size of the gene name
      &cgiError("The name of the gene is too long and nay not be valid.<P>\n");
    }
    push @query_genes,  $fields[0];
  }

  if (scalar(@query_genes) > $max_genes) {
    &Warning("The RSAT web server cannot handle large scale analysis.<P> The analysis has been restricted to the ", $max_genes, " first genes<P>\n");
    @query_genes= splice(@query_genes,0,$max_genes);
    $query_prefix = $max_genes."_first_genes";
  }else{
    if (scalar(@query_genes) == 1) {
      $query_prefix = $query_genes[0];
    } elsif (scalar(@query_genes) <= 10) {
      $query_prefix = join "_", @query_genes;
    }else {
      $query_prefix = scalar(@query_genes)."_genes";
    }
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
#  print QUERY $query->param('queries');
  print QUERY join("\n",@query_genes);
  close QUERY;
  &DelayedRemoval($query_file);
  $parameters .= " -genes ".$query_file;
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

## Convert assembled patterns to PSSM
if ($query->param('to_matrix')) {
  $parameters .= " -to_matrix";
}

## Background model
$bg_model = $query->param('bg_model');
$parameters .= " -bg_model ".$bg_model;

## Output prefix
$parameters .= " -o ".$file_prefix;

## Report the command
print "<PRE>$command $parameters </PRE>" if ($ENV{rsat_echo} >= 1);

$index_file = $result_subdir."/".$query_prefix."_index.html";
my $mail_title = join (" ", "[RSAT]", "footprint-discovery", $query_prefix, $bg_model, $taxon, $organism, &AlphaDate());
&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);

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

