#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push @INC, "$`lib/";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";
require "footprint.lib.pl";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("footprint-discovery result", "results");


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/footprint-discovery";

#$ENV{rsat_echo}=2;

#### read parameters ####
$parameters = " -v 1";

################################################################
## Tasks: some tasks are not supported on the Web interface:
##
## -task network because it requires too much computation if a user
##   introduces all the genes of an organism with the option -sep_genes
##   (required for network inference)
##
## -task operon because the option has not been introduced in the Web
##   form (could be fixed some day).
$tasks .= " -task query_seq,filter_dyads,orthologs,ortho_seq,purge,dyads,maps,gene_index,index";

## Limit the analysis to only the 100 first genes
#$parameters .= " -max_genes 2 ";
my $max_genes = 100;
#my $max_genename_size = 12;

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
#    if ($fields[0] =~ /^>[actg]+$/i){ # Avoid fasta sequences
#      &cgiError("Fasta sequences are not valid as input. Please use gene identifiers (eg: YFL021W) or gene names (eg: GAL4, NIL1).<P>\n");
#    }elsif ($fields[0] =~ /^[actg]+$/i){
#      &cgiError("Sequences format are not valid as input. Please use gene identifiers (eg: YFL021W) or gene names (eg: GAL4, NIL1).<P>\n");      
#    }elsif(length($fields[0])>= $max_genename_size){ # put a threshold on the size of the gene name
#      &cgiError("The name of the gene is too long and nay not be valid.<P>\n");
#    }
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
local $organism_name = "";
unless ($organism_name = $query->param('organism')) {
    &cgiError("You should specify a query organism");
}
unless (%{$supported_organism{$organism_name}}) {
    &cgiError("Organism $org is not supported on this site");
}
$parameters .= " -org $organism_name";


################################################################
#### Taxon
local $taxon = "";
unless ($taxon = $query->param('taxon')) {
    &cgiError("You should specify a taxon");
}
$parameters .= " -taxon $taxon";

################################################################
## Create output directory. This must be done after having read the
## organism and taxon, in order to include these in the path.
$tmp_file_name = join( "_", "footprint-discovery", $taxon, $organism_name, $query_prefix, &AlphaDate());
$result_subdir = $tmp_file_name;
$result_dir = &RSAT::util::make_temp_file("", $result_subdir, 1, 1);
#$result_prefix = "footprint-discovery";
system("mkdir -p $result_dir; chmod 755 $result_dir");

#`mkdir -p $result_dir`;
#$file_prefix = $result_dir."/".$query_prefix;
$query_file = $result_dir."/".$query_prefix."_genes";

# &RSAT::message::Debug("tmp_file_name=".$tmp_file_name,
# 		      "<br>result_subdir=".$result_subdir,
# 		      "<br>result_dir=".$result_dir,
# #		      "<br>file_prefix=".$file_prefix,
# 		      "<br>query_file=".$query_file,
# 		     ) if ($main::verbose >= 10);

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
  $tasks.=",operons";
}

## Dyad filter
if ($query->param('dyads_filter')) {
  $parameters .= " -filter";
} else {
  $parameters .= " -no_filter";
}

## Convert assembled patterns to PSSM
if ($query->param('to_matrix')) {
  $parameters .= " -to_matrix";
}

## Background model
local $bg_model = $query->param('bg_model');
$parameters .= " -bg_model ".$bg_model;

## Add Tasks
$parameters .= $tasks;

## Output prefix
$parameters .= " -o ".$result_dir;

## Report the command
&ReportWebCommand($command." ".$parameters);

#$index_file = $result_subdir."/";
$index_dir = join ("/", $result_dir, $taxon, $organism_name);
$index_file = $index_dir."/";
$index_file .= (&MainIndexFileName())[0];

#$index_file .= join("_", $taxon, $organism_name, "bg", $bg_model, "result_index.html");
my $mail_title = join (" ", "[RSAT]", "footprint-discovery", $query_prefix, $bg_model, $taxon, $organism_name, &AlphaDate());
my $log_file = $result_subdir."/server_log.txt";
&EmailTheResult("$command $parameters", $query->param('user_email'), $log_file, index=>$index_file, title=>$mail_title);

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
<INPUT type="hidden" NAME="organism" VALUE="$organism_name">
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

