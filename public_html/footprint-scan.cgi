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

#$ENV{rsat_echo}=2;

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("footprint-scan result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/footprint-scan";

#### read parameters ####
$parameters = " -v 2 -synthesis  -sep_genes ";

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
      &cgiError("The name of the gene is too long and may not be valid.<P>\n");
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
## Create output directory. This must be done after having read the
## organism and taxon, in order to include these in the path.
$tmp_file_name = join( "_", "footprint-scan", $taxon, $organism_name, $query_prefix, &AlphaDate());
$result_subdir = $tmp_file_name;
$result_dir = &RSAT::util::make_temp_file("", $result_subdir, 1, 1);
$result_prefix = "footprint-discovery";
system("mkdir -p $result_dir ; chmod 755 $result_dir ");
#$tmp_file_name = join( "_", "footprint-scan", $taxon, $organism, $query_prefix, &AlphaDate());
#$result_subdir = $tmp_file_name;
#$result_dir = $TMP."/".$result_subdir;
#$result_dir =~ s|\/\/|\/|g;
#`mkdir -p $result_dir`;
$file_prefix = $result_dir."/".$query_prefix;
$query_file = $file_prefix."_genes";


################################################################
#### Matrix specification
$matrix_file = $result_dir."/input_matrix";
local $input_format = lc($query->param('matrix_format'));

if ($query->param('matrix')) {
  open MAT, "> $matrix_file";
  print MAT $query->param('matrix');
  close MAT;
  &DelayedRemoval($matrix_file);
  ($input_format) = split (/\s+/, $input_format);
  $parameters .= " -m $matrix_file";
} else {
  &RSAT::error::FatalError('You did not enter any data in the matrix box');
}
$parameters .=  " -matrix_format " . $input_format;

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

## Return fields and threshold values for footprin-scan
# &CGI_return_fields();

## Infer operon leader genes
if ($query->param('leaders')) {
  $parameters .= " -infer_operons";
}

## info_lines
if ($query->param('info_lines')) {
  $parameters .= " -info_lines ";
}

################################################################
## Background model

## Markov order
my $markov_order = $query->param('markov_order');
&RSAT::error::FatalError("Markov model should be a Natural number") unless &IsNatural($markov_order);

## Method for specifyung the background model
my $bg_method = $query->param('bg_method');
if ($bg_method eq "bginput") {
  $parameters .= " -bginput";
  $parameters .= " -markov ".$markov_order;

} elsif ($bg_method eq "window") {
  my $window_size = $query->param('window_size');
  &RSAT::message::FatalError(join("\t",$window_size, "Invalid value for the window size. Should be a Natural number." )) unless (&IsNatural($window_size));

  $parameters .= " -window ".$window_size;
  $parameters .= " -markov ".$markov_order;

} elsif ($bg_method eq "bgfile") {
  ## Select pre-computed background file in RSAT genome directory
  my $organism_name = $query->param("organism");
  my $noov = "ovlp";
  my $background_model = $query->param("background");
  my $oligo_length = $markov_order + 1;
  $bg_file = &ExpectedFreqFile($organism_name,
			       $oligo_length, $background_model,
			       noov=>$noov, str=>"-1str");

  $parameters .= " -bgfile ".$bg_file.".gz";

} elsif ($bg_method =~ /upload/i) {
  ## Upload user-specified background file
  my $bgfile = "${TMP}/${tmp_file_name}_bgfile.txt";
  my $upload_bgfile = $query->param('upload_bgfile');
  if ($upload_bgfile) {
    if ($upload_bgfile =~ /\.gz$/) {
      $bgfile .= ".gz";
    }
    my $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
    open BGFILE, ">$bgfile" ||
      &cgiError("Cannot store background file in temp dir.");
    while (<$upload_bgfile>) {
      print BGFILE;
    }
    close BGFILE;
    $parameters .= " -bgfile $bgfile";
    $parameters .= " -bg_format ".$query->param('bg_format');
  } else {
    &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
  }

} else {
  &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
}

## Pseudo-frequency for the background model
if (&IsReal($query->param('bg_pseudo'))) {
  $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
}


## Image format
if ($query->param('img_format')) {
    $image_format = $query->param('img_format');
} else {
    $image_format = $ENV{rsat_img_format} || "png";
}
$parameters .= " -plot_format ".$image_format;
$parameters .= " -map_format ".$image_format;

## Output prefix
$parameters .= " -o ".$file_prefix;

## Report the command
&ReportWebCommand($command." ".$parameters);

$index_file = $result_dir."/".$query_prefix."/".$taxon."/".$organism."/all_matrices_report.html";
my $mail_title = join (" ", "[RSAT]", "footprint-scan", $query_prefix, $bg_model, $taxon, $organism, &AlphaDate());
#&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);
my $log_file = $result_subdir."/server_log.txt";
&EmailTheResult("$command $parameters", $query->param('user_email'), $log_file, index=>$index_file, title=>$mail_title);

print $query->end_html();

exit(0);


################################################################
#
# Pipe the result to other commands
#
# sub PipingForm {
#     my $genes = `cat $result_file`;
#     ### prepare data for piping
#     print <<End_of_form;
# <HR SIZE = 3>
# <TABLE class = 'nextstep'>
# <TR>

# <TD>
# <H3>Next step</H3>
# </TD>

# </tr>
# <tr>
# <TD>
# <FORM METHOD="POST" ACTION="retrieve-seq_form.cgi">
# <INPUT type="hidden" NAME="organism" VALUE="$organism">
# <INPUT type="hidden" NAME="single_multi_org" VALUE="multi">
# <INPUT type="hidden" NAME="seq_label" VALUE="gene identifier + organism + gene name">
# <INPUT type="hidden" NAME="genes" VALUE="selection">
# <INPUT type="hidden" NAME="gene_selection" VALUE="$genes">
# <INPUT type="hidden" NAME="ids_only" VALUE="checked">
# <INPUT type="submit" value="retrieve sequences">
# </FORM>
# </TD>
# </TR>
# </TABLE>
# End_of_form

# }

