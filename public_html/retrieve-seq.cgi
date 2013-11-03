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
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("retrieve-seq result", 'results');


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$prefix = "retrieve-seq";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1,0); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();

$parameters = "";

## ##############################################################
## The query contains only identifiers (no need to load synonyms)
if ($query->param('ids_only')) {
    $parameters .= " -ids_only";
}

################################################################
## Single or multi-genome query
if ($query->param('single_multi_org') eq 'multi') {
    $command = "$SCRIPTS/retrieve-seq-multigenome";

    &cgiMessage(join("<P>",
		     "The computation can take more or less time depending on the taxon size.",
		     "If the answer does not appear in due time, use the option <i>output email</i>"));
} else {
    $command = "$SCRIPTS/retrieve-seq";

}



#### organism
$organism = $query->param('organism');
if (defined($supported_organism{$organism})) {
    $organism_name = $supported_organism{$organism}->{'name'};

    $parameters .= " -org ".$organism unless ($query->param('single_multi_org') eq 'multi'); ## For multi-genome retrieval, the query organism name is passed to motif discovery programs, but it is not necessary
} else {
    &cgiError("Organism '",
	      $organism,
	      "' is not supported on this web site.");
}


### feature type
if ($query->param('feattype')) {
    my ($feattype) = split " ", $query->param('feattype'); ### take the first word
    $parameters .= " -feattype ".$feattype;
}

### sequence type
if ($query->param('sequence_type')) {
    ($seq_type) = split " ", $query->param('sequence_type'); ### take the first word
    $parameters .= " -type ".$seq_type;
}

## Output sequence format
$out_format = lc($query->param('format'));
if ($accepted_output_seq{$out_format}) {
  $parameters .= " -format $out_format";
}

### sequence label
$seq_label = lc($query->param('seq_label'));
if ($seq_label eq 'gene identifier') {
    $parameters .= " -label id";
} elsif ($seq_label eq 'gene name') {
    $parameters .= " -label name";
} elsif ($seq_label eq 'gene identifier + name') {
    $parameters .= " -label id,name";
} elsif ($seq_label eq 'gene identifier + organism + gene name') {
    $parameters .= " -label id,organism_name,name";
} elsif ($seq_label eq 'full identifier') {
    $parameters .= " -label id,name,organism_name,sequence_type,current_from,current_to,ctg,orf_strand,reg_left,reg_right";
} else {
    &cgiError("Invalid option for sequence label '$seq_label'");
}

### limits ###
if (&IsInteger($query->param('from'))) {
    $parameters .= " -from ".$query->param('from');
}
if (&IsInteger($query->param('to'))) {
    $parameters .= " -to ".$query->param('to');
}


### prevent orf overlap ###
if (lc($query->param('noorf')) eq "on") {
  $noorf = 1; ## For the piping form
  $parameters .= " -noorf ";
}

### repeat masking
if (lc($query->param('rm')) eq "on") {
    $parameters .= " -rm ";
}

### prevent orf overlap ###
if (lc($query->param('imp_pos')) eq "on") {
    $imp_pos = 1;
    $parameters .= " -imp_pos ";
}


#### queries ####
if ($query->param('genes') eq "all") {
  ### take all genes as query
  $parameters .= " -all ";
} else {
  my $gene_list_file = $tmp_file_path.".genes";
  push @result_files, ("query genes",$gene_list_file);
  if ($query->param('uploaded_file')) {
    $upload_file = $query->param('uploaded_file');
    #    $gene_list_file = "${TMP}/${tmp_file_name}.genes";
    if ($upload_file =~ /\.gz$/) {
      $gene_list_file .= ".gz";
    }
    $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open GENE_FILE, ">".$gene_list_file ||
      &cgiError("Cannot store gene list file in temporary directory");
    while (<$upload_file>) {
      print GENE_FILE;
    }
    close GENE_FILE;
  } else {
    my $gene_selection = $query->param('gene_selection');
    $gene_selection =~ s/\r/\n/g;
    my @gene_selection = split (/[\n\r]/, $gene_selection);
    if ($gene_selection =~ /\S/) {
      open QUERY, ">".$gene_list_file;
      foreach my $row (@gene_selection) {
	next unless $row =~ /\S/; ## Skip empty rows
	chomp($row); ## Suppress newline character
	$row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option
	print QUERY $row, "\n";
      }
      close QUERY;
    } else {
      &cgiError("You should enter at least one gene identifier in the query box..");
    }
  }
  $parameters .= " -i ".$gene_list_file;
  &DelayedRemoval($gene_list_file);
}
&ReportWebCommand($command." ".$parameters);
$sequence_file = "$tmp_file_path.".$out_format;
push @result_files, ("sequences", $sequence_file);

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

#    $ENV{RSA_OUTPUT_CONTEXT} = "text";
    open RESULT, "$command $parameters |";

    ### print the result
    &PipingWarning();

    ### open the sequence file on the server
    if (open MIRROR, ">$sequence_file") {
	$mirror = 1;
	&DelayedRemoval($sequence_file);
    }

    print "<PRE>";
    while (<RESULT>) {
	print "$_" unless ($query->param('output') =~ /server/i);
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close RESULT;
    close MIRROR if ($mirror);

#    if ($query->param('output') =~ /server/i) {
      &PrintURLTable(@result_files);
#      $result_URL = "$ENV{rsat_www}/tmp/${tmp_file_name}.res";
#      print ("The result is available at the following URL: ", "\n<br>",
#	     "<a href=${result_URL}>${result_URL}</a>",
#	     "<p>\n");
#    }

    ### prepare data for piping
    &PipingFormForSequence();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $sequence_file);
}

print $query->end_html;

exit(0);
