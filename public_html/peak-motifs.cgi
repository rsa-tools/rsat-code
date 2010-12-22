#!/usr/bin/perl

############################################ imports
### redirect error log to a file
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "/tmp/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

############################################ configuration
$command = "$ENV{RSAT}/perl-scripts/peak-motifs";
$output_directory = sprintf "peak-motifs.%s", &AlphaDate();
$output_prefix = "peak-motifs";

$output_path = "$TMP/$output_directory";
$output_path =~ s|\/\/|\/|g;
`mkdir -p $output_path`;
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/perl-scripts" . ":$ENV{RSAT}/python-scripts";
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/bin";

############################################ result page header
### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("peak-motifs result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

### update log file
&UpdateLogFile();

############################################ command line paramters
### read parameters
$parameters = "";

### title
if ($query->param('title')){
	my $title = $query->param('title');
	$title =~ s/\s+/_/g;
	$parameters .= " -title $title ";
}

### peak sequences file
($sequence_file, $sequence_format) = &MultiGetSequenceFile(1, "$output_path/$output_prefix" . "peak_seq", 1);

$parameters .= "-i $sequence_file ";


### tasks

## default
@tasks = ("purge", "seqlen", "composition", "merge_words", "collect_motifs", "timelog", "synthesis");

## motif-disco
my $oligo_params = "";
foreach my $i (6..8){
    if ($query->param('oligo_length'.$i) =~ /on/){
	$oligo_params .= " -l ".$i;
    }
}
$parameters .= $oligo_params;

if ($query->param('oligo-analysis') =~ /on/) {
    push(@tasks, "oligos");
    &FatalError("Select at least one oligo size for oligo-analysis") if ($oligo_params eq "");
}

if ($query->param('dyad-analysis') =~ /on/) {
    push(@tasks, "dyads");
}
if ($query->param('local-word-analysis') =~ /on/) {
    push(@tasks, "local_words");
    &FatalError("Select at least one oligo size for local-word-analysis") if ($oligo_params eq "");
}
if ($query->param('position-analysis') =~ /on/) {
    push(@tasks, "positions");
}

### task specific parameters
if (&IsNatural($query->param('markov'))) {
    $parameters .= " -max_markov ".$query->param('markov')." -min_markov ".$query->param('markov');
}

### restrict the input dataset
if ($query->param('top_sequences')){
    if (&IsNatural($query->param('top_sequences'))) {
	$parameters .= " -top_peaks ".$query->param('top_sequences');
    } else {
	&FatalError("Number of top peaks is incorrect");
    }	
}

if ($query->param('max_seq_len')){
    if (&IsNatural($query->param('max_seq_len'))) {
	$parameters .= "  -max_seq_len ".$query->param('max_seq_len')*2; ## here the program needs the length of the fragments, so x2
    } else {
	&FatalError("Incorrect maximal sequence length. Check your parameters for data restriction");
    }	
}

################################################################
## Motif databases supported on the RSAT Web site
#if ($query->param('compare_motif_db') =~ /on/) {
#    push(@tasks, "motifs_vs_db");

## load the files containing the databases
my ($mat_db_params, @selected_db) = &GetMatrixDBfromBox();
if (scalar(@selected_db) > 0) {
  $parameters .= $mat_db_params;
  push(@tasks, "motifs_vs_db");
}

################################################################
## Custom reference motifs
if ($query->param('ref_motif')) {
  my $refmotif_file = ${TMP}."/".$output_path."/".$output_prefix."ref_motifs.tf";
  my $upload_refmotif = $query->param('ref_motif');
  if ($upload_refmotif) {
    my $type = $query->uploadInfo($upload_refmotif)->{'Content-Type'};
    open FILE, ">$refmotif_file" ||
      &cgiError("Cannot store reference motif file in temp dir.");
    while (<$upload_refmotif>) {
      print FILE;
    }
    close FILE;
    $parameters .= " -ref_motifs PERSONAL_MOTIFS transfac ".${TMP}."/".$output_path."/".$output_prefix."ref_motifs.tf";
  } else {
    &FatalError ("If you want to upload a personal matrix file, you should specify the location of this file on your hard drive with the Browse button");
  }
}

################################################################
## Scan sequences to search motif occurrences (matrix-scan-quick).
if ($query->param('matrix-scan-quick') =~ /on/) {
    # push(@tasks, "scan");

    ## HERE need to add the pval and markov order for the background model for matrix-scan-quick

}

## HERE finish the BED custom track parameters + task
## UCSC custom track
#if ($query->param('bed_custom_track') =~ /on/) {
#       push(@tasks, "?");
#      
#       
#}

### add -task
$parameters .= " -task " . join(",", @tasks);

### output directory
$parameters .= " -outdir $output_path";

### output prefix
$parameters .= " -prefix $output_prefix";

### verbosity
#$parameters .= " -v 1";

### other default parmaters
$parameters .= " -2str -noov -img_format png ";

############################################ display or send result
$index_file = $output_directory."/".$output_prefix."_synthesis.html";
my $mail_title = join (" ", "[RSAT]", "peak-motifs", &AlphaDate());
&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);
# $debug = "$command $parameters 2> $TMP/log.txt";
# print $debug;
# `$debug`;

############################################ result page footer
print $query->end_html;

exit(0);

