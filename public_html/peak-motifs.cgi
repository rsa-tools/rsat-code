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

################################################################
## result page header
### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("peak-motifs result", "results");

&ListParameters() if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

## peak-motifs command
$command = "$ENV{RSAT}/perl-scripts/peak-motifs";

################################################################
## configuration

## We need to create the output directory before starting peak-motif
## for uploading the input sequences and reference motifs
$output_dir_prefix = sprintf "peak-motifs.%s", &AlphaDate();
$output_dir_full_path = &RSAT::util::make_temp_file("", $output_dir_prefix, 1, 1); $output_dir = &ShortFileName($tmp_file_path);
$output_prefix = "peak-motifs";
system("mkdir -p $output_dir_full_path; chmod 755 $output_dir_full_path");


################################################################
## command line paramters
### read parameters
$parameters = " -v 1";

### title
if ($query->param('title')){
  my $title = $query->param('title');

  ## Suppress characters that may cause problems when used in file names
  $title =~ s/\s+/_/g;
  $title =~ s/\//_/g;
  $title =~ s/:/_/g;
  $title =~ s/\'/_/g;
  $title =~ s/\"/_/g;
  $parameters .= " -title '".$title."' ";
} else {
  &RSAT::error::FatalError("You must enter a title for this analysis");
}

## Default tasks
@tasks = ("purge",
	  "seqlen",
	  "composition",
	  "disco",
	  "merge_words",
	  "collect_motifs",
	  "motifs_vs_motifs",
	  "timelog",
	  "archive",
	  "synthesis",
	  "small_summary");

################################################################
## Peak sequences file

## Test sequences
($sequence_file, $sequence_format) = &MultiGetSequenceFile(1, $output_dir_full_path."/".$output_prefix."peak_seq", 1);
$parameters .= "-i $sequence_file ";

### control sequences file
($control_file, $sequence_format) = &MultiGetSequenceFile(2, $output_dir_full_path."/".$output_prefix."control_seq", 0);

if ($control_file) {
  $parameters .= "-ctrl $control_file ";
}


## Number of top peaks
if ($query->param('top_sequences')){
  if (&IsNatural($query->param('top_sequences'))) {
    $parameters .= " -top_peaks ".$query->param('top_sequences');
  } else {
    &FatalError("Number of top peaks is incorrect");
  }
}

## Max peak length (clipping)
if ($query->param('max_seq_len')){
  if (&IsNatural($query->param('max_seq_len'))) {
    $parameters .= "  -max_seq_len ".$query->param('max_seq_len')*2; ## here the program needs the length of the fragments, so x2
  } else {
    &FatalError("Incorrect maximal sequence length. Check your parameters for data restriction");
  }
}


################################################################
## Motif discovery parameters

my @disco_algo =(); ## motif disco algorithms to run


## oligo-analysis
if ($query->param('oligo-analysis') =~ /on/) {
  push(@disco_algo, "oligos");

  ## Markov order is specific to oligo-analysis
  if (&IsInteger($query->param('markov'))) {
    $parameters .= " -min_markov ".$query->param('markov')." -max_markov ".$query->param('markov');
  } elsif ($query->param('markov') eq 'auto') {
    $parameters .= " -markov auto";
  }
}

## dyad-analysis
if ($query->param('dyad-analysis') =~ /on/) {
    push(@disco_algo, "dyads");
}

## position-analysis
if ($query->param('position-analysis') =~ /on/) {
    push(@disco_algo, "positions");
}

## local-word-analysis
if ($query->param('local-word-analysis') =~ /on/) {
    push(@disco_algo, "local_words");
}

if ($query->param('local-word-analysis_dyads') =~ /on/) {
  ## TO BE ADDED WHEN THE PROGRAM WILL BE FASTER
}

## Motif discovery algorithms
if (scalar(@disco_algo) >= 1) {
  $parameters .= " -disco ".join(",",@disco_algo);

  ## Number of motifs per algorithm
  if (&IsNatural($query->param('nmotifs'))) {
    $parameters .= " -nmotifs ".$query->param('nmotifs')." ";
  }

  ## Oligonucleotide lengths (used for oligo-analysis, position-analysis and local-word-analysis)
  my @oligo_lengths =();
  foreach my $i (6..7){
    if ($query->param('oligo_length'.$i) =~ /on/){
      push (@oligo_lengths, $i);
    }
  }
  @oligo_lengths = sort @oligo_lengths;
  &RSAT::error::FatalError("Select at least one oligo size for oligo-analysis") if (scalar(@oligo_lengths) < 1 );
  $parameters .= " -minol ".$oligo_lengths[0]." -maxol ".$oligo_lengths[-1]." ";
}


## Strands
if ($query->param('strand')) {
    $parameters .= " ".$query->param('strand')." ";
}

################################################################
## Compare discovered motifs with motif databases
my ($mat_db_params, @selected_db) = &GetMatrixDBfromBox();
if (scalar(@selected_db) > 0) {
  $parameters .= $mat_db_params;
  push(@tasks, "motifs_vs_db");
}

################################################################
## Custom collection motifs (not reference)
if ($query->param('custom_motif_db')) {
  my $persomotif_file = $output_dir_full_path."/".$output_prefix."_custom_motif_db.tf";

  $upload_persomotif = $query->param('custom_motif_db');
  if ($upload_persomotif) {
    my $type = $query->uploadInfo($upload_persomotif)->{'Content-Type'};
    if ($upload_persomotif =~ /\.gz$/) {
      $refmotif_file .= ".gz";
    }
    open REF, ">$persomotif_file" ||
      &cgiError("Cannot store sequence file in temp dir.");
    while (<$upload_persomotif>) {
#      print "<br>UPLOADING REF MOTIFS \t", $_;
      print REF;
    }
    close REF;
  }

  ## Name for the custom motif database
  my $dbname_perso="";
  if ($query->param('custom_motif_db_name')){
    $dbname_perso = $query->param('custom_motif_db_name');

    ## Suppress characters that may cause problems when used in file names
    $dbname_perso =~ s/\s+/_/g;
    $dbname_perso =~ s/\//_/g;
    $dbname_perso =~ s/:/_/g;
  } else {
    $dbname_perso = "personnal_collection";
  }
  $parameters .= " -motif_db ".$dbname_perso." tf ".$persomotif_file;
  push(@tasks, "motifs_vs_db");
}

################################################################
## Custom reference motifs
if ($query->param('ref_motif')) {
  my $refmotif_file = $output_dir_full_path."/".$output_prefix."_ref_motifs.tf";

  $upload_refmotif = $query->param('ref_motif');
  if ($upload_refmotif) {
    my $type = $query->uploadInfo($upload_refmotif)->{'Content-Type'};
    if ($upload_refmotif =~ /\.gz$/) {
      $refmotif_file .= ".gz";
    }
    open REF, ">$refmotif_file" ||
      &cgiError("Cannot store sequence file in temp dir.");
    while (<$upload_refmotif>) {
      print REF;
    }
    close REF;
  }
  $parameters .= " -ref_motifs ".$refmotif_file;
  push(@tasks, "ref_motifs,motifs_vs_ref");
}

################################################################
## Scan sequences to search motif occurrences (matrix-scan-quick).
if ($query->param('matrix-scan-quick') =~ /on/) {
    push(@tasks, "scan");

    ## HERE need to add the pval for the background model for
    ## matrix-scan-quick

  ## Markov order is specific to oligo-analysis
  if (&IsInteger($query->param('scan_markov'))) {
    $parameters .= " -scan_markov ".$query->param('scan_markov');
  }
}

################################################################
## UCSC custom track
if ($query->param('visualize') eq "galaxy") {
  $parameters .= " -source galaxy ";
}
if ($query->param('visualize') eq "bed_coord") {
  ## upload the coord file
  $upload_coord_file = $query->param('bed_file');
  if ($upload_coord_file) {
    my $type = $query->uploadInfo($upload_coord_file)->{'Content-Type'};
    if ($upload_coord_file =~ /\.gz$/) {
      $upload_coord_file .= ".gz";
    }
    open REF, ">$upload_coord_file" ||
      &cgiError("Cannot store sequence file in temp dir.");
    while (<$upload_coord_file>) {
      print REF;
    }
    close REF;
  }

  unless ($query->param('assembly')){
    &cgiError("The assembly version must be provided when using a BED coordinates file.");
  }

  $parameters .= " -coord ".$query->param('assembly')." ".$upload_coord_file;
}

## Add list of tasks
$parameters .= " -task " . join(",", @tasks);

### output prefix
$parameters .= " -prefix $output_prefix";


### other default parmaters
$parameters .= " -noov -img_format png ";

### output directory
$parameters .= " -outdir ".$output_dir_full_path;

&ReportWebCommand($command." ".$parameters);

################################################################
## Display or send result
$index_file = $output_dir_full_path."/".$output_prefix."_synthesis.html";
#$index_file = $output_dir."/".$output_prefix."_synthesis.html";
my $mail_title = join (" ", "[RSAT]", "peak-motifs", &AlphaDate());
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere", $index_file, title=>"$mail_title",no_email=>1);
} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>"$mail_title");
}

################################################################
## result page footer
print $query->end_html;

exit(0);

