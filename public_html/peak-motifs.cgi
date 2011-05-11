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
## configuration
$command = "$ENV{RSAT}/perl-scripts/peak-motifs";
$output_dir = sprintf "peak-motifs.%s", &AlphaDate();
$output_prefix = "peak-motifs";

## We need to create the output directory before starting peak-motif  for uploading the input sequences and reference motifs
$output_path = $TMP."/".$output_dir;
$output_path =~ s|\/\/|\/|g;
system("mkdir -p $output_path");
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/perl-scripts" . ":$ENV{RSAT}/python-scripts";
# $ENV{'PATH'} = $ENV{'PATH'} . ":$ENV{RSAT}/bin";

################################################################
## result page header
### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("peak-motifs result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

### update log file
&UpdateLogFile();

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
  $parameters .= " -title '".$title."' ";
}

## default
@tasks = ("purge", "seqlen", "composition", "disco", "collect_motifs","timelog", "archive", "synthesis");


### peak sequences file
($sequence_file, $sequence_format) = &MultiGetSequenceFile(1, $output_path."/".$output_prefix."peak_seq", 1);

$parameters .= "-i $sequence_file ";

### control sequences file
($control_file, $sequence_format) = &MultiGetSequenceFile(2, $output_path."/".$output_prefix."control_seq", 0);

if ($control_file){
	$parameters .= "-ctrl $control_file ";
}



### tasks

## motif-disco
my $oligo_params = "";
my @oligo_lengths =();
foreach my $i (6..8){
    if ($query->param('oligo_length'.$i) =~ /on/){
		push (@oligo_lengths, $i);
    }
}
@oligo_lengths = sort @oligo_lengths;
$oligo_params .= " 	-minol ".$oligo_lengths[0]." -maxol ".$oligo_lengths[-1]." ";
$parameters .= $oligo_params;

## motif disco algo
my @disco_algo =();
if ($query->param('oligo-analysis') =~ /on/) {
    push(@disco_algo, "oligos");
    &FatalError("Select at least one oligo size for oligo-analysis") if ($oligo_params eq "");
}

if ($query->param('dyad-analysis') =~ /on/) {
    push(@disco_algo, "dyads");
}
if ($query->param('local-word-analysis') =~ /on/) {
    push(@disco_algo, "local_words");
    &FatalError("Select at least one oligo size for local-word-analysis") if ($oligo_params eq "");
}
if ($query->param('position-analysis') =~ /on/) {
    push(@disco_algo, "positions");
}
if ($query->param('local-word-analysis_dyads') =~ /on/) {
	## to add
}

$parameters .= " -disco ".join(",",@disco_algo);

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
## Custom collection motifs (not reference)
if ($query->param('perso_motif')) {
  my $persomotif_file = $output_path."/".$output_prefix."_perso_motif.tf";

  $upload_persomotif = $query->param('perso_motif');
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
  
  ### name
    my $dbname_perso="";
if ($query->param('perso_motif_name')){
  $dbname_perso = $query->param('perso_motif_name');

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
  my $refmotif_file = $output_path."/".$output_prefix."_ref_motifs.tf";

  $upload_refmotif = $query->param('ref_motif');
  if ($upload_refmotif) {
    my $type = $query->uploadInfo($upload_refmotif)->{'Content-Type'};
    if ($upload_refmotif =~ /\.gz$/) {
      $refmotif_file .= ".gz";
    }
    open REF, ">$refmotif_file" ||
      &cgiError("Cannot store sequence file in temp dir.");
    while (<$upload_refmotif>) {
#      print "<br>UPLOADING REF MOTIFS \t", $_;
      print REF;
    }
    close REF;
  }

#  my $upload_refmotif = $query->param('ref_motif');
#  #  if ($upload_refmotif) {
#  my $type = $query->uploadInfo($upload_refmotif)->{'Content-Type'};
#  open FILE, ">$refmotif_file" ||
#    &cgiError("Cannot store reference motif file in temp dir.");
#  while (<$upload_refmotif>) {
#    print FILE;
#  }
#  close FILE;
#  #    $parameters .= " -ref_motifs PERSONAL_MOTIFS transfac ".${TMP}."/".$output_path."/".$output_prefix."ref_motifs.tf";
  $parameters .= " -ref_motifs ".$refmotif_file;
  push(@tasks, "ref_motifs,motifs_vs_ref");
  #  } else {
  #    &FatalError ("If you want to upload a personal matrix file, you should specify the location of this file on your hard drive with the Browse button");
  #  }
}

################################################################
## Scan sequences to search motif occurrences (matrix-scan-quick).
if ($query->param('matrix-scan-quick') =~ /on/) {
    push(@tasks, "scan");
    ## HERE need to add the pval and markov order for the background model for matrix-scan-quick

}

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

### add -task
$parameters .= " -task " . join(",", @tasks);

### output prefix
$parameters .= " -prefix $output_prefix";


### other default parmaters
$parameters .= " -2str -noov -img_format png ";

### output directory
$parameters .= " -outdir $output_path";

################################################################
## display or send result
$index_file = $output_dir."/".$output_prefix."_synthesis.html";
my $mail_title = join (" ", "[RSAT]", "peak-motifs", &AlphaDate());
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere", $index_file, title=>$mail_title ,no_email=>1);
} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);
}

################################################################
## result page footer
print $query->end_html;

exit(0);

