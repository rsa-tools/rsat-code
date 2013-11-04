#!/usr/bin/perl

## CVS
## added the possibility to specify the expected frequency for each nucleotide separately

if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
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
&RSA_header("Random sequence result", "results");


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/random-seq";
$prefix = "random-seq";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
#$tmp_file_name = sprintf "random-seq.%s", &AlphaDate();
@result_files = ();

$size_limit = 25e+6;

################################################################
#### read parameters ####
$parameters = "";

## template file (optional)
($template_file, $template_format) = &MultiGetSequenceFile(1, "$TMP/$tmp_file_name"."_template.fa", 0);

## a template file has been given
if ($template_file) {
  push @result_files, ("Template file ($template_format)",$template_file);

  ## Compute sequence lengths from the template sequence file
  my $length_file = "$TMP/$tmp_file_name".".lengths";
  push @result_files, ("Sequence lengths",$length_file);

  my $seqlength_cmd = $SCRIPTS."/sequence-lengths -v 1 -i ".$template_file;
  $seqlength_cmd .= " -in_format ".$template_format;
  $seqlength_cmd .= " -o ".$length_file;
  system($seqlength_cmd);

  ## Add the sequence length file as template for random-genome-fragments
  $parameters .= " -template_format len -i ".$length_file;
#  $parameters .= " -template_format fasta -i ".$template_file;

#  $parameters .= " -template_format fasta -i ".$template_file;
} else {
  #### sequence length
  $length = $query->param('length');
  if (&IsNatural($length)) {
    $parameters .= " -l $length ";
  } else {
    &FatalError("Sequence length must be a natural number");
  }


  #### number of repetitions
  $repet = $query->param('repet');
  if (&IsNatural($repet)) {
    $parameters .= " -n $repet";
  } else {
    &FatalError("Repetitions must be a natural number");
  }
  #### check lengths and repetitions
  if ($repet * $length == 0) {
    &FatalError("Sequence length and reperitions must be non-null");
  }
  if ($repet*$length > $size_limit) {
    &FatalError("The web interface does not support queries of this size. Maximum size per query (length * repetitions) = ".$size_limit);
  }
}


#### line width
if (&IsNatural($query->param('lw'))) {
    $parameters .= " -lw ".$query->param('lw');
}

#### output format
$out_format = $query->param('format');
&CheckOutputSeqFormat($out_format);
$parameters .= " -format $out_format ";

#### alphabet
if ($query->param('bg_method') eq "alphabet") {

    $freq{A} = $query->param('Afreq');
    $freq{C} = $query->param('Cfreq');
    $freq{T} = $query->param('Tfreq');
    $freq{G} = $query->param('Gfreq');

    ## Check the values 
    foreach my $letter (keys %freq) {
	unless (&IsReal($freq{$letter})) {
	    &FatalError("Invalid frequency value ".$freq{$letter}." for residue ".$letter);
	}
    }

    ## Print residue frequencies in a file
    $alphabet_file = $tmp_file_path.".alphabet";
    push @result_files, ("residue priors",$alphabet_file);
    open ALPHA, ">".$alphabet_file;
    foreach my $letter (keys %freq) {
	print ALPHA $letter, "\t", $freq{$letter}, "\n";
    }
    close ALPHA;
    $parameters .= " -expfreq ".$alphabet_file;
    &DelayedRemoval($alphabet_file);

## Pre-calibrated Markov models
} elsif (($query->param('bg_method') =~ /upstream/i) ||
	 ($query->param('bg_method') =~ /protein/i)) {
    ### check organism
    unless ($organism = $query->param('organism')) {
	&cgiError("You should specify an organism to use upstream frequency calibration");
    }
    unless (defined(%{$supported_organism{$organism}})) {
	&cgiError("Organism $organism is not supported on this site");
    }
    if ($query->param('bg_method') =~ /protein/i) {
      $oligopept_size = $query->param("oligopept_size");
      unless (&IsNatural($oligopept_size)) {
	&cgiError("Invalid oligopeptide length $oligopept_size");
      }
      $seq_type = "protein"; ## Used for the piping form
      $parameters .= " -bg protein -org $organism -ol $oligopept_size -type protein";
    } else {
      $oligo_size = $query->param("oligo_size");
      unless (&IsNatural($oligo_size)) {
	&cgiError("Invalid oligonucleotide length $oligo_size");
      }
      $seq_type = "dna"; ## Used for the piping form
      $parameters .= " -bg upstream-noorf -org $organism -ol $oligo_size -type dna";
    }
  } elsif ($query->param('bg_method') eq 'file_upload') {

  ## User-specific background model (any Markov order)
  my $bgfile = $tmp_file_path.".bgfile";
  push @result_files, ('Background file', $bgfile);
  my $upload_bgfile = $query->param('upload_bgfile');
  if ($upload_bgfile) {
    ## Support compressed .gz files
    if ($upload_bgfile =~ /\.gz$/) {
      $bgfile .= ".gz";
    }
    $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
    open BG, ">$bgfile" ||
      &cgiError("Cannot store background model file in temp dir.");
    while (<$upload_bgfile>) {
      print BG;
    }
    close BG;
    $parameters .= " -expfreq ".$bgfile;

  } else {
    &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
  }
  }
  

&ReportWebCommand($command." ".$parameters);

## Output file
$sequence_file = $tmp_file_path.".".$out_format;
push @result_files, ("sequence",$sequence_file);

### execute the command ###
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
#if ($query->param('output') eq "display") {

    open RESULT, "$command $parameters |";

    ### print the result ###
    &PipingWarning();
    print '<H4>Result</H4>';

    ### open the mirror file ###
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

    ### prepare data for piping
    &PrintURLTable(@result_files);
    &PipingFormForSequence();
    print "<HR SIZE = 3>";

    &DelayedRemoval($sequence_file);


} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $sequence_file);
}
print $query->end_html;

exit(0);

