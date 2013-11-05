#!/usr/bin/perl

## Get the path for RSAT Perl libraries
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}

#### redirect error log to a file
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "/tmp/RSA_ERROR_LOG.txt";
#    $ERR_LOG = "/rubens/dsk2/jvanheld/rsa/rsa-tools/logs/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}

## Load RSAT librarires
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

## Open result page
&RSA_header("oligo-analysis result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();

## Commands
$oligo_analysis_command = $SCRIPTS."/oligo-analysis";
$convert_seq_command = $SCRIPTS."/convert-seq";
$purge_sequence_command = $SCRIPTS."/purge-sequence";
$prefix = "oligo-analysis";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

#$tmp_file_name = sprintf "oligo-analysis.%s", &AlphaDate();

#### read parameters ####
$parameters = "";
$parameters .= " -sort";

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("Input sequence", $sequence_file);

## Sequence purging seems not to work with proteins -> only apply it
## for DNA sequences
my $sequence_type = lc($query->param('sequence_type'));
if ($sequence_type eq "dna") {
  $purge = $query->param('purge');
}

## Purge sequence if required, and start oligo-analysis command
if ($purge) {
  $purged_seq_file = $sequence_file.".purged";
  push @result_files, ("Purged sequence", $purged_seq_file);

  #### purge sequence option
  #    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format |  $oligo_analysis_command ";
  $command = $purge_sequence_command;
  $command .= " -i ".$sequence_file;
  $command .= " -format ".$sequence_format;
  $command .= " -o ".$purged_seq_file;
  $command .= " -seqtype ".$sequence_type if ($sequence_type eq "dna");
  $command .= "; ".$oligo_analysis_command." -i ".$purged_seq_file." -format fasta ";
} else {
  $command = $oligo_analysis_command." -i ".$sequence_file;
}

### fields to return
if ($query->param('return') eq "table") {
    $parameters .= " -return occ -table"; 
} elsif ($query->param('return') eq "distrib") {
    $parameters .= " -return occ -distrib"; 
} else {
  &CGI_return_fields();
}

### single or both strands
$str = "";
if ($query->param('strand') =~ /single/) {
  $str = " -1str";
  $parameters .= " -1str";
} else {
  $str = " -2str";
  $parameters .= " -2str";
}

### group patterns by pairs of reverse complements
unless ($query->param('grouprc')) {
  $parameters .= " -nogrouprc";
}

### group patterns by pairs of reverse complements
if ($query->param('side') eq 'under-represented') {
  $parameters .= " -under";
} elsif ($query->param('side') eq 'both') {
  $parameters .= " -two_tail";
}

### prevent overlapping matches of the same pattern
if ($query->param('noov')) {
  $parameters .= " -noov";
}

### verbose
$parameters .= " -v 1 ";

### quick count mode
$parameters .= " -quick_if_possible ";

#### sequence type
$parameters .= " -seqtype ".$query->param("sequence_type");

#### oligo size ####
$oligo_length = $query->param('oligo_length') ;
&FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
$parameters .= " -l $oligo_length";

################################################################
## Background model
if ($query->param('bg_method') =~ /background/i) {

  ## Check background subset
  %supported_background = (
			   "upstream"=>1,
			   "upstream-noorf"=>1,
#			   "intergenic"=>1,
			   "protein"=>1,
			  );
  my $background = $query->param("background");
  unless ($supported_background{$background}) {
    &cgiError("$background is not supported as background model");
  }
  $freq_option = " -bg $background";

  ### Organism-specific background model
  if ($query->param('bg_level') eq 'organism') {
    unless ($organism = $query->param('organism')) {
      &cgiError("You should specify an organism to use intergenic frequency calibration");
    }
    unless (%{$supported_organism{$organism}}) {
      &cgiError("Organism $org is not supported on this site");
    }
    $freq_option .= " -org $organism";
  }

  ### Taxon-specific background model
  if ($query->param('bg_level') eq 'taxon') {
    unless ($taxon = $query->param('taxon')) {
      &cgiError("You should specify an taxon to use intergenic frequency calibration");
    }
    &CheckTaxon($taxon);
    $freq_option .= " -taxon $taxon";
  }

} elsif ($query->param('bg_method') eq 'freq_file_upload') {
  ## User-specific expected freqency file (oligos of same sise as analyzed oligos)
  $exp_freq_file = $tmp_file_path.".expfreq";
#  $exp_freq_file = "${TMP}/$tmp_file_name.expfreq";
  push @result_files, ('Expected frequencies', $exp_freq_file);
  $upload_freq_file = $query->param('upload_freq_file');
  if ($upload_freq_file) {
    ## Support compressed .gz files
    if ($upload_freq_file =~ /\.gz$/) {
      $exp_freq_file .= ".gz";
    }
    $type = $query->uploadInfo($upload_freq_file)->{'Content-Type'};
    open FREQ, ">$exp_freq_file" ||
      &cgiError("Cannot store expected frequency file in temp dir.");
    while (<$upload_freq_file>) {
      print FREQ;
    }
    close FREQ;
    $freq_option = " -expfreq $exp_freq_file";
  } else {
    &FatalError ("If you want to upload an expected frequency file, you should specify the location of this file on your hard drive with the Browse button");
  }

} elsif ($query->param('bg_method') eq 'file_upload') {

  ## User-specific background model (any Markov order)
  $bgfile = $tmp_file_path.".bgfile";
  push @result_files, ('Background file', $bgfile);
#  $bgfile = "${TMP}/$tmp_file_name.bgfile";
  $upload_bgfile = $query->param('upload_bgfile');
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
    $freq_option = " -bgfile ".$bgfile;

    ## Background format
    $freq_option .= " -bg_format ".$query->param('bg_format');

  } else {
    &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
  }

} elsif ($query->param('bg_method') =~ /residue frequenc/i) {
  ## Bernoulli background model
  $freq_option = " -bg bernoulli";

} elsif ($query->param('bg_method') =~ /markov/i) {
  ## Markov model calibrated on input sequences
  $freq_option = " -markov";
  if (&IsNatural($query->param('markov_order'))) {
    $freq_option .= " ".$query->param('markov_order');
  }

} elsif ($query->param('bg_method') =~ /lexicon/i) {
  ## Lexicon background model
  $freq_option = " -lexicon";

} elsif ($query->param('bg_method') =~ /equiprobable/i) {
  ## Equiprobable residues
  $freq_option = " -bg equi";

} else {
  $freq_option = "";
}
$parameters .= "$freq_option";

#### pseudo weight
if (&IsReal($query->param('pseudo_freq'))) {
    my $pseudo = $query->param('pseudo_freq');
    $parameters .= " -pseudo $pseudo";
}

#### neighborhood ####
if ($query->param('neighborhood') =~ /N at one position/i) {
  $parameters .= " -oneN";
} elsif ($query->param('neighborhood') =~ /one degenerated position/i) {
  $parameters .= " -onedeg ";
}

$command .= $parameters;

&ReportWebCommand($command." ".$parameters);

&SaveCommand("$command", $tmp_file_path);

if ($query->param('output') =~ /display/i) {
    &PipingWarning();

    ### execute the command ###
    $result_file = $tmp_file_path.".tab";
#    $result_file = "$TMP/$tmp_file_name.res";
    push @result_files, ('oligos', $result_file);

    open RESULT, "$command |";

    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);

    #### oligonucleotide assembly ####
    if (($query->param('return') ne "table") &&
	($query->param('return') ne "distrib") &&
	(&IsReal($query->param('lth_occ_sig')))) {


      ## Assemble the significant patterns with pattern-assembly
      $assembly_file = $tmp_file_path.".asmb";
      push @result_files, ('Assembly', $assembly_file);
      $pattern_assembly_command = $SCRIPTS."/pattern-assembly -v 1 -subst 1 -toppat 50";
      if ($query->param('strand') =~ /single/) {
	$pattern_assembly_command .= " -1str";
      } else {
	$pattern_assembly_command .= " -2str";
      }
      if (&IsNatural($query->param('max_asmb_nb'))) {
	$pattern_assembly_command .= " -max_asmb_nb ".$query->param('max_asmb_nb');
      }
      $pattern_assembly_command .= " -i ".$result_file;
      $pattern_assembly_command .= " -o ".$assembly_file;

      unless ($ENV{RSA_ERROR}) {
	## Assemble the significant patterns
	print "<H2>Pattern assembly</H2>\n";
	print "<PRE>pattern-assembly command: ", &RSAT::util::hide_RSAT_path($pattern_assembly_command), "<P>\n</PRE>" if ($ENV{rsat_echo} >=1);
	system "$pattern_assembly_command";
	open ASSEMBLY, $assembly_file;
	print "<PRE>\n";
	while (<ASSEMBLY>) {
#	  s|$ENV{RSAT}/||g;
	  print;
	}
	print "</PRE>\n";
	close(ASSEMBLY);

	## Convert pattern-assembly result into PSSM
	if ($query->param('to_matrix')) {
	  if ($sequence_type eq "dna") {
	    &MatrixFromPatterns_run();
	  } else {
	    &RSAT::message::Warning("Conversion to matrix is only supported for DNA sequences");
	  }
	}
      }
    }

    &PrintURLTable(@result_files);
    &OligoDyadPipingForm();
    print '<HR SIZE=3>';

} else {
    &EmailTheResult("$command", $query->param('user_email'), $tmp_file_path);
}

print $query->end_html;

exit(0);





