#!/usr/bin/perl

#### redirect error log to a file
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "/tmp/RSA_ERROR_LOG.txt";
#    $ERR_LOG = "/rubens/dsk2/jvanheld/rsa/rsa-tools/logs/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
@result_files = ();

#### TEMPORARY
$oligo_analysis_command = $SCRIPTS."/oligo-analysis";
$convert_seq_command = $SCRIPTS."/convert-seq";
$purge_sequence_command = $SCRIPTS."/purge-sequence";
$tmp_file_name = sprintf "oligo-analysis.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("oligo-analysis result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();

#### read parameters ####
$parameters = "";
$parameters .= " -sort";

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("input sequence",$sequence_file);

$purge = $query->param('purge');
if ($purge) {
    #### purge sequence option
#    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format |  $oligo_analysis_command ";
    $command = "$purge_sequence_command -i $sequence_file -format $sequence_format -o ${sequence_file}.purged ";
    $command .= " -dna" if (lc($query->param('sequence_type')) eq "dna");
    $command .= "; $oligo_analysis_command -i ${sequence_file}.purged -format fasta ";
} else {
    $command = "$oligo_analysis_command -i $sequence_file  ";
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
if ($query->param('freq_estimate') =~ /background/i) {

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

} elsif ($query->param('freq_estimate') =~ /upload/i) {
  ## User-specific background model
  $exp_freq_file = "${TMP}/$tmp_file_name.expfreq";
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

} elsif ($query->param('freq_estimate') =~ /residue frequenc/i) {
  ## Bernoulli background model
  $freq_option = " -bg bernoulli";

} elsif ($query->param('freq_estimate') =~ /markov/i) {
  ## Markov model calibrated on input sequences
  $freq_option = " -markov";
  if (&IsNatural($query->param('markov_order'))) {
    $freq_option .= " ".$query->param('markov_order');
  }

} elsif ($query->param('freq_estimate') =~ /lexicon/i) {
  ## Lexicon background model
  $freq_option = " -lexicon";

} elsif ($query->param('freq_estimate') =~ /equiprobable/i) {
  ## Equiprobable residues
  $freq_option = " -bg equi";

} else {
  $freq_option = "";
} 
$parameters .= "$freq_option";

#### pseudo weight
if (&IsReal($query->param('pseudo_weight'))) {
    my $pseudo = $query->param('pseudo_weight');
    $parameters .= " -pseudo $pseudo";
}

#### neighborhood ####
if ($query->param('neighborhood') =~ /N at one position/i) {
  $parameters .= " -oneN";
} elsif ($query->param('neighborhood') =~ /one degenerated position/i) {
  $parameters .= " -onedeg ";
}

$command .= $parameters;

print "<PRE>command: $command<P>\n</PRE>" if ($ENV{rsat_echo} >=1);

&SaveCommand("$command", "$TMP/$tmp_file_name");

if ($query->param('output') =~ /display/i) {
    &PipingWarning();

    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$command |";

    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, 1);
    close(RESULT);
    push @result_files, ('oligos', $result_file);

    #### oligonucleotide assembly ####
    if (($query->param('return') ne "table") &&
	($query->param('return') ne "distrib") &&
	(&IsReal($query->param('lth_occ_sig')))) {

      ## Pattern-assembly
      $assembly_file = "$TMP/$tmp_file_name.asmb";
      $top_patterns = 50;
      $pattern_assembly_command = $SCRIPTS."/pattern-assembly -v 1 -subst 1 -top ".$top_patterns;
      if ($query->param('strand') =~ /single/) {
	$pattern_assembly_command .= " -1str";
      } else {
	$pattern_assembly_command .= " -2str";
      }
      $pattern_assembly_command .= "  -i $result_file";
      $pattern_assembly_command .= "  -o $assembly_file";

      unless ($ENV{RSA_ERROR}) {

	## Assemble the significant patterns
	print "<H2>Pattern assembly</H2>\n";
	print "<PRE>pattern-assembly command: $pattern_assembly_command<P>\n</PRE>" if ($ENV{rsat_echo} >=1);
	system "$pattern_assembly_command";
	open ASSEMBLY, $assembly_file;
	print "<PRE>\n";
	while (<ASSEMBLY>) {
	  s|$ENV{RSAT}/||g;
	  print;
	}
	print "</PRE>\n";
	close(ASSEMBLY);
	push @result_files, ('assembly', $assembly_file);


	## Convert pattern-assembly result into PSSM
	if ($query->param('to_matrix')) {
	  &MatrixFromPatterns_run();
	}
      }
    }

    &PrintURLTable(@result_files);
    &OligoDyadPipingForm();
    print '<HR SIZE=3>';

} else {
    &EmailTheResult("$command", $query->param('user_email'), $tmp_file_name);
}

print $query->end_html;

exit(0);





