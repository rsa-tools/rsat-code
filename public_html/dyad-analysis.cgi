#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis.cgi,v 1.54 2013/11/03 19:33:31 jvanheld Exp $
#
# Time-stamp: <2003-10-11 00:30:17 jvanheld>
#
############################################################

## Get the path for RSAT Perl libraries
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

## Load RSAT libraries
require "RSA.lib";
require "RSA2.cgi.lib";
require "RSA.disco.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("dyad-analysis result", "results");
&ListParameters if ($ENV{rsat_echo} >= 2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();


## Commands to be used
$dyad_analysis_command = "$SCRIPTS/dyad-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$prefix = "dyad-analysis";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "dyad-analysis.%s", &AlphaDate;

@result_files = ();

#### read parameters ####
$parameters = " -v 1 -quick -sort";

################################################################
## Timeout prevents over-utilization of the server when large files are analyzed
## After 1h of computation, the dyad-anlysis process is killed 
$parameters .= " -timeout 3600 ";

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("Input sequence",$sequence_file);

$purge = $query->param('purge');
if ($purge) {
  #### purge sequence option
  #    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format | $dyad_analysis_command ";
  $purged_seq_file = $sequence_file.".purged";
  push @result_files, ("Purged sequence", $purged_seq_file);
  $command = $purge_sequence_command;
  $command .= " -i ".$sequence_file;
  $command .= " -format ".$sequence_format;
  $command .= " -o ".$purged_seq_file;
  $command .= " -seqtype ".$sequence_type if ($sequence_type eq "dna");
  $command .= "; $dyad_analysis_command -i ".$purged_seq_file;
} else {
  $command= "$dyad_analysis_command -i $sequence_file  ";
}


### sequence file
#($sequence_file,$sequence_format) = &GetSequenceFile();
#$parameters .= " -i $sequence_file -format $sequence_format";

### dyad type
if ($query->param('dyad_type') =~ /invert/) {
  $parameters .= " -type ir";
} elsif ($query->param('dyad_type') =~ /direct/) {
  $parameters .= " -type dr";
} elsif ($query->param('dyad_type') =~ /repeat/) {
  $parameters .= " -type rep";
} else {
  $parameters .= " -type any";
}

### single strand search
if ($query->param('strand') =~ /single/) {
  $parameters .= " -1str";
} else {
  $parameters .= " -2str";
}

### group patterns by pairs of reverse complements
if ($query->param('side') eq 'under-represented') {
  $parameters .= " -under";
} elsif ($query->param('side') eq 'both') {
  $parameters .= " -two_tail";
}

### overlapping matches
if ($query->param('noov')) {
  $parameters .= " -noov";
}

################################################################
## Treat the return fields and thresholds
&CGI_return_fields();

#### oligo size ####
if (&IsNatural($query->param('oligo_size'))) {
  $oligo_length = $query->param('oligo_size') ;
} 
$parameters .= " -l $oligo_length";

#### spacing ####
unless ((&IsNatural($query->param('spacing_from'))) && (&IsNatural($query->param('spacing_from'))))  {
  &cgiError("Invalid spacing specification.");    
}
if ($query->param('spacing_from') >$query->param('spacing_to')) {
  $spacing = $query->param('spacing_to')."-".$query->param('spacing_from') ;   
} else {
  $spacing = $query->param('spacing_from')."-".$query->param('spacing_to') ;   
}
$parameters .= " -spacing $spacing";

################################################################
## Background model
if ($query->param('bg_method') eq 'background') {
  ## Check genome subset
  %supported_background = (
			   "upstream"=>1,
			   "upstream-noorf"=>1,
#			   "intergenic"=>1
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
  $parameters .= " ".$freq_option;

  ### organism
  #if ($organism = $query->param('organism')) {
  #  $parameters .= " -org $organism";
  #}
  #    $parameters .= " -bg ".$query->param('background');
} elsif ($query->param('bg_method') =~ /upload/i) {
  $exp_freq_file = $tmp_file_path.".expfreq";
#  $exp_freq_file = "${TMP}/$tmp_file_name.expfreq";
  push @result_files, ('exp_freq', $exp_freq_file);
  $upload_freq_file = $query->param('upload_freq_file');
  if ($upload_freq_file) {
    if ($upload_file =~ /\.gz$/) {
      $exp_freq_file .= ".gz";
    }
    $type = $query->uploadInfo($upload_freq_file)->{'Content-Type'};
    open FREQ, ">$exp_freq_file" ||
      &cgiError("Cannot store expected frequency file in temp dir.");
    while (<$upload_freq_file>) {
      print FREQ;
    }
    close FREQ;
    $parameters .= " -expfreq $exp_freq_file";
  } else {
    &FatalError ("If you want to upload an expected frequency file, you should specify the location of this file on your hard drive with the Browse button");
  }

} else {
  unless ($query->param('bg_method') eq 'monads') {
    &FatalError("Invalid expected frequency calibration");
  }
}

$command .= $parameters;

## Output file
$result_file = $tmp_file_path.".tab";
push @result_files, ('dyads', $result_file);

&ReportWebCommand($command);

&SaveCommand("$command", $tmp_file_path);

if ($query->param('output') eq "display") {  

  &PipingWarning();

  ### execute the command ###
  open RESULT, "$command | ";

  ### Print result on the web page
  print '<H4>Result</H4>';
  &PrintHtmlTable(RESULT, $result_file, 1);
  close(RESULT);

  #### pattern assembly ####
  if ((&IsReal($query->param('lth_occ_sig'))) && ($query->param('lth_occ_sig')>= -1)) {

    ## Assemble the significant patterns with pattern-assembly
    $assembly_file = $tmp_file_path.".asmb";
#    $assembly_file = "$TMP/$tmp_file_name.asmb";
    push @result_files, ('assembly', $assembly_file);
    $pattern_assembly_command = $SCRIPTS."/pattern-assembly -v 1 -subst 0 -toppat 50";
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

    print "<H2>Pattern assembly</H2>\n";
    print "<PRE>pattern-assembly command: ", &RSAT::util::hide_RSAT_path($pattern_assembly_command), "<P>\n</PRE>" if ($ENV{rsat_echo} >=1);
    system "$pattern_assembly_command";
    open ASSEMBLY, $assembly_file;
    print "<PRE>\n";
    while (<ASSEMBLY>) {
      s|$ENV{RSAT}/||g;
      print;
    }
    print "</PRE>\n";
    close(ASSEMBLY);


    ## Convert pattern-assembly result into PSSM
    if ($query->param('to_matrix')) {
      &MatrixFromPatterns_run();
    }
  }

  &PrintURLTable(@result_files);
  &OligoDyadPipingForm();

} else {
  &EmailTheResult("$command", $query->param('user_email'), $result_file);
}



print $query->end_html;


exit(0);
