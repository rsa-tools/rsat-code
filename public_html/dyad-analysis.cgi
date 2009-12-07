#!/usr/bin/perl
############################################################
#
# $Id: dyad-analysis.cgi,v 1.34 2009/12/07 22:19:33 jvanheld Exp $
#
# Time-stamp: <2003-10-11 00:30:17 jvanheld>
#
############################################################
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
require "RSA.disco.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$dyad_analysis_command = "$SCRIPTS/dyad-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "dyad-analysis.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("dyad-analysis result", "results");
&ListParameters if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = " -v 1 -sort";

################################################################
## Timeout prevents over-utilization of the server when large files are analyzed
## After 1h of computation, the dyad-anlysis process is killed 
$parameters .= " -timeout 3600 ";

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
$purge = $query->param('purge');
if ($purge) {
  #### purge sequence option
  #    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format | $dyad_analysis_command ";
  $command= "$purge_sequence_command -i $sequence_file -format $sequence_format -o $sequence_file.purged ";
  $command .= " -dna" if (lc($query->param('sequence_type')) eq "dna");
  $command .= "; $dyad_analysis_command -i $sequence_file.purged ";
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
if ($query->param('freq_estimate') eq 'background') {
  ## Check genome subset
  %supported_background = (
			   "upstream"=>1,
			   "upstream-noorf"=>1,
			   "intergenic"=>1
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
    unless (defined(%{$supported_organism{$organism}})) {
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
} elsif ($query->param('freq_estimate') =~ /upload/i) {
  $exp_freq_file = "${TMP}/$tmp_file_name.expfreq";
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
  unless ($query->param('freq_estimate') eq 'monads') {
    &FatalError("Invalid expected frequency calibration");
  }
}

$command .= $parameters;

print "<PRE><B>Command:</B> $command </PRE>" if ($ENV{rsat_echo});

&SaveCommand($command, "$TMP/$tmp_file_name");

if ($query->param('output') eq "display") {  

  &PipingWarning();

  ### execute the command ###
  $result_file = "$TMP/$tmp_file_name.res";
  open RESULT, "$command | ";


  ### Print result on the web page
  print '<H4>Result</H4>';
  &PrintHtmlTable(RESULT, $result_file, true);
  close(RESULT);
  @result_files = ('oligos', "$tmp_file_name.res");

  #### pattern assembly ####
  if ((&IsReal($query->param('lth_occ_sig'))) && ($query->param('lth_occ_sig')>= -1)) {

    ## Assemble the significant patterns with pattern-assembly
    $assembly_file = "$TMP/$tmp_file_name.asmb";
    $pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 0 -top 50";
    if ($query->param('strand') =~ /single/) {
      $pattern_assembly_command .= " -1str";
    } else {
      $pattern_assembly_command .= " -2str";
    }
    $pattern_assembly_command .= "  -i $result_file";
    $pattern_assembly_command .= "  -o $assembly_file";

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
    push @result_files, ('assembly', "$tmp_file_name.asmb");


    ## Convert pattern-assembly result into PSSM
    if ($query->param('to_matrix')) {
      &MatrixFromPatterns_run();
    }

    # 	## Convert pattern-assembly result into PSSM 
    # 	if ($query->param('to_matrix')) {
    # 	  $pssm_prefix = $tmp_file_name."_pssm";
    # 	  $sig_matrix_file = $pssm_prefix."_sig_matrices.txt";
    # 	  $pssm_file = $pssm_prefix."_count_matrices.txt";
    # 	  $pssm_command = "$SCRIPTS/matrix-from-patterns -v 1 ".$str;
    # 	  $pssm_command .= " -seq ".$sequence_file;
    # 	  $pssm_command .= " -format $sequence_format";
    # 	  $pssm_command .= " -asmb ".$assembly_file;
    # 	  $pssm_command .= " -uth Pval 0.00025";
    # 	  $pssm_command .= " -bginput -markov 0";
    # 	  $pssm_command .= " -o ".$TMP."/".$pssm_prefix;
    # 	  print "<PRE>command to generate matrices (PSSM): $pssm_command<P>\n</PRE>" if ($ENV{rsat_echo} >=1);
    # 	  system "$pssm_command";
    # 	  push @result_files, ('significance matrices', $sig_matrix_file);
    # 	  push @result_files, ("gibbs matrices", $pssm_prefix."_gibbs_matrices.txt");
    # 	  push @result_files, ('count matrices', $pssm_file);
    # #	  print "<H2>Significance matrices</H2>\n";
    # #	  open SIG, $TMP."/".$sig_matrix_file;
    # #	  print "<PRE>\n";
    # #	  while (<SIG>) {
    # #	    s|$ENV{RSAT}/||g;
    # #	    print;
    # #	  }
    # #	  print "</PRE>\n";
    # #	  close(SIG);
    # 	  print "<H2>Matrices</H2>\n";
    # 	  open PSSM, $TMP."/".$pssm_file;
    # 	  print "<PRE>\n";
    # 	  while (<PSSM>) {
    # 	    s|$ENV{RSAT}/||g;
    # 	    print;
    # 	  }
    # 	  print "</PRE>\n";
    # 	  close(PSSM);
    # 	}
  }

  &PrintURLTable(@result_files);
  &OligoDyadPipingForm();

} elsif ($query->param('output') =~ /server/i) {
  &ServerOutput("$command", $query->param('user_email'));
} else {
  &EmailTheResult("$command", $query->param('user_email'));
}



print $query->end_html;


exit(0);
