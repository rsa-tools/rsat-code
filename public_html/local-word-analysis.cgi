#!/usr/bin/perl

## Get the path for RSAT Perl libraries
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;

#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}

## Load required libraries
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("local-word-analysis result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

@result_files = ();
## Commands to be used
$motif_command = "$ENV{RSAT}/python-scripts/local-word-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$prefix = "local-word-analysis";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

#### read parameters ####
$parameters = "";
#$parameters .= " -sort";

$parameters .= " --sort=-occ_sig";
#$parameters .= " --count=hash ";
#$parameters .= " --count=tree --spacing=1:10";

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("Input sequence", $sequence_file);


### window
if ($query->param('windowtype') =~ /no/ ){
    $parameters .= ' --window=none';
}elsif ($query->param('windowtype') =~ /fixed/){
  if (IsReal($query->param('window_width'))) {
    if ($query->param('window_group')) {
      $parameters .= ' --center=0 --windowgroup=' . $query->param('window_width');                
    }else{
      $parameters .= ' --window=' . $query->param('window_width');        
    }
  }

}else {
  $parameters .= ' --heuristic=slices';
  #variable size
}

if (IsReal($query->param('bg_window_width'))) {
  $parameters .= ' --bgwindow=' . $query->param('bg_window_width');
}

### filters
if (IsReal($query->param('lth_occ'))) {
    $parameters .= " --min occ " . $query->param('lth_occ');
}
if (IsReal($query->param('uth_occ'))) {
    $parameters .= " --max occ " . $query->param('uth_occ');
}
if (IsReal($query->param('lth_occ_E'))) {
    $parameters .= " --min occ_E " . $query->param('lth_occ_E');
}
if (IsReal($query->param('uth_occ_E'))) {
    $parameters .= " --max occ_E " . $query->param('uth_occ_E');
}
if (IsReal($query->param('lth_occ_sig'))) {
    $parameters .= " --min occ_sig " . $query->param('lth_occ_sig');
}
if (IsReal($query->param('uth_occ_sig'))) {
    $parameters .= " --max occ_sig " . $query->param('uth_occ_sig');
}
if (IsReal($query->param('lth_rank'))) {
    $parameters .= " --min rank " . $query->param('lth_rank');
}
if (IsReal($query->param('uth_rank'))) {
    $parameters .= " --max rank " . $query->param('uth_rank');
}

if (IsReal($query->param('lth_w_rank'))) {
    $parameters .= " --min w_rank " . $query->param('lth_w_rank');
}
if (IsReal($query->param('uth_w_rank'))) {
    $parameters .= " --max w_rank " . $query->param('uth_w_rank');
}


### align
if ($query->param('align') =~ /right/) {
  $parameters .= " --right=-1";
}else{
  $parameters .= " --left=1";    
}

    
### single or both strands
$str = "";
if ($query->param('strand') =~ /single/) {
  $str = " -1str";
  $strand = "-1str";
  $parameters .= " --strand=+";
} else {
  $str = " -2str";
  $strand = "-2str";
  $parameters .= " --strand=+-";
}

### group patterns by pairs of reverse complements
#unless ($query->param('grouprc')) {
#  $parameters .= " -nogrouprc";
#} 

### group patterns by pairs of reverse complements
#if ($query->param('side') eq 'under-represented') {
#  $parameters .= " -under";
#} elsif ($query->param('side') eq 'both') {
#  $parameters .= " -two_tail";
#} 

### prevent overlapping matches of the same pattern
if (! $query->param('noov')) {
  $overlap='-ovlp';
  $parameters .= " --overlap";
}else{
  $overlap='-noov';
}

### verbose
$parameters .= " -v 5";

#### sequence type
#$parameters .= " -seqtype ".$query->param("sequence_type");

#### oligo type and size ####

if ($query->param('oligotype') =~ /dyad/i) {
    $oligotype="dyad";
    $oligo_length = $query->param('monad_length');
    $parameters .= " --dyad";
    $parameters .= " -l $oligo_length";
    $spacing_a = $query->param('spacing_a');
    $spacing_b = $query->param('spacing_b');
    &FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
    $parameters .= " --spacing=$spacing_a:$spacing_b";
}else{
    $oligotype="oligo";
    $oligo_length = $query->param('oligo_length');
    &FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
    $parameters .= " -l $oligo_length";
}


#### expected frequency estimation ####
if ($query->param('bg_method') =~ /background/i) {
    %supported_background = (
			      "upstream"=>1,
			      "upstream-noorf"=>1,
#			      "intergenic"=>1
			      );

	### check organism
	unless ($organism = $query->param('organism')) {
	    &cgiError("You should specify an organism to use intergenic frequency calibration");
	}
	unless (%{$supported_organism{$organism}}) {
	    &cgiError("Organism $org is not supported on this site");
	}
	my $background = $query->param("background");
	unless ($supported_background{$background}) {
	    &cgiError("$background is not supported as background model");
	}
    ### Taxon-specific background model
    if ($query->param('bg_level') eq 'taxon') {
      unless ($taxon = $query->param('taxon')) {
        &cgiError("You should specify an taxon to use intergenic frequency calibration");
      }
      &CheckTaxon($taxon);
      $organism = $taxon;
    }

    #$exp_freq_file = "$ENV{RSAT}/public_html/data/genomes/$organism/oligo-frequencies/" . "$oligo_length" . "nt_" . "$background" . "_" . "$organism$overlap$strand.freq.gz";
    $exp_freq_file = &ExpectedFreqFile($organism, $oligo_length, $background, type=>$oligotype, noov=>$overlap, str=>$strand, taxon=>$taxon);

    #print $exp_freq_file;
	#$freq_option = " -bg $background -org $organism";
	$freq_option = " --bgoligo=$exp_freq_file.gz";


  } elsif ($query->param('freq_estimate') =~ /upload/i) {
    $exp_freq_file = $tmp_file_path.".expfreq";
    push @result_files, "Expected frequencies", $exp_freq_file;
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
      $freq_option = " --bgoligo=$exp_freq_file";
    } else {
      &FatalError ("If you want to upload an expected frequency file, you should specify the location of this file on your hard drive with the Browse button");
    }

} elsif ($query->param('freq_estimate') =~ /Equiprobable residues/i) {
  $freq_option = " --markov=0";
} elsif ($query->param('freq_estimate') =~ /markov/i) {
  if (&IsNatural($query->param('markov_order'))) {
      $freq_option .= " --markov=".$query->param('markov_order');
  }
} else {
  $freq_option = "";
} 
$parameters .= "$freq_option";


## command

$purge = $query->param('purge');
if ($purge) {
  $purged_seq_file = $sequence_file.".purged";
  push @result_files, ("Purged sequence",$purged_seq_file);

  #### purge sequence option
  #    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format |  $oligo_analysis_command ";
  $command = $purge_sequence_command;
  $command .= " -i ".$sequence_file;
  $command .= " -format ".$sequence_format;
  $command .= " -o ".$purged_seq_file;
  $command .= " -seqtype ".$sequence_type if ($sequence_type eq "dna");
    #### purge sequence option
#    $command = "$purge_sequence_command -i $sequence_file -format $sequence_format -o ${sequence_file}.purged; "
  $command .=  "; $motif_command -i $purged_seq_file $parameters";

} else {
  $command .=  "; $motif_command -i $sequence_file $parameters";
}

#print '<style> <!-- pre {overflow: auto;} --></style>';

## Output file
$result_file = $tmp_file_path.".tab";
push @result_files, "Result file (tab)", $result_file;

&ReportWebCommand($command);

#&SaveCommand("$command", &RSAT::util::get_pub_temp()."/$tmp_file_name");

if ($query->param('output') =~ /display/i) {

    &PipingWarning();
    
    ### execute the command ###
    open RESULT, "$command  |";

    
    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);
    
    #### oligonucleotide assembly ####
    if (($query->param('return') ne "table") &&
	($query->param('return') ne "distrib") &&
	(&IsReal($query->param('lth_occ_sig')))) {

      ## Assemble the significant patterns with pattern-assembly
      $assembly_file = $tmp_file_path.".asmb";
      push @result_files, "Assembly", $assembly_file;
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

#       ## Pattern-assembly
#       $assembly_file = &RSAT::util::get_pub_temp()."/$tmp_file_name.asmb";
#       $pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 1 -toppat 50";
#       if ($query->param('strand') =~ /single/) {
# 	$pattern_assembly_command .= " -1str";
#       } else {
# 	$pattern_assembly_command .= " -2str";
#       }
#       $pattern_assembly_command .= "  -i $result_file";
#       $pattern_assembly_command .= "  -o $assembly_file";

      unless ($ENV{RSA_ERROR}) {

	## Assemble the significant patterns
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
      }
    }
    &PrintURLTable(@result_files);
    &PipingForm();
    print '<HR SIZE=3>';

} else {
    &EmailTheResult("$command", $query->param('user_email'), $result_file);
}

print $query->end_html;

exit(0);


sub PipingForm {
    ### prepare data for piping

    #### title
    $title = $query->param('title');
    $title =~ s/\"/\'/g;

    #### strand for pattern-assembly
    if ($query->param('strand') =~ /single/) {
	$strand_opt .= " sensitive";
    } else {
	$strand_opt .= " insensitive";
    }
  print <<End_of_form;
<HR SIZE = 3>
<TABLE CLASS = "nextstep" CELLSPACING=0 CELLPADDING=10 BORDER=0 NOWRAP>
<TR>

<TR VALIGN="top" ALIGN="center">
    <Th VALIGN=BOTTOM ALIGN=CENTER COLSPAN=6>
	Next step
    </Th>

</TR>

<td valign=bottom align=center>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="string-based pattern matching (dna-pattern)">
</FORM>
</TD>

<!--
<td valign=bottom align=center>
<b><font color="red">New !</font></b>
<FORM METHOD="POST" ACTION="matrix-scan_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="matrix_file" VALUE="$pssm_file">
<INPUT type="hidden" NAME="matrix_format" VALUE="tab">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="matrix-based pattern matching (matrix-scan)">
</FORM>
</TD>
-->


<td valign=bottom align=center>
<FORM METHOD="POST" ACTION="pattern-assembly_form.cgi">
<INPUT type="hidden" NAME="local_pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="subst" VALUE=1>
<INPUT type="hidden" NAME="maxfl" VALUE=1>
<INPUT type="hidden" NAME="sc" VALUE="auto">
<INPUT type="hidden" NAME="strand" VALUE=$strand_opt>
<INPUT type="submit" value="pattern assembly">
</FORM>
</TD>


</TR>
</TABLE>
End_of_form

}



