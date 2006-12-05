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
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

#### TEMPORARY

$oligo_analysis_command = "$SCRIPTS/oligo-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "oligo-analysis.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("oligo-analysis result");
&ListParameters() if ($ECHO >=2);

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "";
$parameters .= " -sort";

### sequence file
($sequence_file,$sequence_format) = &GetSequenceFile();

$purge = $query->param('purge');
if ($purge) {
    #### purge sequence option
#    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format |  $oligo_analysis_command ";
    $command = "$purge_sequence_command -i $sequence_file -format $sequence_format -o ${sequence_file}.purged;  $oligo_analysis_command -i ${sequence_file}.purged -format fasta ";
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
    if ($query->param('strand') =~ /single/) {
	$parameters .= " -1str";
    } else {
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
$parameters .= " -v";

#### sequence type
$parameters .= " -seqtype ".$query->param("sequence_type");

#### oligo size ####
$oligo_length = $query->param('oligo_length') ;
&FatalError("$oligo_length Invalid oligonucleotide length") unless &IsNatural($oligo_length);
$parameters .= " -l $oligo_length";

#### expected frequency estimation ####
if ($query->param('freq_estimate') =~ /background/i) {
    %supported_background = (
			      "upstream"=>1,
			      "upstream-noorf"=>1,
			      "intergenic"=>1
			      );

#     if (($query->param("proba")) ||
# 	($query->param("ratio")) ||
# 	($query->param("zscore"))
# 	) {
	
	### check organism
	unless ($organism = $query->param('organism')) {
	    &cgiError("You should specify an organism to use intergenic frequency calibration");
	}
	unless (defined(%{$supported_organism{$organism}})) {
	    &cgiError("Organism $org is not supported on this site");
	}
	my $background = $query->param("background");
	unless ($supported_background{$background}) {
	    &cgiError("$background is not supported as background model");
	}

	$freq_option = " -bg $background -org $organism";
#    }

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
	$freq_option = " -expfreq $exp_freq_file";
    } else {
	&FatalError ("If you want to upload an expected frequency file, you should specify the location of this file on your hard drive with the Browse button");
    }

} elsif ($query->param('freq_estimate') =~ /residue frequenc/i) {
  $freq_option = " -bg bernoulli";
} elsif ($query->param('freq_estimate') =~ /markov/i) {
  $freq_option = " -markov";
  if (&IsNatural($query->param('markov_order'))) {
      $freq_option .= " ".$query->param('markov_order');
  }
} elsif ($query->param('freq_estimate') =~ /lexicon/i) {
  $freq_option = " -lexicon";
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

print "<PRE>command: $command<P>\n</PRE>" if ($ECHO >=1);

&SaveCommand("$command", "$TMP/$tmp_file_name");

if ($query->param('output') =~ /display/i) {

    &PipingWarning();
    
    ### execute the command ###
    $result_file = "$TMP/$tmp_file_name.res";
    open RESULT, "$command |";
    
    ### Print result on the web page
    print '<H2>Result</H2>';
    &PrintHtmlTable(RESULT, $result_file, true);
    close(RESULT);
    
    #### oligonucleotide assembly ####
    if (($query->param('return') ne "table") &&
	($query->param('return') ne "distrib") &&
	(&IsReal($query->param('lth_occ_sig')))) {

      ## Pattern-assembly
      $assembly_file = "$TMP/$tmp_file_name.asmb";
      $pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 1";
      if ($query->param('strand') =~ /single/) {
	$pattern_assembly_command .= " -1str";
      } else {
	$pattern_assembly_command .= " -2str";
      }
      $pattern_assembly_command .= "  -i $result_file";
      $pattern_assembly_command .= "  -o $assembly_file";
      
      unless ($ENV{RSA_ERROR}) {
	print "<H2>Pattern assembly</H2>\n";
	system "$pattern_assembly_command";
	open ASSEMBLY, $assembly_file;
	print "<PRE>\n";
	while (<ASSEMBLY>) {
	  s|$RSA/||g;
	  print;
	}
	print "</PRE>\n";
	close(ASSEMBLY);
      }


      ## Convert pattern-assembly result into PSSM
      $pssm_file = "$TMP/$tmp_file_name.pssm";
      $pssm_command = "$SCRIPTS/convert-matrix -v 1 ";
      $pssm_command .= " -in_format assembly -out_format patser";
      $pssm_command .= " -return counts,parameters";
      $pssm_command .= " -i $assembly_file";
      $pssm_command .= " -o $pssm_file";
      
      unless ($ENV{RSA_ERROR}) {
	print "<H2>Position-specific scoring matrices (PSSM)</H2>\n";
	system "$pssm_command";
	open PSSM, $pssm_file;
	print "<PRE>\n";
	while (<PSSM>) {
	  s|$RSA/||g;
	  print;
	}
	print "</PRE>\n";
	close(PSSM);
      }


    }
    
    &PipingForm();
    print '<HR SIZE=3>';
  
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command", $query->param('user_email'), $tmp_file_name);
} else {
    &EmailTheResult("$command", $query->param('user_email'), $tmp_file_name);
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
<TABLE>
<TR>

<TD>
<H3>Next step</H3>
</TD>
<TD>
<FORM METHOD="POST" ACTION="dna-pattern_form.cgi">
<INPUT type="hidden" NAME="title" VALUE="$title">
<INPUT type="hidden" NAME="pattern_file" VALUE="$result_file">
<INPUT type="hidden" NAME="sequence_file" VALUE="$sequence_file">
<INPUT type="hidden" NAME="sequence_format" VALUE="$sequence_format">
<INPUT type="submit" value="pattern matching (dna-pattern)">
</FORM>
</TD>

</TD>
<TD>
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



