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
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

#### TEMPORARY

$oligo_analysis_command = "$SCRIPTS/oligo-analysis";
$convert_seq_command = "$SCRIPTS/convert-seq";
$purge_sequence_command = "$SCRIPTS/purge-sequence";
$tmp_file_name = sprintf "oligo-analysis.%s", &AlphaDate;

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
    $command= "$purge_sequence_command -i $sequence_file -format $sequence_format -o ${sequence_file}.purged;  $oligo_analysis_command -i ${sequence_file}.purged  ";
} else {
    $command= "$oligo_analysis_command -i $sequence_file  ";
}


### fields to return
$return_fields = "";

if ($query->param('return') eq "table") {
    $parameters .= " -return occ -table"; 
} elsif ($query->param('return') eq "distrib") {
    $parameters .= " -return occ -distrib"; 
} else {
    
    ### occurrences
    if ($query->param('occ')) {
	$return_fields .= "occ,";

	### Lower threshold on occurrences
	if (&IsReal($query->param('lth_occ'))) {
	    $parameters .= " -lth occ ".$query->param('lth_occ');
	}

	### Upper threshold on occurrences
	if (&IsReal($query->param('uth_occ'))) {
	    $parameters .= " -uth occ ".$query->param('uth_occ');
	}
    } 
    
    ### frequencies
    if ($query->param('freq')) {
	$return_fields .= "freq,";

	### Lower threshold on frequencies
	if (&IsReal($query->param('lth_observed_freq'))) {
	    $parameters .= " -lth observed_freq ".$query->param('lth_observed_freq');
	}

	### Upper threshold on frequencies
	if (&IsReal($query->param('uth_observed_freq'))) {
	    $parameters .= " -uth observed_freq ".$query->param('uth_observed_freq');
	}
    } 
    
    ### matching sequences
    if ($query->param('mseq')) {
	$return_fields .= "mseq,";
	### threshold on matching sequences
	if ($query->param('lth_mseq') =~ /^\d+$/) {
	    $parameters .= " -thms ".$query->param('lth_mseq');
	}  
    } 
    
    ### observed/expected ratio
    if ($query->param('ratio')) {
	$return_fields .= "ratio,";
    } 
    
    ### rank
    if ($query->param('rank')) {
	$return_fields .= "rank,";
    } 
    
    ### z-score
    if ($query->param('zscore')) {
	$return_fields .= "zscore,";
    } 
    
    ### binomial probabilities
    if ($query->param('proba')) {
	$return_fields .= "proba,";

	### Lower threshold on probabilities
	if ($query->param('lth_occ_pro') =~ /^[\d\.\-+e]+$/i) {
	    $parameters .= " -lth occ_pro ".$query->param('lth_occ_pro');
	}

	### Upper threshold on probabilities
	if ($query->param('uth_occ_pro') =~ /^[\d\.\-+e]+$/i) {
	    $parameters .= " -uth occ_pro ".$query->param('uth_occ_pro');
	}

	### Lower threshold on significance
	if ($query->param('lth_occ_sig') =~ /^-{0,1}[\d\.\-+e]+$/i) {
	    $parameters .= " -lth occ_sig ".$query->param('lth_occ_sig');
	}

	### Upper threshold on significance
	if ($query->param('uth_occ_sig') =~ /^-{0,1}[\d\.\-+e]+$/i) {
	    $parameters .= " -uth occ_sig ".$query->param('uth_occ_sig');
	}
    } 
    
    ### positions
    if ($query->param('pos')) {
	$return_fields .= "pos,";
    } 
    
    $return_fields =~ s/,$//;
    

    if ($return_fields eq "") {
	&cgiError("You should select at least one option in the \"Return\" box.");
    } else {
	$parameters .= " -return $return_fields";
    }
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
  $freq_option = " -a input";
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
	$pattern_assembly_command = "$SCRIPTS/pattern-assembly -v 1 -subst 1";
	if ($query->param('strand') =~ /single/) {
	    $pattern_assembly_command .= " -1str";
	} else {
	    $pattern_assembly_command .= " -2str";
	}
	
	unless ($ENV{RSA_ERROR}) {
	    print "<H2>Pattern assembly</H2>\n";
	    open CLUSTERS, "$pattern_assembly_command -i $result_file |";
	    print "<PRE>\n";
	    while (<CLUSTERS>) {
		s|$RSA/||g;
		print;
	    }
	    print "</PRE>\n";
	    close(CLUSTERS);
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



